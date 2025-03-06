/*
  IOBUFF is a ld_preloaded library that intercepts and buffers sequential read and 
    write to file operations. 

    Copyright (C) 2013  <Jean-Francois St-Pierre> jf.stpierre@calculquebec.ca

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

 */

#include <dlfcn.h>
#include <errno.h>
#include <fcntl.h>
#include <pthread.h>
#include <map>
#include <list>
#include <algorithm>
#include <string>
#include <string.h>


#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stddef.h> 
#include <paths.h> 
#include <errno.h>
//#include <alloca.h>
#ifndef __MAX_ALLOCA_CUTOFF
#define __MAX_ALLOCA_CUTOFF   65536
#endif

#include <sys/resource.h>
#include <sys/socket.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>
#include <asm/unistd.h>

#include <vector>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <sys/un.h>

#include <sys/ipc.h>
#include <sys/shm.h>
#include <sys/sem.h>

#include "iobuff_glob_vars.h"
#include "iobuff_debug_calls.h"
#include "iobuff_debug_fork.h"
#include "iobuff_mpibuff_interface.h"

using namespace std;


////////////////////////////////////////////////////////////////////////////////////////////////////

//  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  // bfile
#ifdef MPIBUFF_SUPPORTED
  bfile::bfile(int fd, int flags_in=0, char *cname=NULL, int cnlength=0, bool fileIsBuffered=true, 
               bool fileIsMPIBuffered=false)   {
#else 
  bfile::bfile(int fd, int flags_in=0, char *cname=NULL, int cnlength=0, bool fileIsBuffered=true) {
#endif
    flags = flags_in;
    isBuffered = fileIsBuffered;
    pthread_mutex_lock(&mutex2);
    buffUID = globBuffUID++;
    pthread_mutex_unlock(&mutex2);
    if (isBuffered)
      buffer= new char[BSIZE];
    else
      buffer = NULL;
    buffer_ref = 1;
    file_begin_buff=0;
    if(cnlength != 0) {
#ifdef IOBUFF_DEBUG_BUFFER
      pthread_mutex_lock(&mutex2);
      ++nb_buff;
      if(nb_buff2++ % 500==0){
        fprintf(debugOut, "IOBUFF: nb_buff %d %d\n", nb_buff, nb_buff2);
      }
      pthread_mutex_unlock(&mutex2);
#endif
      canonical_name = new char[cnlength+1];
      strncpy(canonical_name, cname, cnlength);
      canonical_name[cnlength] = '\0';
    }
    else {
      canonical_name = NULL;
    }
    file_end_buff=0;
    buff_pos=0;
    buff_end=0;
    buff_need_flush=false;
    buff_need_fill=true;

    //lecture de la taille du fichier
    struct stat filestat;
    fstat(fd, &filestat);
    file_end_pos = filestat.st_size;
#ifdef MPIBUFF_SUPPORTED
    isMPIBuffered = fileIsMPIBuffered;
    //printf("IOBUFF: %d %s\n", (int) isMPIBuffered, canonical_name);
    if( isMPIBuffered ){
      mpibuff_pos=0;
      fBufId = _open_mpibuff_file(canonical_name, cnlength);
    }
#endif //MPIBUFF_SUPPORTED
          
  } 
//  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  // bfile
  bfile::~bfile() {
    pid_t lp = getpid();
#ifdef IOBUFF_DEBUG
    if (canonical_name != NULL)
      fprintf(debugOut,"IOBUFF: %d delete bfile %p %s\n", lp, this, canonical_name);
    else
      fprintf(debugOut,"IOBUFF: %d delete bfile %p\n", lp, this);
#endif
    if(isBuffered && buffer != NULL)
      delete[] buffer;
#ifdef IOBUFF_DEBUG_BUFFER
      pthread_mutex_lock(&mutex2);
      --nb_buff;
      if(nb_buff2++ % 500==0){
        fprintf(debugOut, "IOBUFF: nb_buff %d\n", nb_buff);
      }
      pthread_mutex_unlock(&mutex2);
#endif
    if (canonical_name != NULL)
      delete[] canonical_name;
#ifdef IOBUFF_DEBUG
  fprintf(debugOut,"IOBUFF: Exit destructor\n");
#endif
  }
//  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  // // //

// Map des fichiers ouverts (nom cannonique) et des bfile associé
struct cmp_str {
 bool operator()(char const *a, char const *b) {
   return strcmp(a, b) < 0;
 }
};

struct StrCompare : public std::binary_function<const char*, const char*, bool> {
public:
    bool operator() (const char* str1, const char* str2) const
    { return strcmp(str1, str2) < 0; }
};

map <std::string, list<bfile *> > fileMap;


//initialisation des pointeurs sur les fonctions de la libc 
volatile bool g_init_completed = false;
struct Init_lib{
//  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  Init_lib
  Init_lib () {
    g_supported_exec = true;  //Default behavior is to catch all executables
    lp = getpid();
#ifdef IOBUFF_DEBUG
    printf("IOBUFF:%d Init_lib\n",lp);
#endif
    // Check if any executable should be caught by IOBUFF 
    const char* c_env_iobuff_exec_list = getenv("IOBUFF_EXEC_LIST");
    if( c_env_iobuff_exec_list != NULL ) {                
      //Catch all executables by default
      const std::string env_iobuff_exec_list( c_env_iobuff_exec_list );
      const char delimiter = ':';
      size_t previous = 0;
      size_t index = env_iobuff_exec_list.find( delimiter );   
      string this_exec_name = string( canonicalize_file_name("/proc/self/exe") ); 
      g_supported_exec = false;   
      list<string> exec_name_list; 
      char *tmp_c_exec_name;
      string tmp_exec_name; 
      while( index != string::npos ){ 
        tmp_exec_name= env_iobuff_exec_list.substr(previous, index-previous); 
        if(tmp_exec_name.size() > 0 ){
          exec_name_list.push_back( tmp_exec_name );
          previous=index+1;
          index = env_iobuff_exec_list.find( delimiter, previous );
        }
      }
      tmp_exec_name = env_iobuff_exec_list.substr(previous); 
      if( tmp_exec_name.size() > 0 )
        exec_name_list.push_back( tmp_exec_name );
      list<string>::iterator ls_it;
      for(ls_it =exec_name_list.begin(); ls_it != exec_name_list.end(); ls_it++){
        string exec_name = *ls_it;
        string exec_tolower = exec_name;
        std::transform(exec_tolower.begin(), exec_tolower.end(), exec_tolower.begin(), ::tolower); 
        if( string("all").compare(exec_tolower) ==0 ){
          g_supported_exec=true; 
        }
        else if( exec_name.compare(0, 1, string("!"), 0,1) ==0) {
          tmp_c_exec_name = canonicalize_file_name( exec_name.c_str() +1 );
          if (tmp_c_exec_name != NULL) {
            string canon_exec_name = string( tmp_c_exec_name ) ; 
            if( canon_exec_name.compare( this_exec_name )  ==0) {
              // Negation of a path (!/path/exec) removes it from the list, overwriting any other settings
              g_supported_exec = false;
              break;
            }
          }
        }
        else {
          tmp_c_exec_name = canonicalize_file_name( exec_name.c_str() );
          if (tmp_c_exec_name != NULL) {
            string canon_exec_name = string( tmp_c_exec_name ) ; 
            if( canon_exec_name.compare( this_exec_name ) ==0) 
              g_supported_exec=true; 
          }
        } 
        previous=index+1;
        index = env_iobuff_exec_list.find( delimiter, previous );
      }
    }
    real_open =  (int     (*)(const char*, int, mode_t)) dlsym(RTLD_NEXT, "open");
    real_open64 =(int     (*)(const char*, int, mode_t)) dlsym(RTLD_NEXT, "open64");
    real_fopen =  (FILE*  (*)(const char *, const char *)) dlsym(RTLD_NEXT, "fopen"); 
    real_fdopen =  (FILE* (*)(int, const char *))        dlsym(RTLD_NEXT, "fdopen"); 
    real_freopen = (FILE* (*)(const char *, const char *, FILE *)) dlsym(RTLD_NEXT, "freopen"); 
    real_write = (ssize_t (*)(int,const void*,size_t))   dlsym(RTLD_NEXT, "write");
    real_pwrite = (ssize_t (*)(int,const void*,size_t,off_t)) dlsym(RTLD_NEXT, "pwrite");
    real_pwrite64 = (ssize_t (*)(int,const void*,size_t,off64_t))dlsym(RTLD_NEXT, "pwrite64");
    real_pwritev  = (ssize_t (*)(int,const struct iovec *,int,off_t))dlsym(RTLD_NEXT, "pwritev");
    real_writev   = (ssize_t (*)(int,const struct iovec *,int))dlsym(RTLD_NEXT, "writev");
    real_read =  (ssize_t (*)(int,void*,size_t))         dlsym(RTLD_NEXT, "read");
    real_pread =  (ssize_t (*)(int,void*,size_t,off_t))  dlsym(RTLD_NEXT, "pread");
    real_pread64 =(ssize_t (*)(int,void*,size_t,off64_t))dlsym(RTLD_NEXT, "pread64");
    real_preadv   = (ssize_t (*)(int,const struct iovec *,int,off_t))dlsym(RTLD_NEXT, "preadv");
    real_readv   = (ssize_t (*)(int,const struct iovec *,int))dlsym(RTLD_NEXT, "readv");
    real_lseek = (off_t   (*)(int,off_t,int))            dlsym(RTLD_NEXT, "lseek");
    real_llseek =(loff_t  (*)(int,loff_t,int))           dlsym(RTLD_NEXT, "llseek");
    real_lseek64=(off64_t (*)(int,off64_t,int))          dlsym(RTLD_NEXT, "lseek64");
    real_close = (int     (*)(int))                      dlsym(RTLD_NEXT, "close");
    real_fclose = (int    (*)(FILE*))                    dlsym(RTLD_NEXT, "fclose");
    real_dup   = (int     (*)(int))                      dlsym(RTLD_NEXT, "dup");
    real_dup2  = (int     (*)(int,int))                  dlsym(RTLD_NEXT, "dup2");
    real_dup3  = (int     (*)(int,int,int))              dlsym(RTLD_NEXT, "dup3");
    real_fcntl = (int     (*)(int,int,long))             dlsym(RTLD_NEXT, "fcntl");
    real_openat= (int     (*)(int,const char*,int,mode_t))dlsym(RTLD_NEXT, "openat");
    real_socket= (int     (*)(int,int,int))              dlsym(RTLD_NEXT, "socket");
    _init_fork_func();
    if( !g_supported_exec ) {
      //printf("IOBUFF:%d Init_lib exec not supported\n",lp); 
      g_init_completed = true;
      return;
    }
    //else
    //  printf("IOBUFF:%d Init_lib exec IS SUPPORTED\n",lp); 

    char *bsize = getenv("IOBUFF_BUFFER_SIZE");
    if(bsize!=NULL)
      BSIZE=atoi(bsize);
    
#ifdef MPIBUFF_SUPPORTED
    mpiSocketFD = _init_mpibuff_support();
#endif //MPIBUFF_SUPPORTED

    struct rlimit limit;
    getrlimit(RLIMIT_NOFILE, &limit);
    NB_FILES = limit.rlim_cur;

    vectFiles = new bfile*[NB_FILES];
    for (int i=0; i<NB_FILES; i++)
      vectFiles[i] = NULL;
      
    g_init_completed = true;
  } //Init_lib
//  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  
  ~Init_lib() {
    if( !g_supported_exec ) 
      return;            //Nothing to do here
      
    for (int i=0; i<NB_FILES; i++)
      if (vectFiles[i]){
        close (i);
      }
    delete [] vectFiles;
#ifdef MPIBUFF_SUPPORTED
    _close_mpibuff_support(mpiSocketFD);
#endif
  } //~Init_lib
} *objInit; 

////////////////////////////////////////////////////////////////////////////////////////////// _init 
pthread_mutex_t mutexInit = PTHREAD_ADAPTIVE_MUTEX_INITIALIZER_NP;
void _init(void){ 
  pthread_mutex_lock( &mutexInit );
    if( !g_init_completed )  //Recheck under mutex
      objInit = new struct Init_lib;
  pthread_mutex_unlock( &mutexInit );
}//_init 
///////////////////////////////////////////////////////////////////////////////////// spin_wait_init
void spin_wait_init() { 
  while( !g_init_completed ) {
    _init();
  }
}
////////////////////////////////////////////////////////////////////////////////////////// get_bfile
bfile * get_bfile(int fd){
  return vectFiles[fd];
}//get_bfile

/////////////////////////////////////////////////////////////////////////////////////// delete_bfile
void delete_bfile(int fd) {
  _init_debug_delete_bfile(fd);
  if (fd>=0) {
    bfile * lbfile = vectFiles[fd];
    if (lbfile) {
      if (lbfile->buffer_ref == 1){
        delete lbfile;
      }
      else if (lbfile->buffer_ref < 1){
        fprintf(debugOut, "IOBUFF: Buffer ref error for lbfile %p named %s", 
                          lbfile, lbfile->canonical_name);
        exit(1);
      }else{
        lbfile->buffer_ref--;
      }
    }
    vectFiles[fd] = NULL;
  }
  _exit_debug_delete_bfile();
}//delete_bfile
/////////////////////////////////////////////////////////////////////////////////////// flush_buffer
bool flush_buffer(int fd, bfile * lbfile) {
  if(lbfile->buff_need_flush && lbfile->isBuffered) {      
    real_lseek(fd, lbfile->file_begin_buff, SEEK_SET);
    int wrote = real_write(fd, lbfile->buffer, lbfile->buff_end);
    lbfile->buff_need_flush=false;
    //lbfile->file_begin_buff += wrote;
    lbfile->file_begin_buff += lbfile->buff_pos; //New : assume buff_pos != wrote
    real_lseek(fd, lbfile->file_begin_buff, SEEK_SET);
    lbfile->file_end_buff = lbfile->file_begin_buff;
    lbfile->buff_pos=0;
    return true;
  }
  return false;
}//flush_buffer

//////////////////////////////////////////////////////////////////////////////////////// fill_buffer
void fill_buffer(int fd, bfile * lbfile) {
  if (lbfile->buff_need_fill) {
    real_lseek(fd, lbfile->file_begin_buff + lbfile->buff_end, SEEK_SET);
    ssize_t ret ;
#ifdef MPIBUFF_SUPPORTED
    if( lbfile->isMPIBuffered ) 
      ret = mpibuff_read(fd, lbfile->buffer + lbfile->buff_end, BSIZE - lbfile->buff_end);
    else 
#endif
      ret = real_read(fd, lbfile->buffer + lbfile->buff_end, BSIZE - lbfile->buff_end);
    if(ret < 0) {
      fprintf(debugOut,"IOBUFF : READ ERROR 1!!! fd=%d\n",fd);
      //printf("IOBUFF : READ ERROR 1!!! fd=%d\n",fd);
      exit(1);
    }
    lbfile->buff_end += ret;
    lbfile->file_end_buff = lbfile->file_begin_buff + lbfile->buff_end;
    lbfile->buff_need_fill = false;
  }
}//fill_buffer

//////////////////////////////////////////////////////////////////////////////// update_file_end_pos
void update_file_end_pos(bfile * lbfile) {
  if (lbfile->file_end_buff > lbfile->file_end_pos)
    lbfile->file_end_pos = lbfile->file_end_buff;
}//update_file_end_pos

/////////////////////////////////////////////////////////////////////////////////////// exclude_path
bool exclude_path( char** cname ){
  for(int i = 0; i< NB_EXCLUDE; i++) {
    if( strncmp(*cname, EXCLUDE[i], strlen(EXCLUDE[i])) ==0) { 
      free(*cname);
      return true;
    }
  }
  return false;
}//exclude_path

//////////////////////////////////////////////////////////////////////////// get_canonical_file_name
void get_canonical_file_name(const char* pathname, char** cname, ssize_t* cnlen) {
  //if ((cnlen = readlink(pathname, cname, sizeof(cname)-1)) != -1) {
  if ((*cname = canonicalize_file_name(pathname)) != NULL) {
    //printf("IOBUFF: Cname found ,%s,\n",cname);
    *cnlen=strlen(*cname); 
  }
  else{
    fprintf(debugOut, pathname);
    perror("Cname not-found: "); //should not happen 
    *cnlen=strlen(pathname);
    *cname = (char*) malloc(((*cnlen)+1)*sizeof(char));
    strncpy(*cname, pathname, *cnlen);
    (*cname)[*cnlen] = '\0'; 
  }
}//get_canonical_file_name

/////////////////////////////////////////////////////////////////////////////////////// enlist_bfile
void enlist_bfile(bfile * lbfile){ // Keep track of all bfiles opened with one canonical name
  list<bfile *> *bfileList;
  bfileList = &fileMap[std::string(lbfile->canonical_name)];  //List associated with canonical name
  _init_debug_enlist_bfile(bfileList);

#ifdef IOBUFF_DISALLOW_MIXED_OPEN_TYPES
  //Can't mix buffered and unbuffered files
  if (bfileList->size() > 0 ) {
    list<bfile *>::iterator it;
    for (it=bfileList->begin(); it != bfileList->end(); it++){ 
      if ((*it)->isBuffered != lbfile->isBuffered){
        fprintf(debugOut,"IOBUFF: Can't open the same files using fopen and open : %s \n",
                         (*it)->canonical_name);  
        pthread_mutex_unlock(&mutex);
        exit(1);
      }
    }
  }// if (bfileList->size() > 0 )
#endif
  bfileList->push_back(lbfile); 
   
  //Un fichier peut être ouvert en O_RDONLY plusieurs fois, mais pas en O_WRONLY ou O_RDWR
  if (bfileList->size() != 1) {
    list<bfile *>::iterator it;
    for (it=bfileList->begin(); it != bfileList->end(); it++){ 
      if ((*it)->flags & O_WRONLY || (*it)->flags & O_RDWR){
        fprintf(debugOut,"IOBUFF: Can't open the same writable files concurently : %s \n",
                         (*it)->canonical_name); 
        pthread_mutex_unlock(&mutex);
        exit(1);
      }
    }
  }// if (bfileList->size() != 1)
  _exit_debug_enlist_bfile(bfileList);
}//enlist_bfile

/////////////////////////////////////////////////////////////////////////////////////// delist_bfile
void delist_bfile(int fd, bfile* lbfile){ 
  _init_debug_delist_bfile();
  if(lbfile->buffer_ref == 1) { //Last duplicate of the file (dup)
    list<bfile *> *bfileList = &fileMap[std::string(lbfile->canonical_name)];
    _debug_delist_bfile(bfileList);  

    if(bfileList->size() ==1) {  //No other open copy of the file
      int tmp_cnlength = strlen(lbfile->canonical_name);
      char *tmp_cname = new char[tmp_cnlength+1];
      strncpy(tmp_cname, lbfile->canonical_name, tmp_cnlength);
      tmp_cname[tmp_cnlength] = '\0'; 
      bfileList->remove(lbfile); 
      fileMap.erase(std::string(tmp_cname));
      delete [] tmp_cname;
    } 
    else { //File is still open on a different non-dup fd. Only remove this fd. 
      bfileList->remove(lbfile);
    }
    delete_bfile(fd); 
  }//if(lbfile->buffer_ref == 1)
  else {// Some duplicate of the file exist (dup)
    delete_bfile(fd);
  }
  _exit_debug_delist_bfile();
}//delist_bfile

/////////////////////////////////////////////////////////////////////////////////////// mode_to_flag
int mode_to_flag(const char* mode){
  int flags = 0;
  bool rflag, eflag, bflag, wflag, xflag, aflag, plusflag;
  rflag = eflag = bflag = wflag = xflag = aflag = plusflag = false;
  for(int i=0; i< strlen(mode); i++){
    switch(mode[i]){
      case 'r': case 'R': rflag=true; break;
      case 'e': case 'E': eflag=true; break;
      case 'b': case 'B': bflag=true; break;
      case 'w': case 'W': wflag=true; break;
      case 'x': case 'X': xflag=true; break;
      case 'a': case 'A': aflag=true; break;
      case '+': plusflag=true; break;
      default : fprintf(debugOut,"IOBUFF: Unknow f(re)open mode : \"%s\" \n", mode); break;
      }
    }
        
  if (rflag && !plusflag)
    flags = O_RDONLY;
  else if((wflag || aflag) && !plusflag)
    flags = O_WRONLY;
  else if((rflag || wflag || aflag) && plusflag)
    flags = O_RDWR;
  else {
    fprintf(debugOut,"IOBUFF: Unknow f(re)open mode : \"%s\" \n", mode);
    exit(1);
  }
  return flags;
}//mode_to_flag
///////////////////////////////////////////////////////////////////////////////////// debuffer_bfile
void debuffer_bfile(int fd, bfile* lbfile){
  if(!flush_buffer(fd, lbfile) ){                            
    real_lseek(fd, lbfile->file_begin_buff + lbfile->buff_pos, SEEK_SET);
    //real_lseek(fd, lbfile->file_begin_buff, SEEK_SET);
  }
  lbfile->isBuffered=false;
  if(lbfile->buffer != NULL){
    delete[] lbfile->buffer;
    lbfile->buffer = NULL;
  }
}//debuffer_bfile
//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\//\\//\\//\\/
extern "C" { //\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\

///////////////////////////////////////////////////////////////////////////////////////////// open64
  int open64(const char *pathname, int flags, ...) {
    _init_debug_open64();
    spin_wait_init();
    mode_t mode;
    if (flags & O_CREAT) {
      va_list arg;
      va_start(arg, flags);
      mode = va_arg(arg, mode_t);
      va_end(arg);
    }
    if( !g_supported_exec )
      return real_open64(pathname,flags,mode);
    return open(pathname,flags,mode);
  }//open64
  
/////////////////////////////////////////////////////////////////////////////////////////////// open
  int open(const char *pathname, int flags, ...) {
    _init_debug_open();
    spin_wait_init();
    mode_t mode;
    if (flags & O_CREAT) {
      va_list arg;
      va_start(arg, flags);
      mode = va_arg(arg, mode_t);
      va_end(arg);
    }
    if( !g_supported_exec)  //Pass-through
      return real_open(pathname,flags,mode);
      
    int p,fd;
    char *cname;//nom canonique  
    ssize_t cnlen=0; 

    bfile * lbfile;

    fd=real_open(pathname,flags,mode);
    _debug_open(fd, pathname, flags);
    pthread_mutex_lock(&mutex);  
    if (fd != -1 && vectFiles[fd] != NULL){ 
      fprintf(debugOut,"IOBUFF: open ERROR : vectFiles[fd] not NULL : %d %p\n", fd, vectFiles[fd]);
      fprintf(debugOut,"IOBUFF: new name %s \n", canonicalize_file_name(pathname));
      fprintf(debugOut,"IOBUFF: old name %s vectFiles[fd] \n", vectFiles[fd]->canonical_name);
      
      list<bfile *> *bfileList;
      bfileList = &fileMap[  std::string(vectFiles[fd]->canonical_name)  ];
      fprintf(debugOut,"IOBUFF: bflist %p %d vectFiles[fd] %p canonical_name %p\n", bfileList, 
                       bfileList->size(), vectFiles[fd], vectFiles[fd]->canonical_name);
      
      pthread_mutex_unlock(&mutex);  
      exit(1);
    }
    else if (fd != -1 && vectFiles[fd] == NULL) {
      //Extrait le nom canonique du fichier
      get_canonical_file_name(pathname, &cname, &cnlen);   
      if (exclude_path(&cname)){
        pthread_mutex_unlock(&mutex);  
        _exit_debug_open(fd);
        return fd; 
      }
#ifdef MPIBUFF_SUPPORTED
      if(std::find(mpibufFiles.begin(), mpibufFiles.end(), std::string(cname))!=mpibufFiles.end()){
        if ( (flags & O_ACCMODE) != O_RDONLY){
          fprintf(debugOut, "MPIBUF error : Only O_RDONLY is supported\n"); 
          exit(1);
        }
        //vectFiles[fd] = new bfile (fd, flags, cname, cnlen, true, true ); 
        vectFiles[fd] = new bfile (fd, flags, cname, cnlen, false, true ); 
      }
      else
#endif  
        vectFiles[fd] = new bfile (fd, flags, cname, cnlen);
      _debug_open2(fd, pathname, cname, cnlen);
      free(cname);

      //Ajout du fichier à la map des noms de fichier
      enlist_bfile(vectFiles[fd]);  
    } 
    pthread_mutex_unlock(&mutex);
    _exit_debug_open(fd);
    return fd;
  }//open

////////////////////////////////////////////////////////////////////////////////////////////// pread
  ssize_t pread(int fd, void *buf, size_t count, off_t offset){
     _init_debug_pread(fd, count, offset);
    spin_wait_init();
    if( !g_supported_exec )  //Pass-through
      return real_pread(fd, buf, count, offset);
    bfile * lbfile = get_bfile(fd);

    if (!lbfile)
      return real_pread(fd, buf, count, offset);
#ifdef MPIBUFF_SUPPORTED
    else if ( lbfile->isMPIBuffered )  {
      return mpibuff_pread(fd, buf, count, offset);
    }
#endif
    else if ( !lbfile->isBuffered )  
      return real_pread(fd, buf, count, offset);
#ifdef IOBUFF_PREAD_PWRITE_SUPPORTED
    else if(lbfile->flags & O_WRONLY || lbfile->flags & O_RDWR){  
      pthread_mutex_lock(&mutex);                                 //Not supporting pread in writable
        debuffer_bfile(fd, lbfile);                               //file. de-buffering 
      pthread_mutex_unlock(&mutex);
    }
    _exit_debug_pread();
    return real_pread(fd, buf, count, offset);
#else
    fprintf(debugOut, "IOBUFF: pread not supported on IOBUFF files\n");
    exit(1);
#endif
  }//pread
  
//////////////////////////////////////////////////////////////////////////////////////////// pread64
  ssize_t pread64(int fd, void *buf, size_t count, off64_t offset){
     _init_debug_pread(fd, count, offset);
    spin_wait_init();
    if( !g_supported_exec )  //Pass-through
      return real_pread64(fd, buf, count, offset);
    bfile * lbfile = get_bfile(fd);

    if (!lbfile)
      return real_pread64(fd, buf, count, offset);
#ifdef MPIBUFF_SUPPORTED
    else if ( lbfile->isMPIBuffered )  {
      return mpibuff_pread(fd, buf, count, offset);
    }
#endif
    else if ( !lbfile->isBuffered )  
      return real_pread64(fd, buf, count, offset);
#ifdef IOBUFF_PREAD_PWRITE_SUPPORTED
    else if(lbfile->flags & O_WRONLY || lbfile->flags & O_RDWR){  
      pthread_mutex_lock(&mutex);                               //Not supporting pread64 in writable
        debuffer_bfile(fd, lbfile);                             //file. de-buffering 
      pthread_mutex_unlock(&mutex);
    }
    _exit_debug_pread64();
    return real_pread64(fd, buf, count, offset);
#else
      
    fprintf(debugOut, "IOBUFF: pread64 not supported on IOBUFF files\n");
    exit(1);
#endif
  }//pread64

///////////////////////////////////////////////////////////////////////////////////////////// preadv
  ssize_t preadv(int fd, const struct iovec *iov, int iovcnt, off_t offset){
    spin_wait_init();
    if( !g_supported_exec )  //Pass-through
      return real_preadv(fd, iov, iovcnt, offset);
    
    bfile * lbfile = get_bfile(fd);
    if (!lbfile)
      return real_preadv(fd, iov, iovcnt, offset);
      
#ifdef MPIBUFF_SUPPORTED
    if ( !lbfile->isMPIBuffered && !lbfile->isBuffered )  {
#else
    if ( !lbfile->isBuffered )  {
#endif 
      return real_preadv(fd, iov, iovcnt, offset);
    }
#if defined(IOBUFF_PREAD_PWRITE_SUPPORTED) && defined(IOBUFF_READV_WRITEV_SUPPORTED)
    if(lbfile->flags & O_WRONLY || lbfile->flags & O_RDWR){
      pthread_mutex_lock(&mutex);                                //Not supporting preadv in writable
        debuffer_bfile(fd, lbfile);                              //file. de-buffering 
      pthread_mutex_unlock(&mutex);
      return real_preadv(fd, iov, iovcnt, offset);
    }
    pthread_mutex_lock(&mutex);                      //preadv is atomic !! watch for mutex in pread 
    int i;
    off_t curOffset = offset;
    ssize_t ret=0, curRet=0;
    for (i = 0; i < iovcnt; i++){
      curRet = pread(fd, iov[i].iov_base, iov[i].iov_len, curOffset);
      if(curRet == iov[i].iov_len){ //Normal case
        ret += curRet;
        curOffset += curRet;
      }
      else if(curRet< 0){             //Error case 
        pthread_mutex_unlock(&mutex); 
        return curRet;
      }
      else {                          //Incomplete case
        ret += curRet;
        pthread_mutex_unlock(&mutex); 
        return ret;
      }
    }
    pthread_mutex_unlock(&mutex); 
    return ret;
#else
    fprintf(debugOut, "IOBUFF: preadv not supported on IOBUFF files\n");
    exit(1);
#endif 
  }//preadv
////////////////////////////////////////////////////////////////////////////////////// read_no_mutex
  ssize_t read_no_mutex(int fd, void *buf, size_t count){    
    size_t inextra,inbuffer;                        // Needed for func read and readv that
    size_t offset = 0;                              // provide their own mutex
    size_t red = 0;

    bfile * lbfile = get_bfile(fd);

    if (!lbfile){ 
      _exit_debug_read();
      return real_read (fd, buf, count);
    } 
#ifdef MPIBUFF_SUPPORTED
    else if ( lbfile->isMPIBuffered && !lbfile->isBuffered ) { 
      ssize_t ret = mpibuff_read(fd, buf, count);
      _exit_debug_read();
      return ret;
    } 
#endif
    else if ( !lbfile->isBuffered ) { 
      _exit_debug_read();
      return real_read (fd, buf, count);   
    }

    //fait qu'une lecture pour les grande lecture
    if (count > BSIZE) {
      //se place a la position equivalente de lecture du buffer dans le fichier
#ifdef MPIBUFF_SUPPORTED
      if( lbfile->isMPIBuffered ) {
        mpibuff_lseek64(fd, lbfile->file_begin_buff + lbfile->buff_pos, SEEK_SET);
        red = mpibuff_read(fd, buf, count);
      }
      else  
#endif
     {
        real_lseek(fd, lbfile->file_begin_buff + lbfile->buff_pos, SEEK_SET);
        red = real_read(fd, buf, count);
      }
      lbfile->buff_need_fill=true;
      lbfile->file_end_buff = lbfile->file_begin_buff + lbfile->buff_pos + red;
      lbfile->file_begin_buff = lbfile->file_end_buff;
      lbfile->buff_end = 0;
      lbfile->buff_pos = 0;
      _exit_debug_read();
      return red;
    }
    fill_buffer (fd, lbfile);

    if(count <= (lbfile->buff_end - lbfile->buff_pos)) {
      inbuffer = count;
      inextra = 0;                                                     
    }
    else {                                                                                           
      inbuffer = (lbfile->buff_end - lbfile->buff_pos);
      inextra = count - inbuffer;
    }
      
    if( inbuffer>0 ) {
      memcpy((char*)buf, lbfile->buffer+lbfile->buff_pos, inbuffer);
      lbfile->buff_pos += inbuffer;
      offset += inbuffer;
      red += inbuffer;
    }

    if (inextra > 0) {
      //ecrit le buffer s'il a ete modifie par une ecriture
      if(!flush_buffer(fd, lbfile) ){ 
        //lecture du prochain buffer
#ifdef MPIBUFF_SUPPORTED
        if( lbfile->isMPIBuffered ) 
          mpibuff_lseek64(fd, lbfile->file_end_buff, SEEK_SET);
        else 
#endif
          real_lseek(fd, lbfile->file_end_buff, SEEK_SET);
      }
      ssize_t ret ;
#ifdef MPIBUFF_SUPPORTED
      if( lbfile->isMPIBuffered ) 
        ret = mpibuff_read(fd, lbfile->buffer, BSIZE);
      else 
#endif
        ret = real_read(fd, lbfile->buffer, BSIZE);
      if(ret < 0){
        fprintf(debugOut,"IOBUFF : READ ERROR 2!!! \n");
        //printf("IOBUFF : READ ERROR 2!!! \n");
        exit(1);
      }
      lbfile->buff_need_fill=false;
      lbfile->file_begin_buff = lbfile->file_end_buff;
      lbfile->file_end_buff += ret;
      lbfile->buff_end = ret;
  
      if (inextra > ret)
        inextra = ret;
      memcpy((char*)buf+offset, lbfile->buffer, inextra);
      lbfile->buff_pos = inextra;
      red += inextra;
    }
    _exit_debug_read();
    return red;
  }//read_no_mutex
  
/////////////////////////////////////////////////////////////////////////////////////////////// read
  ssize_t read(int fd, void *buf, size_t count){
    _init_debug_read(fd, buf, count);
    spin_wait_init();
    if( !g_supported_exec )  //Pass-through
      return real_read (fd, buf, count);

#ifdef IOBUFF_MUTEX_ON_ALL_OPS
    pthread_mutex_lock(&mutex);
#endif

    ssize_t ret = read_no_mutex(fd, buf, count);
    
#ifdef IOBUFF_MUTEX_ON_ALL_OPS
    pthread_mutex_unlock(&mutex);
#endif
    return ret;
  }//read
///////////////////////////////////////////////////////////////////////////////////////////// readv
  ssize_t readv(int fd, const struct iovec *iov, int iovcnt){
    spin_wait_init();
    if( !g_supported_exec )  //Pass-through
      return real_readv(fd, iov, iovcnt);
      
#ifdef IOBUFF_READV_WRITEV_SUPPORTED
    bfile * lbfile = get_bfile(fd);
    if (!lbfile || !(lbfile->isBuffered)){
      return real_readv(fd, iov, iovcnt);
    } 

    pthread_mutex_lock(&mutex);                      //readv is atomic !! need write_no_mutex 
    int i;
    ssize_t ret=0, curRet=0;
    for (i = 0; i < iovcnt; i++){
      curRet = read_no_mutex(fd, iov[i].iov_base, iov[i].iov_len);
      if(curRet == iov[i].iov_len){ //Normal case
        ret += curRet;
      }
      else if(curRet< 0){             //Error case 
        pthread_mutex_unlock(&mutex); 
        return curRet;
      }
      else {                          //Incomplete case
        ret += curRet;
        pthread_mutex_unlock(&mutex); 
        return ret;
      }
    }
    pthread_mutex_unlock(&mutex); 
    return ret;
#else 
    fprintf(debugOut, "IOBUFF: readv not supported on IOBUFF files\n");
    exit(1);
#endif //IOBUFF_READV_WRITEV_SUPPORTED

  }//readv
  
///////////////////////////////////////////////////////////////////////////////////////////// pwrite
  ssize_t pwrite(int fd, const void *buf, size_t count, off_t offset) {
     _init_debug_pwrite(fd, count, offset);
    spin_wait_init();
    if( !g_supported_exec )  //Pass-through
      return real_pwrite(fd, buf, count, offset);
    bfile * lbfile = get_bfile(fd);

    if (!lbfile || !(lbfile->isBuffered))  
      return real_pwrite(fd, buf, count, offset);
#ifdef IOBUFF_PREAD_PWRITE_SUPPORTED
    else {
      pthread_mutex_lock(&mutex);                                  //pwrite in buffered file not
        debuffer_bfile(fd, lbfile);                                //supported. de-buffering 
      pthread_mutex_unlock(&mutex);
    }
    return real_pwrite(fd, buf, count, offset);
#else

    fprintf(debugOut, "IOBUFF: pwrite not supported on IOBUFF files\n");
    exit(1);
#endif
  }//pwrite
  
/////////////////////////////////////////////////////////////////////////////////////////// pwrite64
  ssize_t pwrite64(int fd, const void *buf, size_t count, off64_t offset) {
     _init_debug_pwrite64(fd, count, offset);
    spin_wait_init();
    if( !g_supported_exec )  //Pass-through
      return real_pwrite64(fd, buf, count, offset);
      
    bfile * lbfile = get_bfile(fd);

    if (!lbfile || !(lbfile->isBuffered))  
      return real_pwrite64(fd, buf, count, offset);
#if defined(IOBUFF_PREAD_PWRITE_SUPPORTED) && defined(IOBUFF_READV_WRITEV_SUPPORTED)
    else {
      pthread_mutex_lock(&mutex);                                  //pwrite64 in buffered file not
        debuffer_bfile(fd, lbfile);                                //supported. de-buffering 
      pthread_mutex_unlock(&mutex);
    }
    return real_pwrite(fd, buf, count, offset);
#else

    fprintf(debugOut, "IOBUFF: pwrite64 not supported on IOBUFF files\n");
    exit(1);
#endif
  }//pwrite64
//////////////////////////////////////////////////////////////////////////////////////////// pwritev
  ssize_t pwritev(int fd, const struct iovec *iov, int iovcnt, off_t offset){
    spin_wait_init();
    if( !g_supported_exec )  //Pass-through
      return real_pwritev(fd, iov, iovcnt, offset);      
      
    bfile * lbfile = get_bfile(fd);
    if (!lbfile || !(lbfile->isBuffered))  
      return real_pwritev(fd, iov, iovcnt, offset);      
#if defined(IOBUFF_PREAD_PWRITE_SUPPORTED) && defined(IOBUFF_READV_WRITEV_SUPPORTED)
    else {
      pthread_mutex_lock(&mutex);                                  //pwritev in buffered file not
        debuffer_bfile(fd, lbfile);                                //supported. de-buffering 
      pthread_mutex_unlock(&mutex);
    }
    return real_pwritev(fd, iov, iovcnt, offset);      
#else
    fprintf(debugOut, "IOBUFF: pwritev not supported on IOBUFF files\n");
    exit(1);
#endif
  }// pwritev
  
///////////////////////////////////////////////////////////////////////////////////// write_no_mutex
  ssize_t write_no_mutex(int fd, const void *buf, size_t count) {
    size_t inbuffer;                      // Need for func. that already provide their
    size_t inextra;                       // own mutex (writev and write) 
    size_t offset=0;
    size_t wrote=0;
    bfile * lbfile = get_bfile(fd);

    //fait qu'une lecture pour les grandes lectures
    if (count > BSIZE) {
      //ecriture du buffer avant de faire la grosse ecriture
      if(!flush_buffer(fd, lbfile) ){
        real_lseek(fd, lbfile->file_begin_buff, SEEK_SET);
      }
      wrote = real_write(fd, buf, count);
      lbfile->buff_need_flush=0;
      lbfile->file_end_buff += wrote;
      update_file_end_pos(lbfile);
      lbfile->file_begin_buff = lbfile->file_end_buff;
      lbfile->buff_end = 0;
      lbfile->buff_pos = 0;
      return wrote;
    }

    inbuffer = BSIZE - lbfile->buff_pos;
    if( count<=inbuffer) {
      inbuffer = count;
      inextra = 0;
    }
    else {
      inextra = count - inbuffer;
    }

    if(inbuffer>0) {
      memcpy(lbfile->buffer+lbfile->buff_pos,(char*)buf,inbuffer);
      offset += inbuffer;
      lbfile->buff_pos += inbuffer;
      wrote += inbuffer;
  
      if(lbfile->buff_end<lbfile->buff_pos) {
        lbfile->buff_end = lbfile->buff_pos;
        lbfile->file_end_buff = lbfile->file_begin_buff + lbfile->buff_end;
        update_file_end_pos(lbfile);
      }
    }

    if( inextra>0 ) {
      //ecritue du buffer courant
      real_lseek(fd, lbfile->file_begin_buff, SEEK_SET);
      long ret = real_write( fd, lbfile->buffer, lbfile->buff_end);
      if(ret<=0) {
        return ret;
      }
      lbfile->file_begin_buff += ret;
      lbfile->file_end_buff = lbfile->file_begin_buff;
      
      lbfile->buff_need_fill = true;
      memcpy(lbfile->buffer,(char*)buf+offset,inextra);
      wrote+=inextra;
      lbfile->buff_pos = inextra;
      lbfile->buff_end = lbfile->buff_pos;
      //la fin du fichier peu augmenter
      lbfile->file_end_buff += ret;
      update_file_end_pos(lbfile);
    }

    if (wrote>0)
      lbfile->buff_need_flush = true;
    return wrote;
  }
  
////////////////////////////////////////////////////////////////////////////////////////////// write
  ssize_t write(int fd, const void *buf, size_t count) {
    spin_wait_init();
    if( !g_supported_exec )  //Pass-through
      return real_write (fd, buf, count);
      
    bfile * lbfile = get_bfile(fd);
    if (!lbfile || !(lbfile->isBuffered)){
      return real_write (fd, buf, count);
    }
#ifdef IOBUFF_MUTEX_ON_ALL_OPS
    pthread_mutex_lock(&mutex); 
#endif

    ssize_t ret = write_no_mutex(fd, buf, count);
    
#ifdef IOBUFF_MUTEX_ON_ALL_OPS
    pthread_mutex_unlock(&mutex);
#endif
    return ret;
  }//write
  
///////////////////////////////////////////////////////////////////////////////////////////// writev
  ssize_t writev(int fd, const struct iovec *iov, int iovcnt){
    spin_wait_init();
    if( !g_supported_exec )  //Pass-through
      return real_writev(fd, iov, iovcnt);
      
#ifdef IOBUFF_READV_WRITEV_SUPPORTED
    bfile * lbfile = get_bfile(fd);
    if (!lbfile || !(lbfile->isBuffered)){
      return real_writev(fd, iov, iovcnt);
    } 

    pthread_mutex_lock(&mutex);                      //writev is atomic !! need write_no_mutex 
    int i;
    ssize_t ret=0, curRet=0;
    for (i = 0; i < iovcnt; i++){
      curRet = write_no_mutex(fd, iov[i].iov_base, iov[i].iov_len);
      if(curRet == iov[i].iov_len){   //Normal case
        ret += curRet;
      }
      else if(curRet< 0){             //Error case 
        pthread_mutex_unlock(&mutex); 
        return curRet;
      }
      else {                          //Incomplete case
        ret += curRet;
        pthread_mutex_unlock(&mutex); 
        return ret;
      }
    }
    pthread_mutex_unlock(&mutex); 
    return ret;
#else 
    fprintf(debugOut, "IOBUFF: writev not supported on IOBUFF files\n");
    exit(1);
#endif //IOBUFF_READV_WRITEV_SUPPORTED

  }//writev
////////////////////////////////////////////////////////////////////////////////////////////// lseek
  off_t lseek(int fd, off_t offset, int whence) throw() {
    _init_debug_lseek(fd, offset, whence);
    spin_wait_init();
    if( !g_supported_exec )  //Pass-through
      return  real_lseek(fd,offset,whence); 
    return lseek64(fd,offset,whence); 
  }//lseek
  
  
///////////////////////////////////////////////////////////////////////////////////////////// llseek
  loff_t llseek(int fd, loff_t offset, int whence) throw() {
    spin_wait_init();
    if( !g_supported_exec )  //Pass-through
      return  real_llseek(fd, offset, whence) ;
    //return real_llseek(fd,offset,whence);
    fprintf(debugOut,"IOBUFF: llseek not implemented yet!!!! aborting.\n"); 
    exit(1);
  }//llseek
  
//////////////////////////////////////////////////////////////////////////////////////////// lseek64
  off64_t lseek64(int fd, off64_t offset, int whence) throw() {
    _init_debug_lseek64(fd, offset, whence);
    spin_wait_init();
    if( !g_supported_exec )  //Pass-through
      return  real_lseek64(fd,offset,whence); 

#ifdef IOBUFF_MUTEX_ON_ALL_OPS
      pthread_mutex_lock(&mutex);
#endif
    off64_t trans,off;
    bfile * lbfile = get_bfile(fd);
   
    if (!lbfile){ 
#ifdef IOBUFF_MUTEX_ON_ALL_OPS
      pthread_mutex_unlock(&mutex);
#endif
      return real_lseek64(fd,offset,whence); 
    } 
#ifdef MPIBUFF_SUPPORTED
    else if ( lbfile->isMPIBuffered ) { 
#ifdef IOBUFF_MUTEX_ON_ALL_OPS
      pthread_mutex_unlock(&mutex);
#endif
      return mpibuff_lseek64(fd,offset,whence); 
    } 
#endif
    else if ( !lbfile->isBuffered ) { 
#ifdef IOBUFF_MUTEX_ON_ALL_OPS
      pthread_mutex_unlock(&mutex);
#endif
      return real_lseek64(fd,offset,whence); 
    }
    
    switch (whence) {
    case SEEK_SET:
      trans = offset;
      break;
    case SEEK_CUR:
      trans = lbfile->file_begin_buff + lbfile->buff_pos + offset;
      if (trans < 0) {
        errno = EINVAL;
#ifdef IOBUFF_MUTEX_ON_ALL_OPS
      pthread_mutex_unlock(&mutex);
#endif
        return -1;
      }
      break;
    case SEEK_END:
      trans = lbfile->file_end_pos + offset;
      break;
    }
  
    //position est dans le fichier
    if ( trans < lbfile->file_end_pos ) {
      //avant ou apres le buffer
      if( trans < lbfile->file_begin_buff || trans >= lbfile->file_end_buff ) {
        flush_buffer(fd, lbfile);
  
        lbfile->file_begin_buff = lbfile->file_end_buff = trans;
        lbfile->buff_pos = 0;
        lbfile->buff_end = 0;
        lbfile->buff_need_fill=true;
      }
      else { //dans le buffer 
        off = trans - lbfile->file_begin_buff;
        lbfile->buff_pos=off;
      }
    }
    //on est apres la fin du fichier
    else {
      flush_buffer(fd, lbfile);
    
      lbfile->file_begin_buff = lbfile->file_end_buff = trans;
      lbfile->buff_pos = 0;
      lbfile->buff_end = 0;
      lbfile->buff_need_fill=false;
    }
    
#ifdef IOBUFF_MUTEX_ON_ALL_OPS
      pthread_mutex_unlock(&mutex);
#endif
    return trans;
  }//lseek64
  
//////////////////////////////////////////////////////////////////////////////////////////////// dup
  int dup(int oldfd) throw() {
    _init_debug_dup(oldfd);
    spin_wait_init();
    if( !g_supported_exec )  //Pass-through
      return  real_dup(oldfd);
      
    int fd;
    bfile * lbfile  = get_bfile(oldfd);
    fd=real_dup(oldfd);

    if (fd != -1) {
      if (lbfile){
        lbfile->buffer_ref++;
      }
      if (vectFiles[fd] != NULL){ 
        fprintf(debugOut,"IOBUFF: dup ERROR : vectFiles[fd] not NULL : %d %p\n", fd, vectFiles[fd]);
        exit(1);
      }
      vectFiles[fd] = lbfile;
    }
    _exit_debug_dup(fd);
    return fd;
  }//dup

/////////////////////////////////////////////////////////////////////////////////////////////// dup2
  int dup2(int oldfd, int newfd) throw() {
    _init_debug_dup2(oldfd, newfd);
    spin_wait_init();
    if( !g_supported_exec )  //Pass-through
      return  real_dup2(oldfd,newfd);
      
    int fd;
    list<bfile *> *bfileList;

    bfile * lbfile_new  = vectFiles[newfd];
    bfile * lbfile_old  = vectFiles[oldfd];
    //if ((lbfile_old && !lbfile_old->isBuffered) || (lbfile_new && !lbfile_new->isBuffered)){ 
    //  fprintf(debugOut,"IOBUFF: dup2 ERROR : not supported on files open with fopen : %d %p %d %p\n", newfd, lbfile_new,oldfd, lbfile_old);
    //  exit(1);
    //}
      
    //on se prepare a fermer ce file descriptor 
    if (lbfile_new)
      flush_buffer(newfd, lbfile_new);
    
    fd=real_dup2(oldfd,newfd);

    if (fd != -1 && oldfd != newfd) {

      //si les deux fd pointent sur des fichiers différents (sinon rien a faire)
      if (lbfile_old != lbfile_new) {
        //detruit le fdnew. Le fichier est deja fermer par dup2
        if (lbfile_new && newfd>=0){
          pthread_mutex_lock(&mutex);
          delist_bfile(newfd, lbfile_new); 
          pthread_mutex_unlock(&mutex); 
        }//if (lbfile_new && newfd>=0){...
        //si c'est un fichier qu'on "buffe", ajouter une reference
        if (lbfile_old){
          lbfile_old->buffer_ref++;
        }
        //pointe vers le meme buffer
        vectFiles[newfd] = lbfile_old;
        //vectFiles[oldfd] = NULL;  //????
      }//if (lbfile_old != lbfile_new) ... 
    }//if (fd != -1 && oldfd != newfd) ...
    _exit_debug_dup2(fd);
    return fd;
  }//dup2


/////////////////////////////////////////////////////////////////////////////////////////////// dup3
  int dup3(int oldfd, int newfd, int flags) throw() {
    _init_debug_dup3( oldfd, newfd);
    spin_wait_init();
    if( !g_supported_exec )  //Pass-through
      return real_dup3(oldfd,newfd,flags);
      
    int fd;

    // Comportement non-définit lorsque oldfd = newfd:
    if(oldfd == newfd)
      return real_dup3(oldfd,newfd,flags);
    list<bfile *> *bfileList;
    
    bfile * lbfile_new  = get_bfile(newfd);
    bfile * lbfile_old  = get_bfile(oldfd);
   // if ((lbfile_old && !lbfile_old->isBuffered) || (lbfile_new && !lbfile_new->isBuffered)){ 
   //   fprintf(debugOut,"IOBUFF: dup3 ERROR : not supported on files open with fopen : %d %p %d %p\n", newfd, lbfile_new,oldfd, lbfile_old);
   //   exit(1);
   // }
    //on se prepare a fermer ce file descriptor
    if (lbfile_new)
      flush_buffer(newfd, lbfile_new); 
    fd=real_dup3(oldfd,newfd,flags);
    
    if (fd != -1 && oldfd != newfd) {
      
      //si les deux fd pointent sur des fichiers différents (sinon rien a faire, perror...)
      if (lbfile_old != lbfile_new) {
        //detruit le fdnew. Le fichier est deja fermer par dup3
        if (lbfile_new && newfd>=0){
          pthread_mutex_lock(&mutex); 
          delist_bfile(newfd, lbfile_new); 
          pthread_mutex_unlock(&mutex); 
        }
        //si c'est un fichier qu'on "buffe", ajouter une reference
        if (lbfile_old){
          lbfile_old->buffer_ref++;
        }
        //pointe vers le meme buffer
        vectFiles[newfd] = lbfile_old;
      }//if (lbfile_old != lbfile_new) ...
    }//if (fd != -1 && oldfd != newfd) ...
    _exit_debug_dup3(fd);
    return fd;
  }//dup3

////////////////////////////////////////////////////////////////////////////////////////////// fcntl
  int fcntl(int fd, int cmd, ...) {
    _init_debug_fcntl();
    spin_wait_init();
    long farg=0;
#ifdef F_DUPFD_CLOEXEC
    if (cmd == F_DUPFD || cmd == F_DUPFD_CLOEXEC || cmd == F_SETFD || cmd == F_SETFL || cmd == F_GETLK || cmd == F_SETLK || cmd == F_SETLKW)
#else
    if (cmd == F_DUPFD || cmd == F_SETFD || cmd == F_SETFL || cmd == F_GETLK || cmd == F_SETLK || cmd == F_SETLKW)
#endif
    {
      va_list parg;
      va_start(parg, cmd);
      farg = va_arg(parg, long);
      va_end(parg);
    }
    if( !g_supported_exec )  //Pass-through
      return real_fcntl(fd, cmd, farg);

    int ret;
    bfile * lbfile = get_bfile (fd);
    _debug_fcntl(fd, cmd, farg);

    ret=real_fcntl(fd, cmd, farg);

#ifdef F_DUPFD_CLOEXEC
    if ((cmd == F_DUPFD || cmd == F_DUPFD_CLOEXEC) && ret != -1)
#else
    if ((cmd == F_DUPFD ) && ret != -1)
#endif
    {
      if (lbfile){
        lbfile->buffer_ref++;
      }
      if (vectFiles[ret] != NULL){ 
        fprintf(debugOut,"IOBUFF: fcntl error F_DUPFD or F_DUPFD_CLOEXEC : vectFiles[ret] != NULL");
        exit(1);
      }
      vectFiles[ret] = lbfile;
    }
    else if (cmd == F_GETLK) {
      fprintf(debugOut,"IOBUFF: fcntl not implemented for F_GETLK !!!! aborting.\n");
      //printf("IOBUFF: fcntl not implemented for F_GETLK !!!! aborting.\n");
      exit(1);
    }
    else if (cmd == F_SETLK) {
      struct flock *sfl = (struct flock *)farg;
      switch (sfl->l_type) {
      case F_RDLCK:
        break;
      case F_WRLCK:
        fprintf(debugOut,"IOBUFF: F_WRLCK not yet implemented !!!!!!!\n");
        //printf("IOBUFF: F_WRLCK not yet implemented !!!!!!!\n");
        exit(1);
        break;
      case F_UNLCK:
        break;
      }
    }
    else if (cmd == F_SETLKW)
    {
      fprintf(debugOut,"IOBUFF: fcntl not implemented for F_SETLKW !!!! aborting.\n");
      //printf("IOBUFF: fcntl not implemented for F_SETLKW !!!! aborting.\n");
      exit(1);
    }
    return ret;
  }//fcntl
/////////////////////////////////////////////////////////////////////////////////////////////
  /*
    int mknod(const char *pathname, mode_t mode, dev_t dev){
    printf("IOBUFF: mknod not implemented yet!!!! aborting.\n");
    exit(1);
    }
  
    void *mmap(void *addr, size_t length, int prot, int flags,int fd, off_t offset){
    printf("IOBUFF: mmap not implemented yet!!!! aborting.\n");
    exit(1);
    }
    void *mmap64(void *addr, size_t length, int prot, int flags,int fd, off64_t offset){
    printf("IOBUFF: mmap64 not implemented yet!!!! aborting.\n");
    exit(1);
    }
    int munmap(void *addr, size_t length){
    printf("IOBUFF: munmap not implemented yet!!!! aborting.\n");
    exit(1);
    }
  */
///////////////////////////////////////////////////////////////////////////////////////////// openat
  int openat(int dirfd, const char *pathname, int flags, ...){
    _init_debug_openat(dirfd);
    spin_wait_init();
    mode_t mode;
    
    if (flags & O_CREAT) {
      va_list arg;
      va_start(arg, flags);
      mode = va_arg(arg, mode_t);
      va_end(arg);
    }
    if( !g_supported_exec )  //Pass-through
      return real_openat(dirfd, pathname, flags, mode);
      
    int p,fd;
    char *cname; //nom canonique 
    ssize_t cnlen=0;
    list<bfile *> *bfileList;
    bfile * lbfile;
     
    fd=real_openat(dirfd, pathname, flags, mode);
 
    pthread_mutex_lock(&mutex);
    if (fd != -1 && vectFiles[fd] != NULL){ 
      fprintf(debugOut,"IOBUFF: open ERROR : vectFiles[fd] not NULL : %d %p\n", fd, vectFiles[fd]);
      fprintf(debugOut,"IOBUFF: new name %s \n", canonicalize_file_name(pathname));
      fprintf(debugOut,"IOBUFF: old name %s vectFiles[fd] \n", vectFiles[fd]->canonical_name);
      
      list<bfile *> *bfileList;
      bfileList = &fileMap[  std::string(vectFiles[fd]->canonical_name)  ];
      fprintf(debugOut,"IOBUFF: bflist %p %d vectFiles[fd] %p canonical_name %p\n", bfileList, 
                       bfileList->size(), vectFiles[fd], vectFiles[fd]->canonical_name);
      
      pthread_mutex_unlock(&mutex);
      exit(1);
    }
    else if (fd != -1 && vectFiles[fd] ==NULL) {
      char *fdpath;
      if(dirfd == AT_FDCWD){
        char *dirname= get_current_dir_name();
        fdpath = (char*)malloc(strlen(dirname) + strlen(pathname) +4);
        sprintf(fdpath,"%s/%s",dirname,pathname);
        free(dirname);
      }
      else{
        fdpath = (char*)malloc(strlen(pathname) +30);
        sprintf(fdpath,"/proc/self/fd/%d/%s",dirfd,pathname);
      }
      get_canonical_file_name(fdpath, &cname, &cnlen);   
      free(fdpath);
      if (exclude_path(&cname)){
        pthread_mutex_unlock(&mutex);
        return fd;
      }
#ifdef MPIBUFF_SUPPORTED
      if(std::find(mpibufFiles.begin(), mpibufFiles.end(), std::string(cname))!=mpibufFiles.end()){
        if ( (flags & O_ACCMODE) != O_RDONLY){
          fprintf(debugOut, "MPIBUF error : Only O_RDONLY is supported\n"); 
          exit(1);
        }
        //vectFiles[fd] = new bfile (fd, flags, cname, cnlen, true, true ); 
        vectFiles[fd] = new bfile (fd, flags, cname, cnlen, false, true ); 
      }
      else
#endif
        vectFiles[fd] = new bfile (fd, flags, cname, cnlen);
      free(cname);
      _debug_openat(fd, pathname, cname, cnlen, flags);
      //Ajout du fichier à la map des noms de fichier
      enlist_bfile(vectFiles[fd]);
    }
    pthread_mutex_unlock(&mutex);
    _exit_debug_openat(fd);
    return fd;
  }//openat

///////////////////////////////////////////////////////////////////////////////////////////// socket
  int socket(int domain, int type, int protocol) throw() {
    _init_debug_socket(domain, type, protocol);
    spin_wait_init();
    if( !g_supported_exec )  //Pass-through
      return real_socket(domain, type, protocol);
    int fd;
    bfile * lbfile;
    fd=real_socket(domain, type, protocol);
    _exit_debug_socket(fd);
    return fd;
  }//socket

  
////////////////////////////////////////////////////////////////////////////////////////////// close
  int close(int fd) {
    spin_wait_init();
    if( !g_supported_exec )  //Pass-through
      return real_close(fd);
    if (fd>=0) {
      pthread_mutex_lock(&mutex);
      bfile * lbfile = get_bfile(fd);
      if (lbfile) {
        flush_buffer(fd, lbfile); 
        delist_bfile(fd, lbfile);
        vectFiles[fd]=NULL;
      }
      pthread_mutex_unlock(&mutex);
    } 
    _exit_debug_close();
    return real_close(fd);
  }//close 
  
////////////////////////////////////////////////////////////////////////////////////////////// fopen
  FILE *fopen( const char * filename, const char * mode){
    _init_debug_fopen();
    spin_wait_init();
    if( !g_supported_exec )  //Pass-through
      return real_fopen( filename, mode );
      
    FILE *ret = real_fopen( filename, mode );
    if (ret != NULL) {
      int fd = fileno(ret); 
      _debug_fopen(fd); 
      pthread_mutex_lock(&mutex);    
      if (fd != -1 && vectFiles[fd] != NULL){ 
        fprintf(debugOut,"IOBUFF: fopen ERROR : vectFiles[fd] != NULL : %d %p\n", fd,vectFiles[fd]); 
        pthread_mutex_unlock(&mutex);
        exit(1);
      } 
      else if (fd != -1 && vectFiles[fd] ==NULL) {
        char *cname; //nom canonique 
        ssize_t cnlen=0;  
        get_canonical_file_name(filename, &cname, &cnlen);    
        if (exclude_path(&cname)){ 
          pthread_mutex_unlock(&mutex);
          return ret;
        }
        int flags= mode_to_flag(mode); 
#ifdef MPIBUFF_SUPPORTED
        if(std::find(mpibufFiles.begin(), mpibufFiles.end(), string(cname))!=mpibufFiles.end()){
          if ( (flags & O_ACCMODE) != O_RDONLY){
            fprintf(debugOut, "MPIBUF error : Only O_RDONLY is supported\n"); 
            exit(1);
          }
          vectFiles[fd] = new bfile (fd, flags, cname, cnlen, false, true ); 
        }
        else
#endif
          vectFiles[fd] = new bfile (fd, flags, cname, cnlen, false); 
        free(cname); 
        
        //Ajout du fichier à la map des noms de fichier 
        enlist_bfile(vectFiles[fd]);  
      } 
      pthread_mutex_unlock(&mutex);
      _exit_debug_close(fd);
    }
    return ret;
  }//fopen
  
  
///////////////////////////////////////////////////////////////////////////////////////////// fdopen
  FILE *fdopen( int fd, const char * mode){
    spin_wait_init();
    if( !g_supported_exec )  //Pass-through
      return real_fdopen(fd, mode);
      
#ifdef IOBUFF_FDOPEN_SUPPORTED
    bfile * lbfile = get_bfile(fd); 
    if (lbfile) { 
      pthread_mutex_lock(&mutex);
      if (lbfile->isBuffered){
        if(!flush_buffer(fd, lbfile) ){ 
          real_lseek(fd, lbfile->file_begin_buff + lbfile->buff_pos, SEEK_SET);
        }
        lbfile->isBuffered=false;
        if(lbfile->buffer != NULL){
          delete[] lbfile->buffer;
          lbfile->buffer = NULL;
        }
      }  
      pthread_mutex_unlock(&mutex);
    }
    FILE* ret = real_fdopen(fd, mode);
    return ret;
#else
  fprintf(debugOut,"IOBUFF: fdopen not supported\n");
  exit(1);
  return NULL;
#endif
  }//fdopen
//////////////////////////////////////////////////////////////////////////////////////////// freopen
  
  FILE *freopen( const char * filename, const char * mode, FILE *stream){
    _init_debug_freopen();
    spin_wait_init();
    if( !g_supported_exec )  //Pass-through
      return real_freopen( filename, mode, stream );
      
    int fd = fileno(stream); 
    if (fd>=0) { //First deallocate the bfile of the stream
      pthread_mutex_lock(&mutex);
      bfile * lbfile = get_bfile(fd); 
      _debug_freopen(fd);
      if (lbfile) {
        if(lbfile->buffer_ref == 1) { //Last duplicate of the file (dup) 
          list<bfile *> *bfileList= &fileMap[std::string(lbfile->canonical_name)]; 
          if(bfileList->size() ==1) {  //No other open copy of the file
            int tmp_cnlength = strlen(lbfile->canonical_name);
            char *tmp_cname = new char[tmp_cnlength+1];
            strncpy(tmp_cname, lbfile->canonical_name, tmp_cnlength);
            tmp_cname[tmp_cnlength] = '\0'; 
            bfileList->remove(lbfile); 
            fileMap.erase(std::string(tmp_cname));
            delete [] tmp_cname;
          }
          else { //File is still open on a different non-dup fd. Only remove this fd. 
            bfileList->remove(lbfile);
            delete_bfile(fd); 
          }
        delete_bfile(fd);
        }//if(lbfile->buffer_ref == 1)
        else {// Some duplicate of the file exist (dup)
          delete_bfile(fd);
        }
        vectFiles[fd]=NULL; 
      }
      pthread_mutex_unlock(&mutex); 
    }

    FILE* ret = real_freopen( filename, mode, stream );
    if (ret != NULL) {  //Then open the new bfile
      fd = fileno(ret); 

      pthread_mutex_lock(&mutex);    
      if (fd != -1 && vectFiles[fd] != NULL){
        fprintf(debugOut,"IOBUFF: freopen ERROR: vectFiles[fd] != NULL : %d %p\n",fd,vectFiles[fd]); 
        pthread_mutex_unlock(&mutex);
        exit(1);
      } 
      else if (fd != -1 && vectFiles[fd] ==NULL) {
        char *cname; //nom canonique 
        ssize_t cnlen=0;  
        get_canonical_file_name(filename, &cname, &cnlen);    
        if (exclude_path(&cname)){
          pthread_mutex_unlock(&mutex);
          return ret;
        }
        int flags= mode_to_flag(mode); 
#ifdef MPIBUFF_SUPPORTED
        if(std::find(mpibufFiles.begin(), mpibufFiles.end(), string(cname))!=mpibufFiles.end()){
          if ( (flags & O_ACCMODE) != O_RDONLY){
            fprintf(debugOut, "MPIBUF error : Only O_RDONLY is supported\n"); 
            exit(1);
          }
          vectFiles[fd] = new bfile (fd, flags, cname, cnlen, false, true ); 
        }
        else
#endif
          vectFiles[fd] = new bfile (fd, flags, cname, cnlen, false); 
        free(cname); 
        
        //Ajout du fichier à la map des noms de fichier 
        enlist_bfile(vectFiles[fd]);  
      }
      pthread_mutex_unlock(&mutex);
      _exit_debug_freopen(fd);
    }
  return ret;
  }
  
///////////////////////////////////////////////////////////////////////////////////////////// fclose
  int fclose(FILE* fp){
    _init_debug_fclose();
    spin_wait_init();
    if( !g_supported_exec )  //Pass-through
      return real_fclose(fp);

    int fd = fileno(fp); 
    
    list<bfile *> *bfileList;
    if (fd>=0) {
      _debug_fclose(fd);
      pthread_mutex_lock(&mutex);
      bfile * lbfile = get_bfile(fd); 
      if (lbfile) { 
        if(lbfile->buffer_ref == 1) { //Last duplicate of the file (dup) 
          bfileList = &fileMap[std::string(lbfile->canonical_name)]; 
          if(bfileList->size() ==1) {  //No other open copy of the file
            int tmp_cnlength = strlen(lbfile->canonical_name);
            char *tmp_cname = new char[tmp_cnlength+1];
            strncpy(tmp_cname, lbfile->canonical_name, tmp_cnlength);
            tmp_cname[tmp_cnlength] = '\0'; 
            bfileList->remove(lbfile); 
            fileMap.erase(std::string(tmp_cname));
            delete [] tmp_cname;
          } 
          else { //File is still open on a different non-dup fd. Only remove this fd. 
            bfileList->remove(lbfile);
            delete_bfile(fd); 
          }    
        delete_bfile(fd);
        }//if(lbfile->buffer_ref == 1)
        else {// Some duplicate of the file exist (dup)
          delete_bfile(fd);
        }
        vectFiles[fd]=NULL; 
      }  
      pthread_mutex_unlock(&mutex); 
      _exit_debug_fclose(fd);
    }
    return real_fclose(fp);
  }
  
}//extern c


