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
#ifndef _IOBUFF_GLOB_VARS
#define _IOBUFF_GLOB_VARS
#include <sys/uio.h>

bool g_supported_exec; //Current exec file is supported by iobuff

using namespace std;
#define IOBUFF_FDOPEN_SUPPORTED
//#define IOBUFF_DISALLOW_MIXED_OPEN_TYPES
#define IOBUFF_PREAD_PWRITE_SUPPORTED
#define IOBUFF_READV_WRITEV_SUPPORTED

#if defined(IOBUFF_FDOPEN_SUPPORTED) || defined(IOBUFF_READV_WRITEV_SUPPORTED)
#define IOBUFF_MUTEX_ON_ALL_OPS
#endif

#ifdef MPIBUFF_SUPPORTED
enum MPIBCommType{MPIB_OPEN=1, MPIB_READ=2, MPIB_CLOSE=3, MPIB_REGISTER=4, MPIB_USE_SMS=5};
std::vector<string> mpibufFiles;
int mpiSocketFD;

#define MPIBUFF_MULTI_THREADED

#endif //MPIBUFF_SUPPORTED

//#define MPIB_LOCK_DBG

//#define IOBUFF_TRACK
#ifdef IOBUFF_DEBUG
int nb_buff = 0;
int nb_buff2 = 0;
#define IOBUFF_DEBUG_ENTRY
//#define IOBUFF_DEBUG_CLOSE 
//#define IOBUFF_MAP_DUMP
//#define IOBUFF_DEBUG_BUFFER 
//#define IOBUFF_DEBUG_FORK
FILE *debugOut = stderr;
#else
FILE *debugOut = stderr;
#endif //IOBUFF_DEBUG

pid_t lp;
int lockState = 0;
int NB_FILES=32768; 
//buffer size en bytes (defaut 4M)
//int BSIZE = 4 * 1024 * 1024;
long globBuffUID = 0;
unsigned int BSIZE = 400  * 1024;

#ifdef MPIBUFF_SUPPORTED
unsigned int MPIBUFF_SMS_SIZE = BSIZE*10;
#endif

#ifdef IOBUFF_ALLOW_DEV
  const int NB_EXCLUDE = 6; // 6 to allow /dev/shm
  char EXCLUDE[NB_EXCLUDE][10] = { "/bin\0", "/lib\0", "/lib64\0", "/opt\0", "/proc\0", "/sys\0" };
#else 
  const int NB_EXCLUDE = 7; 
  char EXCLUDE[NB_EXCLUDE][10] = { "/bin\0", "/lib\0", "/lib64\0", "/opt\0", "/proc\0", "/sys\0", "/dev\0" };
#endif

//mutex pour les structure de donnees
pthread_mutex_t mutex = PTHREAD_ADAPTIVE_MUTEX_INITIALIZER_NP ;
pthread_mutex_t mutex2 = PTHREAD_ADAPTIVE_MUTEX_INITIALIZER_NP ; 
pthread_mutex_t mpibuff_mutex = PTHREAD_ADAPTIVE_MUTEX_INITIALIZER_NP ;
//pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER ;
///////////////////////////////////////////////////////////////////////////////////////////////////
//Declarations si nécessaire
void get_canonical_file_name(const char* , char** , ssize_t* );
#ifdef MPIBUFF_SUPPORTED
extern "C" {
ssize_t mpibuff_read(int , void *, size_t );
}
#endif //MPIBUFF_SUPPORTED


class bfile {
public:
  char  *buffer;
  char  *canonical_name;  //nom cannonique du fichier (unique)
  int    flags;
  int    buffer_ref;      //nombre de référence sur le buffer (fct: dup, dup2, etc)
  size_t file_begin_buff; //position du debut du buffer dans le fichier
  size_t file_end_buff;   //position de la fin du buffer dans le fichier
  size_t file_end_pos;    //position de la fin du fichier (nb de bytes)
  size_t buff_pos;        //position d'ecriture ou lecture dans le buffer
  size_t buff_end;        //nombre de bytes lus dans le buffer
  bool   buff_need_flush; //le buffer a ete modifie par une ecriture
  bool   buff_need_fill;  //le buffer a ete rempli partiellement par un write et 
                          //peut etre rempli completement par un read
  long buffUID;
  bool isBuffered;        //File is buffered or just refferenced
#ifdef MPIBUFF_SUPPORTED
  bool isMPIBuffered;     //MPI Buffered file
  int  fBufId;            //MPIBUF file id
  size_t mpibuff_pos;      //position de lecture dans le fichier
#endif
//  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  // bfile
#ifdef MPIBUFF_SUPPORTED
  bfile(int fd, int flags_in, char *cname, int cnlength, bool fileIsBuffered, bool fileIsMPIBuffered) ;
#else 
  bfile(int fd, int flags_in, char *cname, int cnlength, bool fileIsBuffered) ;
#endif

  ~bfile();
};  //bfile 
//  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  // bfile
bfile **vectFiles;


//pointeur sur les vrais fonctions
ssize_t (*real_write)(int,const void*,size_t)=NULL;
ssize_t (*real_pwrite)(int,const void*,size_t,off_t)=NULL;
ssize_t (*real_pwrite64)(int,const void*,size_t,off64_t)=NULL;
ssize_t (*real_pwritev)(int,const struct iovec *,int,off_t)=NULL;
ssize_t (*real_writev)(int,const struct iovec *,int)=NULL;
ssize_t (*real_read)(int,void*,size_t)=NULL;
ssize_t (*real_pread)(int,void*,size_t,off_t)=NULL;
ssize_t (*real_pread64)(int,void*,size_t,off64_t)=NULL;
ssize_t (*real_preadv)(int,const struct iovec *,int,off_t)=NULL;
ssize_t (*real_readv)(int,const struct iovec *,int)=NULL;
off_t   (*real_lseek)(int,off_t,int)=NULL;
loff_t  (*real_llseek)(int, loff_t, int)=NULL;
off64_t (*real_lseek64) (int, off64_t, int)=NULL;
int     (*real_open)(const char*,int,mode_t)=NULL;
int     (*real_open64)(const char*,int,mode_t)=NULL;
FILE*   (*real_fopen)( const char *, const char *)=NULL;
FILE*   (*real_fdopen)( int, const char *)=NULL;
FILE*   (*real_freopen)( const char *, const char *, FILE *)=NULL;
int     (*real_close)(int)=NULL;
int     (*real_fclose)(FILE*)=NULL;
int     (*real_dup)(int)=NULL;
int     (*real_dup2)(int,int)=NULL;
int     (*real_dup3)(int,int,int)=NULL;
int     (*real_fcntl)(int,int,long)=NULL;
int     (*real_openat)(int,const char*,int,mode_t)=NULL;
int     (*real_socket)(int,int,int)=NULL;


#endif  // _IOBUFF_GLOB_VARS