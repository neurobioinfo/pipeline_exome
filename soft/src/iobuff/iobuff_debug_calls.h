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
#ifndef _IOBUFF_DEBUG_CALLS
#define _IOBUFF_DEBUG_CALLS

/**************
 * Usage notes :
 *   None of the functions define here do anything by default.
 *   The body of each function contains some #ifdef ... #endif regions that
 *   must be defined at compile time to enable the different debug segments.
 *
 **************/

#include "iobuff_glob_vars.h"

///////////////////////////////////////////////////////////////////////////////////////// track_file
#ifdef IOBUFF_TRACK
void track_file(int fd, const char* msg) {
  if (fd >=0 && vectFiles[fd] != NULL && strcmp(vectFiles[fd]->canonical_name, 
    "/home/myHome/someFileName") ==0)
    fprintf(debugOut, "IOBUFF: track %s", msg);
} //track_file
#endif

//////////////////////////////////////////////////////////////////////////////////////////////
void  _init_debug_delete_bfile(int fd) {
#ifdef IOBUFF_DEBUG_ENTRY
   fprintf(debugOut, "IOBUFF: Entry delete_bfile\n");
#endif
#ifdef IOBUFF_DEBUG
   pid_t lp2 = (pid_t)syscall(__NR_gettid);
   if(fd>=0 && vectFiles[fd]){
     fprintf(debugOut, "IOBUFF:%d delete_bfile fd %d buffer_ref %d\n",lp2, fd, vectFiles[fd]->buffer_ref);
   }
#endif 
}
//////////////////////////////////////////////////////////////////////////////////////////////
void  _exit_debug_delete_bfile() {
#ifdef IOBUFF_DEBUG_ENTRY
   fprintf(debugOut, "IOBUFF: Exit delete_bfile\n");
#endif
}
//////////////////////////////////////////////////////////////////////////////////////////////
void  _init_debug_enlist_bfile(list<bfile *> *bfileList){
#ifdef IOBUFF_DEBUG_ENTRY
    fprintf(debugOut, "IOBUFF Entry: enlist_bfile\n");
#endif 
#ifdef IOBUFF_DEBUG
  fprintf(debugOut,"IOBUFF: bflist %p %d\n",bfileList, bfileList->size());
#endif
}
//////////////////////////////////////////////////////////////////////////////////////////////
void  _exit_debug_enlist_bfile(list<bfile *> *bfileList){
#ifdef IOBUFF_DEBUG
  fprintf(debugOut,"IOBUFF: bflist %p %d\n",bfileList, bfileList->size());
#endif
#ifdef IOBUFF_DEBUG_ENTRY
    fprintf(debugOut, "IOBUFF Exit: enlist_bfile\n");
#endif 
}
//////////////////////////////////////////////////////////////////////////////////////////////
void _init_debug_delist_bfile(){
#ifdef IOBUFF_DEBUG_ENTRY
    fprintf(debugOut, "IOBUFF Entry: delist_bfile\n");
#endif 
}
//////////////////////////////////////////////////////////////////////////////////////////////
void _debug_delist_bfile(list<bfile *> *bfileList){
#ifdef IOBUFF_DEBUG_DELIST
    pid_t lp3 = (pid_t)syscall(__NR_gettid);       
    fprintf(debugOut,"IOBUFF: bflist %p\n", bfileList);
    list<bfile*>::iterator it = find( bfileList->begin(), bfileList->end(), lbfile);
    fprintf(debugOut,"IOBUFF:%d delist fd  \"%s\" %p FileList %p size %d\n", lp3,fd, 
            lbfile->canonical_name,lbfile, bfileList,bfileList->size());
    if (it == bfileList->end() ){
      fprintf(debugOut,"IOBUFF:%d BFile not found  %p\n",lp3,  lbfile);
      for(it = bfileList->begin(); it != bfileList->end(); ++it)
        fprintf(debugOut,"IOBUFF:%d BFile from list %p %s\n",lp3, (*it),  (*it)->canonical_name); 
      fprintf(debugOut,"IOBUFF delist:%d Unocking %d\n",lp2, --lockState); 
      pthread_mutex_unlock(&mutex);
      exit(1);
    }
#endif
}
//////////////////////////////////////////////////////////////////////////////////////////////
void _exit_debug_delist_bfile(){
#ifdef IOBUFF_DEBUG_ENTRY
    fprintf(debugOut, "IOBUFF Exit: delist_bfile\n");
#endif 
}
//////////////////////////////////////////////////////////////////////////////////////////////
void _init_debug_open64(){
#ifdef IOBUFF_DEBUG_ENTRY
  fprintf(debugOut, "IOBUFF Entry: open64\n");
#endif
}
//////////////////////////////////////////////////////////////////////////////////////////////
void _init_debug_open(){
#ifdef IOBUFF_DEBUG_ENTRY
  fprintf(debugOut, "IOBUFF Entry: open\n");
#endif
#ifdef IOBUFF_DEBUG
  pid_t lp2 = (pid_t)syscall(__NR_gettid);
#endif
#ifdef IOBUFF_MAP_DUMP
  map <std::string, list<bfile *>>::iterator ifts;
  fprintf(debugOut, "OPEN: ");
  for (ifts=fileMap.begin(); ifts!=fileMap.end(); ++ifts){
    fprintf(debugOut, "(%p)-%p ", &(ifts->first), &(ifts->second));
  }
  fprintf(debugOut, "\n");
#endif
}
//////////////////////////////////////////////////////////////////////////////////////////////
void _debug_open(int fd, const char* pathname, int flags){
#ifdef IOBUFF_TRACK
  track_file(fd, "pre open");
#endif
#ifdef IOBUFF_DEBUG
  pid_t lp2 = (pid_t)syscall(__NR_gettid);
  char cflags[12];
  sprintf(cflags,"%d",flags);
  if ((flags & O_ACCMODE) == O_WRONLY) strcpy(cflags, "O_WRONLY");
  else if ( (flags & O_ACCMODE) == O_RDONLY) strcpy(cflags, "O_RDONLY");
  else if ( (flags & O_ACCMODE) == O_RDWR) strcpy(cflags, "O_RDWR");
  fprintf(debugOut,"IOBUFF:%d open ,%s, fd %d, flag %s, %p\n", lp2, pathname, fd, cflags, 
                   vectFiles[fd]);
#endif
}
//////////////////////////////////////////////////////////////////////////////////////////////
void _debug_open2(int fd, const char *pathname, const char *cname, int cnlen){
#ifdef IOBUFF_DEBUG
  //printf("IOBUFF:%d Canonical name ,%s, sb %d, sa %d, len %d\n",lp, 
  //        vectFiles[fd]->canonical_name, strlen(pathname), strlen(cname), cnlen);
  pid_t lp3 = (pid_t)syscall(__NR_gettid);
  printf("IOBUFF:%d open fd %d, Cname ,%s, sb %d, sa %d, len %d %p\n",lp3, fd, 
          vectFiles[fd]->canonical_name, strlen(pathname), strlen(cname), cnlen, vectFiles[fd]);
#endif
  
}
//////////////////////////////////////////////////////////////////////////////////////////////
void _exit_debug_open(int fd){
#ifdef IOBUFF_TRACK
    track_file(fd, "open");
#endif
#ifdef IOBUFF_DEBUG_ENTRY
    fprintf(debugOut, "IOBUFF Exit: open\n");
#endif 
}
//////////////////////////////////////////////////////////////////////////////////////////////
void _init_debug_mpibuff_lseek64(int fd, int offset, int whence){
#ifdef IOBUFF_DEBUG_ENTRY
    fprintf(debugOut, "IOBUFF Entry: mpibuff_lseek64 %d %d %d", fd, (int) offset, whence);
#endif  
}
//////////////////////////////////////////////////////////////////////////////////////////////
void _init_debug_pread(int fd, size_t count, off_t offset){
#ifdef IOBUFF_DEBUG_ENTRY
    fprintf(debugOut, "IOBUFF Entry: pread\n");
#endif
#ifdef IOBUFF_DEBUG_PREAD_PWRITE
    pid_t lp2 = (pid_t)syscall(__NR_gettid);
    printf("IOBUFF:%d pread %d %d %d\n", lp2, fd, (int)count, (int)offset);
#endif
}
//////////////////////////////////////////////////////////////////////////////////////////////
void _exit_debug_pread(){
#ifdef IOBUFF_DEBUG_ENTRY
    fprintf(debugOut, "IOBUFF Exit: pread\n");
#endif 
}
//////////////////////////////////////////////////////////////////////////////////////////////
void _init_debug_pread64(int fd, size_t count, off64_t offset){
#ifdef IOBUFF_DEBUG_ENTRY
    fprintf(debugOut, "IOBUFF Entry: pread64\n");
#endif
#ifdef IOBUFF_DEBUG_PREAD_PWRITE
    pid_t lp2 = (pid_t)syscall(__NR_gettid);
    printf("IOBUFF:%d pread64 %d %d %d\n", lp2, fd, count, offset);
#endif
}
//////////////////////////////////////////////////////////////////////////////////////////////
void _exit_debug_pread64(){
#ifdef IOBUFF_DEBUG_ENTRY
    fprintf(debugOut, "IOBUFF Exit: pread64\n");
#endif 
}
//////////////////////////////////////////////////////////////////////////////////////////////
void _init_debug_read(int fd, void *buf, size_t count){
#ifdef IOBUFF_DEBUG_ENTRY
  fprintf(debugOut, "IOBUFF Entry: read\n");
#endif
#ifdef IOBUFF_TRACK
    track_file(fd, "read");
#endif
}
//////////////////////////////////////////////////////////////////////////////////////////////
void _exit_debug_read(){
#ifdef IOBUFF_DEBUG_ENTRY
    fprintf(debugOut, "IOBUFF Exit: read\n");
#endif 
}
//////////////////////////////////////////////////////////////////////////////////////////////
void _init_debug_pwrite(int fd, size_t count, off_t offset){
#ifdef IOBUFF_DEBUG_ENTRY
    fprintf(debugOut, "IOBUFF Entry: pwrite\n");
#endif
#ifdef IOBUFF_DEBUG_PREAD_PWRITE
    pid_t lp2 = (pid_t)syscall(__NR_gettid);
    printf("IOBUFF:%d pwrite %d %d %d\n", lp2, (int)fd, (int)count, (int)offset);
#endif
  
}
//////////////////////////////////////////////////////////////////////////////////////////////
void _init_debug_pwrite64(int fd, size_t count, off64_t offset){
#ifdef IOBUFF_DEBUG_ENTRY
    fprintf(debugOut, "IOBUFF Entry: pwrite64\n");
#endif
#ifdef IOBUFF_DEBUG_PREAD_PWRITE
    pid_t lp2 = (pid_t)syscall(__NR_gettid);
    printf("IOBUFF:%d pwrite64 %d %d %d\n", lp2, (int)fd, (int)count, (int)offset);
#endif
}
//////////////////////////////////////////////////////////////////////////////////////////////
void _init_debug_write(int fd, size_t count){
#ifdef IOBUFF_DEBUG_ENTRY
  fprintf(debugOut, "IOBUFF Entry: write\n");
#endif
#ifdef IOBUFF_TRACK
    track_file(fd, "write");
#endif
  
}
//////////////////////////////////////////////////////////////////////////////////////////////
void _init_debug_lseek(int fd, off_t offset, int whence) {
#ifdef IOBUFF_DEBUG_ENTRY
  fprintf(debugOut, "IOBUFF Entry: lseek\n");
#endif
#ifdef IOBUFF_DEBUG_LSEEK
  pid_t lp2 = (pid_t)syscall(__NR_gettid);
    printf("IOBUFF:%d lseek %d %d %d\n", lp2, fd, offset, whence);
#endif

}
//////////////////////////////////////////////////////////////////////////////////////////////
void _init_debug_lseek64(int fd, off64_t offset, int whence) {
#ifdef IOBUFF_DEBUG_ENTRY
  fprintf(debugOut, "IOBUFF Entry: lseek64\n");
#endif
#ifdef IOBUFF_DEBUG_LSEEK
  pid_t lp2 = (pid_t)syscall(__NR_gettid);
    printf("IOBUFF:%d lseek64 %d %d %d %x\n", lp2,fd, offset, whence, lbfile);
#endif
#ifdef IOBUFF_TRACK
    track_file(fd, "lseek64");
#endif
}
//////////////////////////////////////////////////////////////////////////////////////////////
void _init_debug_dup(int oldfd) {
#ifdef IOBUFF_DEBUG_ENTRY
  fprintf(debugOut, "IOBUFF Entry: dup\n");
#endif
#ifdef IOBUFF_DEBUG
  pid_t lp2 = (pid_t)syscall(__NR_gettid);
    printf("IOBUFF:%d dup %d\n",lp2, oldfd);
#endif
#ifdef IOBUFF_TRACK
    track_file(oldfd, "old fd dup");
#endif
}
//////////////////////////////////////////////////////////////////////////////////////////////
void _exit_debug_dup(int fd) {
#ifdef IOBUFF_TRACK
    track_file(fd, "new fd dup");
#endif
}
//////////////////////////////////////////////////////////////////////////////////////////////
void _init_debug_dup2(int oldfd, int newfd) {
#ifdef IOBUFF_DEBUG_ENTRY
  fprintf(debugOut, "IOBUFF Entry: dup2\n");
#endif
#ifdef IOBUFF_DEBUG
  pid_t lp2 = (pid_t)syscall(__NR_gettid);
    printf("IOBUFF:%d dup2 %d %d\n",lp2, oldfd, newfd);
#endif
#ifdef IOBUFF_TRACK
    track_file(oldfd, "old fd dup2");
#endif
}
//////////////////////////////////////////////////////////////////////////////////////////////
void _exit_debug_dup2(int fd) {
#ifdef IOBUFF_TRACK
    track_file(fd, "new fd dup2");
#endif
}
//////////////////////////////////////////////////////////////////////////////////////////////
void _init_debug_dup3(int oldfd, int newfd) {
#ifdef IOBUFF_DEBUG_ENTRY
  fprintf(debugOut, "IOBUFF Entry: dup3\n");
#endif
#ifdef IOBUFF_DEBUG
  pid_t lp2 = (pid_t)syscall(__NR_gettid);
    printf("IOBUFF:%d dup3 %d %d\n", lp2, oldfd, newfd);
#endif
#ifdef IOBUFF_TRACK
    track_file(oldfd, "old fd dup3");
    track_file(newfd, "new fd pre dup3");
#endif
}
//////////////////////////////////////////////////////////////////////////////////////////////
void _exit_debug_dup3(int fd) {
#ifdef IOBUFF_TRACK
    track_file(fd, "new fd dup3");
#endif
}
//////////////////////////////////////////////////////////////////////////////////////////////
void _init_debug_fcntl(){
#ifdef IOBUFF_DEBUG_ENTRY
  fprintf(debugOut, "IOBUFF Entry: fcntl\n");
#endif
  
}
//////////////////////////////////////////////////////////////////////////////////////////////
void _debug_fcntl(int fd, int cmd, long farg){
#ifdef IOBUFF_MAP_DUMP
    map <std::string, list<bfile *>>::iterator ifts;
    fprintf(debugOut, "FCNTL: ");
    for (ifts=fileMap.begin(); ifts!=fileMap.end(); ++ifts){
      fprintf(debugOut, "(%p)-%p ", &(ifts->first), &(ifts->second));
    }
    fprintf(debugOut, "\n");
#endif
#ifdef IOBUFF_DEBUG
  pid_t lp2 = (pid_t)syscall(__NR_gettid);
  printf("IOBUFF:%d fcntl %d %d %d ", lp2, fd, cmd, farg);
  switch(cmd){
    case F_DUPFD: printf("F_DUPFD\n"); break;
    case F_GETFD: printf("F_GETFD\n"); break;
    case F_SETFD: printf("F_SETFD\n"); break;
    case F_SETFL: printf("F_SETFL\n"); break;
    case F_GETLK: printf("F_GETLK\n"); break;
    case F_SETLK: printf("F_SETLK\n"); break;
    case F_SETLKW: printf("F_SETLKW\n"); break;
#ifdef F_DUPFD_CLOEXEC
    case F_DUPFD_CLOEXEC: printf("F_DUPFD_CLOEXEC\n"); break;
#endif
    default : printf("UNKNOWN\n"); 
  }
#endif
#ifdef IOBUFF_TRACK
    track_file(fd, "fcntl");
    switch(cmd){
    case F_DUPFD: track_file(fd,"F_DUPFD\n"); break;
    case F_GETFD: track_file(fd,"F_GETFD\n"); break;
    case F_SETFD: track_file(fd,"F_SETFD\n"); break;
    case F_SETFL: track_file(fd,"F_SETFL\n"); break;
    case F_GETLK: track_file(fd,"F_GETLK\n"); break;
    case F_SETLK: track_file(fd,"F_SETLK\n"); break;
    case F_SETLKW: track_file(fd,"F_SETLKW\n"); break;
#ifdef F_DUPFD_CLOEXEC
    case F_DUPFD_CLOEXEC: track_file(fd,"F_DUPFD_CLOEXEC\n"); break;
#endif
    default : track_file(fd,"UNKNOWN\n"); 
  }
#endif
}
//////////////////////////////////////////////////////////////////////////////////////////////
void _init_debug_openat(int dirfd){
#ifdef IOBUFF_DEBUG_ENTRY
  fprintf(debugOut, "IOBUFF Entry: openat\n");
#endif
#ifdef IOBUFF_DEBUG
  pid_t lp2 = (pid_t)syscall(__NR_gettid);
    printf("IOBUFF:%d openat %d\n",lp2, dirfd);
#endif
#ifdef IOBUFF_MAP_DUMP
    map <std::string, list<bfile *>>::iterator ifts;
    fprintf(debugOut, "OPENAT: ");
    for (ifts=fileMap.begin(); ifts!=fileMap.end(); ++ifts){
      fprintf(debugOut, "(%p)-%p ", &(ifts->first), &(ifts->second));
    }
    fprintf(debugOut, "\n");
#endif
}
//////////////////////////////////////////////////////////////////////////////////////////////
void _debug_openat(int fd, const char* pathname, const char* cname, int cnlen, int flags){
#ifdef IOBUFF_TRACK
  track_file(fd, "post real_openat");
#endif
#ifdef IOBUFF_DEBUG
  fprintf(debugOut,"IOBUFF openat: cname %s",cname);  
  char cflags[12] = "";
  pid_t lp2 = (pid_t)syscall(__NR_gettid);
  if (flags & O_WRONLY) strcpy(cflags, "O_WRONLY");
  else if ( flags & O_RDONLY) strcpy(cflags, "O_RDONLY");
  else if ( flags & O_RDWR) strcpy(cflags, "O_RDWR");
  printf("IOBUFF:%d openat ,%s, fd %d, flag %s\n", lp2, pathname, fd, cflags);
  printf("IOBUFF:%d openat Canonical name ,%s, sb %d, sa %d, len %d\n", lp2, 
         vectFiles[fd]->canonical_name), strlen(pathname), strlen(cname), cnlen;
#endif 
}
//////////////////////////////////////////////////////////////////////////////////////////////
void _exit_debug_openat(int fd) {
#ifdef IOBUFF_TRACK
  track_file(fd, "openat");
#endif 
#ifdef IOBUFF_DEBUG_ENTRY
  fprintf(debugOut, "IOBUFF Exit: openat\n");
#endif
}
//////////////////////////////////////////////////////////////////////////////////////////////
void _init_debug_socket(int domain, int type, int protocol){
#ifdef IOBUFF_DEBUG_ENTRY
  fprintf(debugOut, "IOBUFF Entry: socket\n");
#endif
#ifdef IOBUFF_DEBUG
  pid_t lp2 = (pid_t)syscall(__NR_gettid);
    printf("IOBUFF:%d socket %d %d %d\n", lp2,domain, type, protocol);
#endif
}
//////////////////////////////////////////////////////////////////////////////////////////////
void _exit_debug_socket(int fd){
#ifdef IOBUFF_TRACK
    track_file(fd, "socket");
#endif
#ifdef IOBUFF_DEBUG_ENTRY
  fprintf(debugOut, "IOBUFF Exit: socket\n");
#endif
}
//////////////////////////////////////////////////////////////////////////////////////////////
void _init_debug_close(int fd){
#ifdef IOBUFF_DEBUG_ENTRY
  fprintf(debugOut, "IOBUFF Entry: close\n");
#endif
#ifdef IOBUFF_DEBUG
  pid_t lp2 = (pid_t)syscall(__NR_gettid);
  fprintf(debugOut, "IOBUFF:%d close fd %d\n", lp2,fd);
#endif 
#ifdef IOBUFF_TRACK
    track_file(fd, "close");
#endif
#ifdef IOBUFF_MAP_DUMP
  map <std::string, list<bfile *>>::iterator ifts;
  fprintf(debugOut, "CLOSE: ");
  for (ifts=fileMap.begin(); ifts!=fileMap.end(); ++ifts){
    fprintf(debugOut, "(%p)-%p ", &(ifts->first), &(ifts->second));
  }
  fprintf(debugOut, "\n");
#endif
}
//////////////////////////////////////////////////////////////////////////////////////////////
void  _exit_debug_close(){
#ifdef IOBUFF_DEBUG_ENTRY
  fprintf(debugOut, "IOBUFF Exit: close\n");
#endif
}
//////////////////////////////////////////////////////////////////////////////////////////////
void _init_debug_fopen(){
#ifdef IOBUFF_DEBUG_ENTRY
  fprintf(debugOut, "IOBUFF Entry: fopen\n");
#endif
   
}
//////////////////////////////////////////////////////////////////////////////////////////////
void _debug_fopen(int fd){
#ifdef IOBUFF_TRACK
  track_file(fd, "pre fopen");
#endif
}
//////////////////////////////////////////////////////////////////////////////////////////////
void  _exit_debug_close(int fd){
#ifdef IOBUFF_TRACK
      track_file(fd, "fopen");
#endif 
#ifdef IOBUFF_DEBUG_ENTRY
    fprintf(debugOut, "IOBUFF Exit: fopen\n");
#endif 
}
//////////////////////////////////////////////////////////////////////////////////////////////
void _init_debug_freopen( ){
#ifdef IOBUFF_DEBUG_ENTRY
    fprintf(debugOut, "IOBUFF Entry: freopen\n");
#endif 
}
//////////////////////////////////////////////////////////////////////////////////////////////
void  _debug_freopen(int fd){
#ifdef IOBUFF_TRACK
    track_file(fd, "freopen");
#endif 
}
//////////////////////////////////////////////////////////////////////////////////////////////
void  _exit_debug_freopen(int fd){
#ifdef IOBUFF_TRACK
    track_file(fd, "freopen");
#endif 
#ifdef IOBUFF_DEBUG_ENTRY
    fprintf(debugOut, "IOBUFF Exit: freopen\n");
#endif 
}
//////////////////////////////////////////////////////////////////////////////////////////////
void  _init_debug_fclose(){
#ifdef IOBUFF_DEBUG_ENTRY
  fprintf(debugOut, "IOBUFF Entry: fclose\n");
#endif
}
//////////////////////////////////////////////////////////////////////////////////////////////
void _debug_fclose(int fd){
#ifdef IOBUFF_TRACK
    track_file(fd, "fclose");
#endif
#ifdef IOBUFF_DEBUG
  pid_t lp2 = (pid_t)syscall(__NR_gettid);
  fprintf(debugOut,"IOBUFF fclose:%d  \n",lp2);
#endif 
}
//////////////////////////////////////////////////////////////////////////////////////////////
void _exit_debug_fclose(int fd){
#ifdef IOBUFF_TRACK
    track_file(fd, "fclose");
#endif 
#ifdef IOBUFF_DEBUG_ENTRY
    fprintf(debugOut, "IOBUFF Exit: fclose\n");
#endif 
}
#endif // _IOBUFF_DEBUG_CALLS