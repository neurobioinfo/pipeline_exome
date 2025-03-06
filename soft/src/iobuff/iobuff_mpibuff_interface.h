#ifndef _IOBUFF_MPIBUFF_INTERFACE
#define _IOBUFF_MPIBUFF_INTERFACE

#ifdef MPIBUFF_SUPPORTED
int shm_id;  //shared mem segment id
int sem_id;  //semaphore id 
char* shm_addr; //shared mem address
const int local_sem_rd = 0;
const int dist_sem_wr  = 1;
const int local_sem_wr = 2;
const int dist_sem_rd  = 3;
long l0= 0;
long l1= 0;
long l2= 0;
long u1= 0;
long u3= 0;

//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ lock_sem
void lock_sem(int i){
  struct sembuf sb;
#ifdef MPIB_LOCK_DBG
  fprintf(debugOut, "L%d\n", i);
#endif
  sb.sem_num = i;
  sb.sem_op = -1;
  sb.sem_flg = 0;//SEM_UNDO;
  if(semop(sem_id, &sb, 1) != 0 ) 
    {fprintf(stderr,"IOBUFF: Can't lock semaphore %d ", i); perror(""); exit(1);}; 
} 

//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ unlock_sem
void unlock_sem(int i){
  struct sembuf sb;
#ifdef MPIB_LOCK_DBG
  fprintf(debugOut, "U%d\n", i);
#endif
  sb.sem_num = i;
  sb.sem_op = 1;
  sb.sem_flg = 0;//SEM_UNDO;
  if( semop(sem_id, &sb, 1) != 0 ) 
    {fprintf(stderr,"IOBUFF: Can't unlock semaphore %d ", i); perror(""); exit(1);}; 
}
#if !defined(__GNU_LIBRARY__) || defined(_SEM_SEMUN_UNDEFINED)
union semun
{
int val; // value for SETVAL
struct semid_ds* buf; // buffer for IPC_STAT, IPC_SET
unsigned short* array; // array for GETALL, SETALL
struct seminfo* __buf; // buffer for IPC_INFO
};
#endif //!defined

//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ shm_write
ssize_t shm_write(char* buf, size_t count){
  size_t wrote=0; 
  //fprintf(stderr, "L1 %ld\n", l1++);
  //lock_sem(dist_sem_wr);      //L1
  while (wrote < count){
    //fprintf(stderr, "L2 %ld\n", l2++);
    lock_sem(local_sem_wr);   //L2
    size_t to_write = MPIBUFF_SMS_SIZE < count-wrote ? MPIBUFF_SMS_SIZE : count-wrote; 
    memcpy(shm_addr, buf+wrote, to_write);
    wrote += to_write;
    //fprintf(stderr, "U3 %ld\n", u3++);
    unlock_sem(dist_sem_rd);  //U3
  }
  //fprintf(stderr, "U1 %ld\n", u1++);
  //unlock_sem(dist_sem_wr);      //U1
  return wrote;
}
//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ shm_read
ssize_t shm_read(char* buf, size_t count){
  size_t rd=0; 
  while (rd < count){
    //fprintf(stderr, "L0 %ld\n", l0++);
    lock_sem(local_sem_rd);   //L0
    size_t to_read = MPIBUFF_SMS_SIZE < count-rd ? MPIBUFF_SMS_SIZE : count-rd; 
    memcpy(buf+rd, shm_addr, to_read);
    rd += to_read;
    //fprintf(stderr, "U1 %ld\n", u1++);
    unlock_sem(dist_sem_wr);  //U1
  }
  return rd;
}
//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ shm_read2
ssize_t shm_read2(char* buf1, size_t count1, char* buf2, size_t count2){
  size_t rd=0, rd2=0; 
  size_t count = count1+count2;
  size_t to_read1, to_read2;
  while (rd < count){
    //fprintf(stderr, "L0 %ld, c1 %ld, c2 %ld, rd %ld, c %ld\n", l0, count1, count2, rd, count);
    lock_sem(local_sem_rd);   //L0
    to_read1 = to_read2 = 0;
    if(rd < count1){
      to_read1 = MPIBUFF_SMS_SIZE < count1-rd ? MPIBUFF_SMS_SIZE : count1-rd; 
      memcpy(buf1+rd, shm_addr, to_read1);
      rd += to_read1;
      count2 = ((long*)buf1)[0];
      count = count1+count2;
    }
    if(rd >= count1){
      to_read2 = MPIBUFF_SMS_SIZE < count-rd ? MPIBUFF_SMS_SIZE-to_read1 : count-rd; 
      memcpy(buf2+rd2, shm_addr+to_read1, to_read2);
      rd += to_read2;
      rd2 += to_read2;
    }
    //fprintf(stderr, "U1 %ld\n", u1++);
    unlock_sem(dist_sem_wr);  //U1
  } 
  return rd2;
}
//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ _open_mpibuff_file
int _open_mpibuff_file(char *canonical_name, int cnlength){
  long msg_buf[4];
  int fBufId;
#ifdef MPIBUFF_MULTI_THREADED
  pthread_mutex_lock(&mpibuff_mutex);
#endif
  //Enregistre le fichier auprès de MPIBUF
  msg_buf[0] = MPIB_OPEN; msg_buf[1] = cnlength+1;
  msg_buf[2] = 0;         msg_buf[3] = 0;
  ssize_t n = shm_write((char*) msg_buf, 4*sizeof(long)); 
      
  n = shm_write((char*) canonical_name, cnlength+1);  

  n = shm_read((char*) (&fBufId), sizeof(int)); 
  //fprintf(debugOut, "IOBUFF: fid %d str_fid %s\n",fBufId, (char*) &fBufId); exit(1);
#ifdef MPIBUFF_MULTI_THREADED
  pthread_mutex_unlock(&mpibuff_mutex);
#endif
  return fBufId;
}

//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ _init_mpibuff_support
int _init_mpibuff_support() {
  int mpiSocketFD = -1;
  const char* c_env_mpibuff_list = getenv("MPIBUF_FILE_LIST");
  if( c_env_mpibuff_list != NULL ) {
    // Some files are MPI-Bufferable. Get their canonical names
    const std::string env_mpibuff_list( c_env_mpibuff_list );
    const char delimiter = ':';
    size_t previous = 0;
    size_t index = env_mpibuff_list.find( delimiter );
    char *cname; ssize_t cnlen;
    while( index != string::npos ){
      get_canonical_file_name(env_mpibuff_list.substr(previous, index-previous).c_str(), 
        &cname, &cnlen);
      mpibufFiles.push_back( std::string(cname) );
      previous=index+1;
      index = env_mpibuff_list.find( delimiter, previous );
    } 
    get_canonical_file_name(env_mpibuff_list.substr(previous).c_str(), &cname, &cnlen);
    mpibufFiles.push_back( std::string(cname) );
    // Open the socket to the local mpibuf instance.
    struct sockaddr_un serv_addr; 
    mpiSocketFD = real_socket(AF_UNIX, SOCK_STREAM, 0);
    if( mpiSocketFD < 0 ){ 
      fprintf(debugOut, "IOBUFF Error : Could not create socket");
      exit(1);
    } 
    memset(&serv_addr, '0', sizeof(serv_addr));  
    const char *addrChar = getenv("MPIBUF_ADDRESS");
    serv_addr.sun_family = AF_UNIX;
    if(addrChar != NULL)
      strcpy(serv_addr.sun_path, addrChar);
    else
      strcpy(serv_addr.sun_path, "/tmp/mpibuff_sock");
    
    if( connect(mpiSocketFD, (struct sockaddr *)&serv_addr, sizeof(serv_addr)) < 0 ){
      perror("IOBUFF Error : Could not create socket \n");
      exit(1);
    } 
    
    if( mpiSocketFD != -1 ){
      long msg_buf[4];
      msg_buf[0] = MPIB_USE_SMS;
      char *sms_size = getenv("MPIBUF_SMS_SIZE");
      if(sms_size!=NULL)
        MPIBUFF_SMS_SIZE =atoi(sms_size);
      else 
        MPIBUFF_SMS_SIZE = BSIZE*10;
      msg_buf[1] = MPIBUFF_SMS_SIZE;
      int shm_id = shmget(IPC_PRIVATE, MPIBUFF_SMS_SIZE, 0666);
      if (shm_id < 0){ perror("IOBUFF Error: Can't create shared mem segment"); exit(1); }
      msg_buf[2] = shm_id;
      shm_addr= (char*)  shmat(shm_id, NULL, 0);
      
      sem_id = semget(IPC_PRIVATE, 4, SHM_R|SHM_W);  //Create a 2-value semaphore
      if (sem_id < 0){ perror("IOBUFF Error: Can't create semaphore"); exit(1); }
      msg_buf[3] = sem_id;
      union semun arg;
      unsigned short initv[4] = {0,0,0,0};  //Both semaphore start locked
      arg.array = initv;
      semctl(sem_id, 0, SETALL, arg);
  
      int n = real_write(mpiSocketFD, msg_buf, 4*sizeof(long));
      //lock_sem(local_sem);
      shmctl(shm_id, IPC_RMID, NULL); //Mark it for destruction
      //unlock_sem(dist_sem);
    } 
  } //End MPIBUF init section
  return mpiSocketFD;
}
//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ _close_mpibuff_support
void _close_mpibuff_support(int mpiSocketFD) {
  if( mpiSocketFD != -1 ){
    long msg_buf[4];
    msg_buf[0] = MPIB_CLOSE;
    int n = shm_write((char*) msg_buf, 4*sizeof(long));
    n = shmdt(shm_addr);
    if (n == -1) { perror("IOBUFF: Error detaching shared mem seg"); exit(1); }
    n = semctl(sem_id, 0, IPC_RMID, NULL);
    if (n == -1) { perror("IOBUFF: Error releasing semaphore"); exit(1); }
    close( mpiSocketFD );
  }
}

//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ mpibuff_pread
ssize_t mpibuff_pread(int fd, void *buf, size_t count, off_t offset){
  int remoteFID= vectFiles[fd]->fBufId;
  long msg_buf[4];
  long n;
  double time_spent;
  struct timeval tv1, tv2;
  msg_buf[0] = MPIB_READ; msg_buf[1] = remoteFID; msg_buf[2] = offset; msg_buf[3] = count;
#ifdef MPIBUFF_MULTI_THREADED
    pthread_mutex_lock(&mpibuff_mutex);
#endif
  n = shm_write((char*) msg_buf, 4*sizeof(long));
      
  ssize_t readSize=0;
  long readCountAndErr[2];  //read count and errno of last mpi read
  
  //\\n = shm_read((char*) readCountAndErr, 2*sizeof(long));
  n = shm_read2((char*) readCountAndErr, 2*sizeof(long), (char*)(buf), count);
  if( readCountAndErr[0] < 0 ){
    errno = readCountAndErr[1];
    perror("MPIBUF: Distant read error from server"); close(mpiSocketFD); exit(1); 
  }
  readSize = n ;//- 2*sizeof(long);
#ifdef MPIBUFF_MULTI_THREADED
    pthread_mutex_unlock(&mpibuff_mutex);
#endif
  return readSize;
} //mpibuff_pread

//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ mpibuff_read
ssize_t mpibuff_read(int fd, void *buf, size_t count){
  ssize_t n = mpibuff_pread(fd, buf, count, vectFiles[fd]->mpibuff_pos);
  if (n > 0)
    vectFiles[fd]->mpibuff_pos += n;
  return n;
}

//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ mpibuff_lseek64
off64_t mpibuff_lseek64(int fd, off64_t offset, int whence) throw() {
  _init_debug_mpibuff_lseek64(fd, offset, whence);
  switch (whence) {
    case SEEK_SET:
      vectFiles[fd]->mpibuff_pos = offset;
      break;
    case SEEK_CUR:
      vectFiles[fd]->mpibuff_pos += offset;
      break;
    case SEEK_END:
      vectFiles[fd]->mpibuff_pos = vectFiles[fd]->file_end_pos + offset;
      break;
  }
  if (vectFiles[fd]->mpibuff_pos < 0) {
    errno = EINVAL;
    return -1;
  }
  else
    return vectFiles[fd]->mpibuff_pos;
}
#endif //MPIBUFF_SUPPORTED

#endif //_IOBUFF_MPIBUFF_INTERFACE
