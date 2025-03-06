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
#ifndef _IOBUFF_DEBUG_FORK
#define _IOBUFF_DEBUG_FORK


#ifdef IOBUFF_DEBUG_FORK
pid_t (*real_fork)(void)=NULL;
int (*real_execl)(const char *, const char *)=NULL;
int (*real_execlp)(const char *, const char *)=NULL;
int (*real_execle)(const char *, const char *, char **const )=NULL;
int (*real_execv)(const char *, char **const )=NULL;
int (*real_execvp)(const char *, char **const )=NULL;
int (*real_execvpe)(const char *, char **const, char **const )=NULL; 
#endif
///////////////////////////////////////////////////////////////////////////////////////////////////
void _init_fork_func(){
#ifdef IOBUFF_DEBUG_FORK
  real_fork  = (pid_t   (*)(void))                     dlsym(RTLD_NEXT, "fork");
  real_execl = (int     (*)(const char *, const char *))dlsym(RTLD_NEXT, "execl");
  real_execlp= (int     (*)(const char *, const char *))dlsym(RTLD_NEXT, "execlp");
  real_execle= (int     (*)(const char *, const char *, char **const ))dlsym(RTLD_NEXT, "execle");
  real_execv = (int     (*)(const char *, char **const ))dlsym(RTLD_NEXT, "execv");
  real_execvp= (int     (*)(const char *, char **const ))dlsym(RTLD_NEXT, "execvp");
  real_execvpe=(int     (*)(const char *, char **const, char **const ))dlsym(RTLD_NEXT, "execvpe");
#endif

}
 
#ifdef IOBUFF_DEBUG_FORK
///////////////////////////////////////////////////////////////////////////////////////////////////
  pid_t fork(void) {
#ifdef IOBUFF_DEBUG_ENTRY
  fprintf(debugOut, "IOBUFF Entry: fork\n");
#endif
    spin_wait_init();
    pid_t ret = real_fork();
    fprintf(debugOut,"IOBUFF: fork pid_t %d\n", ret);
    return ret;
  }//fork
///////////////////////////////////////////////////////////////////////////////////////////////////
  int execl(const char *path, const char *arg, ...){
#ifdef IOBUFF_DEBUG_ENTRY
  fprintf(debugOut, "IOBUFF Entry: execl\n");
#endif
    spin_wait_init();
#define INITIAL_ARGV_MAX 1024
    size_t argv_max = INITIAL_ARGV_MAX;
    const char *initial_argv[INITIAL_ARGV_MAX];
    const char **argv = initial_argv;
    va_list args;
    argv[0] = arg;
    pid_t lp2 = (pid_t)syscall(__NR_gettid);
    fprintf(debugOut,"IOBUFF: exec l pid %d\n", lp2);
    va_start (args, arg);
    unsigned int i = 0;
    while (argv[i++] != NULL){
      if (i == argv_max) {
        argv_max *= 2;
        const char **nptr = (const char **)realloc (argv == initial_argv ? NULL : 
                             argv, argv_max * sizeof (const char *));
        if (nptr == NULL){
          if (argv != initial_argv)
            free (argv);
          return -1;
        }
        if (argv == initial_argv)
        /* We have to copy the already filled-in data ourselves.  */
          memcpy (nptr, argv, i * sizeof (const char *)); 
        argv = nptr;
      }
      argv[i] = va_arg (args, const char *);
    }
    va_end (args);
    int ret = execve (path, (char *const *) argv, __environ);
    if (argv != initial_argv)
      free (argv);
    return ret;
  }//execl

///////////////////////////////////////////////////////////////////////////////////////////////////
  int execlp(const char *file, const char *arg, ...){
#ifdef IOBUFF_DEBUG_ENTRY
  fprintf(debugOut, "IOBUFF Entry: execlp\n");
#endif
    spin_wait_init();
#define INITIAL_ARGV_MAX 1024
    size_t argv_max = INITIAL_ARGV_MAX;
    const char *initial_argv[INITIAL_ARGV_MAX];
    const char **argv = initial_argv;
    va_list args;
    argv[0] = arg;
    pid_t lp2 = (pid_t)syscall(__NR_gettid);
    fprintf(debugOut,"IOBUFF: exec lp pid %d\n", lp2);
    va_start (args, arg);
    unsigned int i = 0;
    while (argv[i++] != NULL) {
      if (i == argv_max) {
        argv_max *= 2;
        const char **nptr = (const char **)realloc (argv == initial_argv ? NULL : 
                               argv, argv_max * sizeof (const char *));
        if (nptr == NULL) {
          if (argv != initial_argv)
            free (argv);
          return -1;
        }
        if (argv == initial_argv)
        /* We have to copy the already filled-in data ourselves.  */
          memcpy (nptr, argv, i * sizeof (const char *));
        argv = nptr;
      }
      argv[i] = va_arg (args, const char *);
    }
    va_end (args);
    int ret = execvp (file, (char *const *) argv);
    if (argv != initial_argv)
      free (argv);
    return ret;    
  }//execlp
///////////////////////////////////////////////////////////////////////////////////////////////////
  int execle(const char *path, const char *arg, ...){
#ifdef IOBUFF_DEBUG_ENTRY
  fprintf(debugOut, "IOBUFF Entry: execle\n");
#endif
    spin_wait_init();
#define INITIAL_ARGV_MAX 1024
    size_t argv_max = INITIAL_ARGV_MAX;
    const char *initial_argv[INITIAL_ARGV_MAX];
    const char **argv = initial_argv;
    va_list args;
    argv[0] = arg;
    pid_t lp2 = (pid_t)syscall(__NR_gettid);
    fprintf(debugOut,"IOBUFF: exec le pid %d\n", lp2);
    va_start (args, arg);
    unsigned int i = 0;
    while (argv[i++] != NULL){
      if (i == argv_max) {
        argv_max *= 2;
        const char **nptr = (const char **)realloc (argv == initial_argv ? NULL : 
                              argv, argv_max * sizeof (const char *));
        if (nptr == NULL) {
          if (argv != initial_argv)
            free (argv);
          return -1;
        }
        if (argv == initial_argv)
        /* We have to copy the already filled-in data ourselves.  */
          memcpy (nptr, argv, i * sizeof (const char *));
        argv = nptr;
      }
      argv[i] = va_arg (args, const char *);
    }
    const char *const *envp = va_arg (args, const char *const *);
    va_end (args);
    int ret = execve (path, (char *const *) argv, (char *const *) envp);
    if (argv != initial_argv)
    free (argv);  
  }//execle
///////////////////////////////////////////////////////////////////////////////////////////////////
  int execv(const char *path, char *const argv[]){
#ifdef IOBUFF_DEBUG_ENTRY
  fprintf(debugOut, "IOBUFF Entry: execv\n");
#endif
    spin_wait_init();
    pid_t lp2 = (pid_t)syscall(__NR_gettid);
    fprintf(debugOut,"IOBUFF: exec v pid %d\n", lp2);
    return execve (path, argv, __environ);
  }//execv
///////////////////////////////////////////////////////////////////////////////////////////////////
  int execvp(const char *file, char *const argv[]){
#ifdef IOBUFF_DEBUG_ENTRY
  fprintf(debugOut, "IOBUFF Entry: execvp\n");
#endif
    spin_wait_init();
    pid_t lp2 = (pid_t)syscall(__NR_gettid);
    fprintf(debugOut,"IOBUFF: execvp pid %d\n", lp2);
    return execvpe (file, argv, __environ);
  }//execvp
///////////////////////////////////////////////////////////////////////////////////////////////////
  /* The file is accessible but it is not an executable file.  Invoke
      the shell to interpret it as a script.  */
  static void scripts_argvb (const char *file, char *const argv[], int argc, char **new_argv);
  static void scripts_argvb (const char *file, char *const argv[], int argc, char **new_argv) {
   /* Construct an argument list for the shell.  */
   new_argv[0] = (char *) _PATH_BSHELL;
   new_argv[1] = (char *) file;
   while (argc > 1) {
     new_argv[argc] = argv[argc - 1];
     --argc;
   }
  }//script_argvb
///////////////////////////////////////////////////////////////////////////////////////////////////
  int execvpe(const char *file, char *const argv[], char *const envp[]){
#ifdef IOBUFF_DEBUG_ENTRY
  fprintf(debugOut, "IOBUFF Entry: execvpe\n");
#endif
    spin_wait_init();
    pid_t lp2 = (pid_t)syscall(__NR_gettid);
    fprintf(debugOut,"IOBUFF: exec vpe pid %d\n", lp2);
    if (*file == '\0') {
      /* We check the simple case first. */
      errno = ENOENT;
      return -1;
    } 
    if (strchr (file, '/') != NULL){
      /* Don't search when it contains a slash.  */
      execve (file, argv, envp);
      
      if (errno == ENOEXEC){
        /* Count the arguments.  */
        int argc = 0;
        while (argv[argc++]) ;
        
        size_t len = (argc + 1) * sizeof (char *);
        char **script_argv;
        void *ptr = NULL; 
        if (len<= __MAX_ALLOCA_CUTOFF)
          script_argv = alloca (len);
        else
          ptr = malloc (len);
          script_argv = (char**)ptr;
        if (script_argv != NULL) {
           scripts_argvb (file, argv, argc, script_argv);
           execve (script_argv[0], script_argv, envp); 
           free (ptr);
         }
      }
    }
    else {
       size_t pathlen;
       size_t alloclen = 0;
       char *path = getenv ("PATH");
       if (path == NULL) {
         pathlen = confstr (_CS_PATH, (char *) NULL, 0);
         alloclen = pathlen + 1;
       }
       else
         pathlen = strlen (path);
    
       size_t len = strlen (file) + 1;
       alloclen += pathlen + len + 1;
    
       char *name;
       char *path_malloc = NULL;
       if (alloclen<= __MAX_ALLOCA_CUTOFF)
         name = alloca (alloclen);
       else{
         path_malloc = name = (char *) malloc (alloclen);
         if (name == NULL)
           return -1;
       }
    
       if (path == NULL){
         /* There is no `PATH' in the environment.
            The default search path is the current directory
            followed by the path `confstr' returns for `_CS_PATH'.  */
         path = name + pathlen + len + 1;
         path[0] = ':';
         (void) confstr (_CS_PATH, path + 1, pathlen);
       }
    
       /* Copy the file name at the top.  */
       name = (char *) memcpy (name + pathlen + 1, file, len);
       /* And add the slash.  */
       *--name = '/';
    
       char **script_argv = NULL;
       void *script_argv_malloc = NULL;
       bool got_eacces = false;
       char *p = path;
       do {
         char *startp;
      
         path = p;
         p = strchrnul (path, ':');
      
         if (p == path)
           /* Two adjacent colons, or a colon at the beginning or the end
              of `PATH' means to search the current directory.  */
           startp = name + 1;
         else
           startp = (char *) memcpy (name - (p - path), path, p - path);
      
         /* Try to execute this name.  If it works, execve will not return. */
         execve (startp, argv, envp);
      
         if (errno == ENOEXEC){
           if (script_argv == NULL){
             /* Count the arguments.  */
             int argc = 0;
             while (argv[argc++])
               ;
             size_t arglen = (argc + 1) * sizeof (char *);
             if ((alloclen + arglen) <= __MAX_ALLOCA_CUTOFF)
               script_argv = alloca (arglen);
             else{
               script_argv_malloc = malloc (arglen);
               script_argv = (char **)script_argv_malloc;
             }
             if (script_argv == NULL){
               /* A possible EACCES error is not as important as
                  the ENOMEM.  */
               got_eacces = false;
               break;
             }
             scripts_argvb (startp, argv, argc, script_argv);
           }
           execve (script_argv[0], script_argv, envp);
         }
         switch (errno){
           case EACCES:
             got_eacces = true;
           case ENOENT:
           case ESTALE:
           case ENOTDIR:
           case ENODEV:
           case ETIMEDOUT:
            break;
           default:
             return -1;
           }
       } while (*p++ != '\0');
       /* We tried every element and none of them worked.  */
       if (got_eacces)
         /* At least one failure was due to permissions, so report that
            error.  */
          errno = EACCES;
    
       free (script_argv_malloc);
       free (path_malloc);
     }
    /* Return the error from the last attempt (probably ENOENT).  */
    return -1;
  }
#endif


#endif //_IOBUFF_DEBUG_FORK
