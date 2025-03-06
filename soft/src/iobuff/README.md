IOBUFF : A ld_preloaded library to buffer read/write file operations and reduce
         the total number of IOPS.

Compiling :
  make

Using :
  LD_PRELOAD=/path/to/libiobuffcpp.so /path/to/executable


Usage notes :

  IOBUFF does not allow a file to be opened in Write or Read-Write mode at the same
  time that it is opened in any mode. Default behavior is to crash with an error message.
  A file can be opened any number of times in Read-only mode.
  
  IOBUFF does not detect if a file open in Write / Read-Write mode is also opened by a 
  second process. This case can lead to data coruption and it is to be avoided. 
  
  The library can be parametrized by changing the values of key environment variables:

  IOBUFF_EXEC_LIST (default value = "ALL") : Allow to specify which executables can be
    buffered through IOBUFF and which should not (using the negation !) .
    ex.:
    To allow only ./foo , /bin/bar and bat   :
      export IOBUFF_EXEC_LIST=./foo:/bin/bar:\`which bat\`

    To disallow only ./foo   :
      export IOBUFF_EXEC_LIST=all:\!./foo
    or
      export IOBUFF_EXEC_LIST="all:!./foo"
    
  IOBUFF_BUFFER_SIZE (default = 400*1024) : Specifies the number of bytes the library
    should prefetch (read) or buffer (write)
