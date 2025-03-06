all:
	g++ -g -O0 -shared -ldl -fPIC -fpermissive -DIOBUFF_DEBUG iobuff.cpp -DMPIBUFF_SUPPORTED -olibiobuff_debug_mpism.so;
	g++ -g -O0 -shared -ldl -fPIC -fpermissive -DIOBUFF_DEBUG iobuff.cpp -olibiobuff_debug.so;
	g++ -O3 -shared -ldl -fPIC -DMPIBUFF_SUPPORTED iobuff.cpp -olibiobuff_mpism.so ;
	g++ -O3 -shared -ldl -fPIC iobuff.cpp -olibiobuff.so ;
	g++ -O3 -shared -ldl -fPIC -DIOBUFF_ALLOW_DEV iobuff.cpp -olibiobuff_dev_shm.so ;
#	g++ -g -O0 -shared -ldl -fPIC -fpermissive -DIOBUFF_DEBUG_ENTRY iobuff.cpp -olibiobuff_entry.so
