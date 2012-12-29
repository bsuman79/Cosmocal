# mac
CC = g++

CFLAGS = -I/usr/local/include
CLIBS = -L/usr/local/lib -lgsl -lgslcblas  -lm

.cpp.o :
	$(CC) $(CFLAGS)  -c $<

cosmo: cosmo_driver.o
	$(CC) $(CFLAGS) -o cosmo_driver.exe cosmo_driver.o $(CLIBS)

clean : 
	/bin/rm -f *.o 

