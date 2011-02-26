CC    = gcc
flags = -std=c99 -Wall -pedantic -o3 -g -pg
libs  = -lm -lfftw3

all:
	$(CC) $(flags) -o km ./stack.c ./RMQ_succinct.c ./RMQ_succinct.h ./sais.c ./stack.h ./sais.h ./km.c ./km.h ./km_FFT.c $(libs)
	$(CC) $(flags) -o make_input ./makeInput.c $(libs)
	

