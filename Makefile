CC    = gcc
flags = -std=c99 -Wall -pedantic -o3 -g -pg
libs  = -lm -lfftw3

all:
	$(CC) $(flags) -o km ./RMQ_succinct.c ./RMQ_succinct.h ./sais.c ./sais.h ./km.c ./km.h ./km_FFT.c $(libs)

