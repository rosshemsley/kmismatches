CC    = gcc
flags = -Wall -pedantic -std=c99 -03
libs  = -lm -lfftw3

all:
	$(CC) -o km ./km.c ./km.h ./km_FFT.c $(libs)

