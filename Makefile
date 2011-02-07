CC    = gcc
flags = -std=c99 -Wall -pedantic -o3
libs  = -lm -lfftw3

all:
	$(CC) $(flags) -o km ./km.c ./km.h ./km_FFT.c $(libs)

