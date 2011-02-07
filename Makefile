CC    = gcc
flags = -std=c99 -Wall -pedantic -g
libs  = -lm -lfftw3

all:
	$(CC) $(flags) -o km ./sais.c ./sais.h ./km.c ./km.h ./km_FFT.c $(libs)

