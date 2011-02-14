CC    = g++
flags = -std=c99 -Wall -pedantic -o2 -pg -g
libs  = -lm -lfftw3

all:
	$(CC) $(flags) -o km ./sais.c ./sais.h ./km.c ./km.h ./km_FFT.c $(libs)

