CC    = gcc
flags = -std=c99 -Wall -pedantic -o3 -pg -g
libs  = -lm -lfftw3

harness: src/harness.c src/sp_km.c src/sp_km_unbounded_matcher.c src/stack.c src/RMQ_succinct.c src/RMQ_succinct.h src/sais.c src/stack.h src/sais.h src/km.c src/km.h src/km_FFT.c 
	cd src/ && $(CC) $(flags) -o ../harness harness.c sp_km.c sp_km_unbounded_matcher.c stack.c RMQ_succinct.c  sais.c km.c km_FFT.c $(libs)
	cd src/ && $(CC) $(flags) -o ../make_input makeInput.c $(libs)
	

