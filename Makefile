CC    = gcc
flags = -std=c99 -Wall -pedantic -o3 -g -DNDEBUG
libs  = -lm -lfftw3

all:
	make harness -B
	make test_maker -B

harness: src/harness.c src/sp_km.c src/sp_km_unbounded_matcher.c src/stack.c src/RMQ_succinct.c src/RMQ_succinct.h src/loadTest.h src/esa.c src/sais.c src/stack.h src/sais.h src/km.c src/km.h src/km_FFT.c  src/breaks.c  src/breaks.h src/esa.h
	cd src/ && $(CC) $(flags) -o ../harness harness.c loadTest.c sp_km.c breaks.c sp_km_unbounded_matcher.c stack.c RMQ_succinct.c  sais.c esa.c km.c km_FFT.c $(libs)

	
test_maker: src/makeInput.c
	cd src/ && $(CC) $(flags) -o ../make_input makeInput.c $(libs)
clean:
	rm harness
	rm make_input
	rm gmon.out
