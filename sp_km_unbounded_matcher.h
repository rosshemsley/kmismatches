/*
 * unbounded.h
 *
 *  Created on: 29 Jul 2010
 *      Author: ben
 */

#ifndef UNBOUNDED_H_
#define UNBOUNDED_H_
#include <fftw3.h>
#include "sp_km.h"

#define NUM_CHARS 256
#define BITS_PER_CHAR CHAR_BIT
#define FREQUENT_CHAR -1
#define NOT_IN_PATTERN -2
#define INFREQUENT_CHAR -3

struct charAndPosition{
	char c;
	int index;
};

struct position{
	int index;
	short charType;
};

void done();

void sp_km_unbounded_kmismatch(char *text, char *pattern, int n, int m,int k,
		int *numMatches,struct SP_KM_MATCHING_POSITIONS *listOfMatches, unsigned int flags);

void maskTextAndPattern(char *text, char *pattern,int n, int m, double *maskedText,double *maskedPattern,char keep);

void computeNumMatchesWithFFT(double *t, double *maskedPattern,int n, int m,
		int transformSize,int *matches, fftw_plan *forward, fftw_plan *inverse);

void createLookupTable(char *pattern,int m,struct charAndPosition *sortedPattern,struct position *positionLookup,int frequentThreshold,int *numFrequentChars);

void findMismatchesWithoutFiltering(char *text, char *pattern, int n, int m,int k,
		int *numMatches,struct SP_KM_MATCHING_POSITIONS *listOfMatches,unsigned flags,
		struct charAndPosition *sortedPattern,struct position *positionLookup);

int compareCharAndPosition(const void *a, const void *b);

#endif /* UNBOUNDED_H_ */
