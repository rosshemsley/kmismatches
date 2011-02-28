/*
 * unbounded.c
 *
 *  Created on: 29 Jul 2010
 *      Author: ben
 *      Perform k-mismatches, ignoring the value of k (i.e. report Hamming
 *      distance for every alignment
 *      Time Complexity: O(n SQRT(m log m))
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <fftw3.h>
#include <string.h>
#include "sp_km_unbounded_matcher.h"

int THE_COUNT = 0;
/*
 * Find k-mismatches not bounded by k.
 * @param text A character array containing the text
 * @param pattern A character array containing the pattern
 * @param n the length of the text
 * @param m the length of the pattern
 * @param k The threshold, k. If k < 0 or >= m hamming distance at all
 * alignments will be found
 * @param results. An unsigned array to store the hamming distance in. Must be
 * of size n-m+1. If NULL, results are printed to stdout.
 * @param numMatches an unsigned pointer to a location to store the number
 * of k-mismatches found. Note that if k < 0 or >= m this is just n-m.
 */
void sp_km_unbounded_kmismatch(char *text, char *pattern, int n, int m,int k,
		int *numMatches,struct SP_KM_MATCHING_POSITIONS *listOfMatches, unsigned int flags){
	int i;
	if( k < 0){
		k = m+1;
	}
	*numMatches = 0;

	struct charAndPosition *sortedPattern = (struct charAndPosition *)malloc(sizeof(struct charAndPosition) * m);
	struct position *positionLookup = (struct position *)malloc(sizeof(struct position)*NUM_CHARS);

	//Set all characters to be non-existant in the pattern
	for(i=0;i<NUM_CHARS;i++){
		positionLookup[i].charType = NOT_IN_PATTERN;
	}
	int frequentThreshold = (int)sqrt((double)m*log2((double)m));
	
	printf("Frequent threshold: %d\n", frequentThreshold);
	//printf("Frequent threshold: %d\n",frequentThreshold);
//	frequentThreshold = 2;
	createLookupTable(pattern,m,sortedPattern,positionLookup,frequentThreshold,NULL);
	findMismatchesWithoutFiltering(text,pattern,n,m,k,numMatches,listOfMatches,flags,sortedPattern,positionLookup);
	free(sortedPattern);
	free(positionLookup);
}

/**
 * Comparator for qsort.
 */
int compareCharAndPosition(const void *a, const void *b){
	return ((struct charAndPosition *)a)->c - ((struct charAndPosition *)b)->c;
}

/*
 * Create a table to look up the list of infrequent chars in a charAndPosition array
 * @param pattern the pattern
 * @param m the length of the pattern
 * @param sortedPattern an array of charAndPosition elements in which to store
 * the sorted pattern
 * @param positionLookup an initialized array of size NUM_CHARS, with all charType values
 * set to NON_IN_PATTERN.
 * @param frequentThreshold if a character occurs >= frequentThreshold, it will be considered frequent
 * and have it's position lookup set to FREQUENT_CHAR
 * @param numFrequentChars -- a pointer to store the number of frequent characters found. May be NULL, in which case it won't be updated
 */
void createLookupTable(char *pattern,int m,struct charAndPosition *sortedPattern,struct position *positionLookup,int frequentThreshold, int *numFrequentChars){
	int i;
	for(i=0;i<m;i++){
		sortedPattern[i].c=pattern[i];
		sortedPattern[i].index=i;
	}
	qsort(sortedPattern,m,sizeof(struct charAndPosition),compareCharAndPosition);

	if(numFrequentChars != NULL){
		(*numFrequentChars) = 0;
	}
	char prev = sortedPattern[0].c;
	int numChars = 1;
	for(i=1;i<m;i++){
		if(prev==sortedPattern[i].c){
			numChars++;
		}else{
			positionLookup[(int)sortedPattern[i-1].c].index = i-numChars;
			if(numChars < frequentThreshold){
				positionLookup[(int)sortedPattern[i-1].c].charType = INFREQUENT_CHAR;
			}else{
				positionLookup[(int)sortedPattern[i-1].c].charType = FREQUENT_CHAR;
				if(numFrequentChars != NULL){
					(*numFrequentChars)++;
				}
			}
			numChars = 1;
			prev=sortedPattern[i].c;
		}
	}

	//Test last position
	positionLookup[(int)sortedPattern[i-1].c].index = i-numChars;
	if(numChars < frequentThreshold){
		positionLookup[(int)sortedPattern[i-1].c].charType = INFREQUENT_CHAR;
	}else{
		positionLookup[(int)sortedPattern[i-1].c].charType = FREQUENT_CHAR;
	}
}

 inline void maskTextAndPattern(char *text, char *pattern,int n, int m, double *maskedText,double *maskedPattern,char keep){
	int i;
	for(i=0;i<n;i++){
		maskedText[i] = (text[i]==keep) ? 1.0 : 0.0;
	}
	for(i=0;i<m;i++){
		maskedPattern[i] = (pattern[m-i-1]==keep) ? 1.0 : 0.0;
	}
}

inline void computeNumMatchesWithFFT(double *t, double *p,int n, int m,
		int transformSize,int *matches, fftw_plan *forward, fftw_plan *inverse){

   printf("USING FFT\n");

   printf("n: %d m: %d\n",n, m);

	//TODO: Should probably move these out of the function. No point allocating
	//and de-allocating the memory over and over for each character
	double *textSubString = (double *) fftw_malloc(sizeof(double) * transformSize);
	double *DFTofPattern = (double *) fftw_malloc(sizeof(double) * transformSize);
	double *DFTofText = (double *) fftw_malloc(sizeof(double) * transformSize);
	double *pTimesT = (double *) fftw_malloc(sizeof(double) * transformSize);
	double *pTimesT_PR = (double *) fftw_malloc(sizeof(double) * transformSize);
	int i;

	if(*forward==NULL){
		*forward = fftw_plan_r2r_1d(transformSize,p,DFTofPattern,FFTW_R2HC,FFTW_ESTIMATE);
		fftw_execute(*forward);
	}else{
		fftw_execute_r2r(*forward,p,DFTofPattern);
	}
	int start = 0;
	//Perform transforms on sub-strings of the text values
	while(start <= n-m){
		memcpy(textSubString,t+start,sizeof(double)*transformSize);
		fftw_execute_r2r(*forward,textSubString,DFTofText);

		/* Multiply the point representations*/

		pTimesT_PR[0] = DFTofPattern[0] * DFTofText[0];


      THE_COUNT++;
		if(transformSize % 2==0){

			for(i=1;i<transformSize/2;i++){
				pTimesT_PR[i] = DFTofPattern[i] * DFTofText[i] - DFTofPattern[transformSize-i] * DFTofText[transformSize-i];
				pTimesT_PR[transformSize-i] = DFTofPattern[transformSize-i] * DFTofText[i] + DFTofPattern[i] * DFTofText[transformSize-i];
			}
			pTimesT_PR[i] = DFTofPattern[i] * DFTofText[i];
		}else{

			for(i=1;i<=transformSize/2;i++){
				pTimesT_PR[i] = DFTofPattern[i] * DFTofText[i] - DFTofPattern[transformSize-i] * DFTofText[transformSize-i];
				pTimesT_PR[transformSize-i] = DFTofPattern[transformSize-i] * DFTofText[i] + DFTofPattern[i] * DFTofText[transformSize-i];
			}

		}

		/* Convert back to a coefficient representation */

		//On first iteration, need to create the inverse plan
		if(*inverse==NULL){
			*inverse = fftw_plan_r2r_1d(transformSize,pTimesT_PR,pTimesT,FFTW_HC2R,FFTW_ESTIMATE);
			fftw_execute(*inverse);
		}else{
			fftw_execute_r2r(*inverse,pTimesT_PR,pTimesT);
		}

		for(i=0;i<=transformSize-m && (i+start)<=n-m;i++){
			//printf("i+start: %d+%d=%d, value=%d\n",i,start,i+start,(int) ((pTimesT[m+i-1]/transformSize)+0.5));
			//plus 0.5 to allow rounding via truncation
			matches[i+start] += (int)((pTimesT[m+i-1]/transformSize)+0.5);
		}

		start +=transformSize-m +1;
	}
	fftw_free(textSubString);
	fftw_free(DFTofPattern);
	fftw_free(DFTofText);
	fftw_free(pTimesT);
	fftw_free(pTimesT_PR);
}

void done()
{
   printf("THE COUNT: %d\n", THE_COUNT);
}

/**
 * Find the mismatches without filtering. This is the complete method in the case where
 * the complexity is unbounded by k, but only applies in one case in the bounded method.
 * @TODO: Move this elsewhere; it doesn't apply just to unbounded. Perhaps the same for the above
 * @param text A character array containing the text
 * @param pattern A character array containing the pattern
 * @param n the length of the text
 * @param m the length of the pattern
 * @param k The threshold, k. If k < 0 or >= m hamming distance at all
 * alignments will be found
 * @param results. An unsigned array to store the hamming distance in. Must be
 * of size n-m+1. If NULL, results are printed to stdout.
 * @param numMatches an unsigned pointer to a location to store the number
 * of k-mismatches found. Note that if k < 0 or >= m this is just n-m.
 * @param sortedPattern. An array of charAndPosition structs, that is sorted but contain
 * indices of original location in the pattern
 * @param positionLookup. An int array of indexes into the sortedPattern saying where to find
 * the start of infrequent characters, or NOT_IN_PATTERN or FREQUENT_CHAR.
 */
void findMismatchesWithoutFiltering(char *text, char *pattern, int n, int m,int k,
		int *numMatches,struct SP_KM_MATCHING_POSITIONS *listOfMatches,unsigned flags,
		struct charAndPosition *sortedPattern,struct position *positionLookup){

	int i;
	int transformSize = 2*m;
	if(transformSize < 2048 && n > 4096){
		transformSize = 2048;
	}

	//Use double as we will compute Fourier Transform of this
	double *maskedPattern = (double *)fftw_malloc(sizeof(double)*transformSize);
	//Zero the top "half" of the masked pattern -- this will never be touched
	for(i=m;i<transformSize;i++){
		maskedPattern[i] = 0.0;
	}

	//Repeat for text, except we mask potential "overflow" of chunks past n
	int overflowed = n + (transformSize - m);
	double *maskedText = (double *)fftw_malloc(sizeof(double)*overflowed);
	for(i=n;i<overflowed;i++){
		maskedText[i] = 0.0;
	}

	//We will re-use plans across multiple masked patterns, so define these here
	fftw_plan forward = NULL;
	fftw_plan inverse =NULL;
	if(!fftw_import_system_wisdom()){
		printf("Failed to read system wisdom!\n");
	}

	//Create an array that we will use to add matches for each character at
	//each alignment
	int *matches = (int *) calloc(n-m+1,sizeof(int));
	int positionInPattern;
	short charType;
	int numFFTMatches = 0;
	double numInfrequentComparisons = 0.0;
	for(i=0;i<n;i++){
		charType = positionLookup[(int)text[i]].charType;
		positionInPattern = positionLookup[(int)text[i]].index;
		if(charType==FREQUENT_CHAR){
		

			maskTextAndPattern(text,pattern,n,m,maskedText,maskedPattern,text[i]);
						printf("Frequent char %c\n",text[i]);
			computeNumMatchesWithFFT(maskedText,maskedPattern,n,m,transformSize,matches,&forward,&inverse);
			numFFTMatches++;
			//This will compute matches for ALL of the same character in the text, so now mark this as not occuring in the pattern
			positionLookup[(int)text[i]].charType=NOT_IN_PATTERN;
		}else if(charType ==INFREQUENT_CHAR){
			//printf("Infrequent char %c\n",text[i]);
			//Loop though all of the up to O(threshold) characters that are infrequent
			while(positionInPattern < m && sortedPattern[positionInPattern].index <= i && sortedPattern[positionInPattern].c == text[i]){
//				printf("i-sortedPattern[positionInPattern].index+1: %d\n",i-sortedPattern[positionInPattern].index+1);
//				printf("i: %d\n",i);
//				printf("sortedPattern[positionInPattern].index: %d\n",sortedPattern[positionInPattern].index);
				if(i-sortedPattern[positionInPattern].index < n-m+1){
					matches[i-sortedPattern[positionInPattern].index]++;
				}
				positionInPattern++;
				numInfrequentComparisons++;
			}
		}
	}
	//printf("Calculated %d matches using Fourier Transforms\n",numFFTMatches);
	//printf("Performed %d (*c) comparisons for infrequent characters\n",(int)numInfrequentComparisons);
	//We now have matches, so we can compute (and output if necessary) mis-matches

	int hamDistance;
	for(i=0;i<n-m+1;i++){
		hamDistance = m - matches[i];
		if(hamDistance <= k){
			(*numMatches)++;
			if(listOfMatches!=NULL){
				sp_km_addToListOfMatches(listOfMatches,i,hamDistance);
			}else{
				printf("Hamming distance at text[%d]: %d\n",i,hamDistance);
			}
			if(flags & SP_KM_FIRST_MATCH_ONLY){
				return;
			}
		}
	}

   done();
	free(matches);
	fftw_free(maskedPattern);
	fftw_free(maskedText);
	fftw_destroy_plan(forward);
	fftw_destroy_plan(inverse);
	fftw_forget_wisdom();
	fftw_cleanup();

}

