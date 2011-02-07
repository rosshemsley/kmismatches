/******************************************************************************/

#ifndef km_ross_H_
#define km_ross_H_
#include <fftw3.h>

/******************************************************************************/

// Maximum size of the alphabet. 
#define ALPHABET_SIZE (1 << CHAR_BIT)

/******************************************************************************/

typedef struct _MATCH
{
   int            k;
   int            index;
   struct _MATCH *next;
} MATCH;

#define ALPHABET_SIZE (1 << CHAR_BIT)

/******************************************************************************
* Functions.
*******************************************************************************/
int sp_km_unbounded_kmismatch      ( const char   *text, 
                                     const char   *pattern, 
                                     int           n, 
                                     int           m,
                                     int           k,
		                               int          *matches,
		                               unsigned int  flags    );
//----------------------------------------------------------------------------//
void sp_km_FFT_match_symbol        ( char         *A,
                                     char          s, 
                                     const char   *text, 
                                     const char   *pattern, 
                                     int           n, 
                                     int           m           );		                               
/******************************************************************************/
#endif
/******************************************************************************/
