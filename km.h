/******************************************************************************/

#ifndef km_ross_H_
#define km_ross_H_
#include <fftw3.h>

/******************************************************************************/

// Maximum size of the alphabet. 
#define ALPHABET_SIZE (1 << CHAR_BIT)

/******************************************************************************/

// This is a container for an Extended Suffix Array.
typedef struct _ESA
{
   int   n;
   int *SA;
   int *SAi;
   int *LCP;
   int *up;
   int *down;
   int *accross;
} ESA;

typedef struct _pTriple
{
   int j;
   int l;
}  pTriple;

#define ALPHABET_SIZE (1 << CHAR_BIT)

/******************************************************************************
* Functions.

*******************************************************************************/
/*
void sp_km_unbounded_kmismatch      ( const char   *text, 
                                     const char   *pattern, 
                                     int           n, 
                                     int           m,
                                     int           k,
		                               int          *matches,
		                               unsigned int  flags    );
		                               */
//----------------------------------------------------------------------------//
void sp_km_FFT_match_symbol        ( char         *A,
                                     char          s, 
                                     const char   *text, 
                                     const char   *pattern, 
                                     int           n, 
                                     int           m           );		
//----------------------------------------------------------------------------//                                     
void k_mismatches_case2(  const char *text, 
                          const char *pattern,
                          const int  *frequency_table,
                                int   k,
                                int   n,
                                int   m,
                                int  *matches          );          
                                
                                
void markMatches(                     const int*      lookup,        
                                 const int*      lookup_matrix,
                                       int       l,
                                       int*      matches,
                                 const char*     text,
                                       int       n,       
                                       int       m                             );                                
                                
                                
//int extendInterval(pTriple *P, const char *text, const char *pattern, int n, int m, const ESA *esa);
                                                                                       
/******************************************************************************/
#endif
/******************************************************************************/
