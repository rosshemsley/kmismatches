#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <fftw3.h>
#include <assert.h>
#include <string.h>

// #define TEST
/*******************************************************************************
*
*
*
********************************************************************************
* GLOBAL VARIABLES. 
*
* These are used for caching purposes. In particular, we store caches of the 
* input length so that we can re-use plans.
*******************************************************************************/

// Start off N as an impossible value.
int      N = -1;

double * x;
double *_p;
double *_t;
double *_r;

// FFT Plans which we can re-use.
fftw_plan plan_pattern_FFT;
fftw_plan plan_text_FFT;
fftw_plan plan_FFT_inverse;

/******************************************************************************/

// Mask the text into a buffer of doubles.
// We assume that the text is the same length as the buffer for now.

void maskText(double *r, char symbol, const char *text, int n, int N)
{
   int i;
   
   // Assume for now that the text is always the same slength as the buffer.
   for (i=0; i<n && i<N; i++)
      r[i] = (text[i] == symbol) ? 1 : 0;
      
   for (i=n; i<N; i++)
      r[i] = 0;

}

/******************************************************************************/

// Mask the pattern and put it into a buffer of type double. 
// We assume that the pattern is half the length of the buffer.
// To perform matching, we need the pattern to be reversed.

void maskPattern(double *r, char symbol, const char *pattern, int N)
{
   int i;
   
   // Copy the pattern into the bottom half of the buffer.
   for (i=0; i<N/2; i++)
      r[N/2-i-1] = (pattern[i] == symbol) ? 1 : 0;
      
   // Fill the top half of the buffer with zeros.
   for (i=N/2; i<N; i++)
      r[i] = 0;

}

/******************************************************************************/

// This will multiply two vectors in half-complex notation, when they are of
// length N and N is even.

void multiply_half_complex(       double  *r, 
                            const double  *a, 
                            const double  *b, 
                                  int      N  )
{
   int i;
   
   // The first entry is real.
   r[0] = a[0] * b[0];

   // TODO: make this faster using clever complex multiplcation.
   for (i=1; i<N/2; i++)
   {
     // Multiply the FFTW half-complex vectors.
     r[i]   = a[i]   * b[i] - a[N-i] * b[N-i];
     r[N-i] = a[N-i] * b[i] + a[i]   * b[N-i];
   }

   // If the length is even, the middle entry is also real.
   if (N % 2 == 0)
     r[N/2] = a[N/2] * b[N/2];
}

/******************************************************************************/

void printArr(int n, double *d, int r)
{ 
  int i;
  for (i=0; i<n; i++)
  {
   if (!r)
    printf("%2.0f", d[i]);  
   else  
      if (i<n/2)
         printf("%2.0f", d[n/2-i-1]);  
      else 
         printf(" 0");
  }
  printf("\n");
}


/******************************************************************************/

// Use the FFT to find all the matches of this particular symbol.
// We assume that the matches array is allocated to length n-m+1 in order
// to store all the matches.

void match_with_FFT(        int  *matches, 
                            char  symbol,
                      const char *text, 
                      const char *pattern, 
                            int   n,
                            int   m        )
{
   int i,j;
   
   // ---- PRECOMPUTATION ---- //
   
   // Store these vectors as globals, and re-use them to save re-allocating
   // memory on every run.
   if (2*m != N)
   {
      // We perform all our FFT manipulations with vectors twice as long as 
      // the length of the pattern.
      N = 2*m;
      
      // A temporary variable to store intermediate vectors.
      x = realloc( x, sizeof(double) * N );

      // the FFT representation of the pattern, text, and pattern*text.
      _p = realloc( _p, sizeof(double) * N );
      _t = realloc( _t, sizeof(double) * N );
      _r = realloc( _r, sizeof(double) * N );
      
      fftw_destroy_plan( plan_pattern_FFT );    
      fftw_destroy_plan( plan_text_FFT    );    
      fftw_destroy_plan( plan_FFT_inverse );    
      
      plan_pattern_FFT = fftw_plan_r2r_1d(N,  x, _p, FFTW_R2HC, FFTW_ESTIMATE);
      plan_text_FFT    = fftw_plan_r2r_1d(N,  x, _t, FFTW_R2HC, FFTW_ESTIMATE);
      plan_FFT_inverse = fftw_plan_r2r_1d(N, _r,  x, FFTW_HC2R, FFTW_ESTIMATE);
   }
   
   // ---- END OF PRECOMPUTATION ---- //

   // This will copy the pattern into the buffer x.
   maskPattern(x, symbol, pattern, N);
   
   // Compute the FFT of x (the pattern) and put it in _p.
   fftw_execute(plan_pattern_FFT);

   for (i=0; i<n; i += m)
   {
      // Copy the text into the buffer x.
      maskText(x, symbol, text+i, n-i,  N);

      // FFT x (this block of text) and store it in _t.
	   fftw_execute(plan_text_FFT);

	   // Point-wise multiply the two vectors _p and _t.
      multiply_half_complex(_r, _t, _p, N);
      
      // FFTW invert _r and put it into x.
	   fftw_execute(plan_FFT_inverse);

      // x now contains the matches.
      for (j=0; j<N/2; j++)
      {
         if (i+j > n-m+1) break;
         
         // Add to matches array.
         // The matches we are interested in start half way into the array.
         matches[i+j] += (int)(x[j+N/2-1]/N + 0.5);
      }
   }
} 

/*******************************************************************************
* UNIT TESTING
*
* Compile with TEST to run unit tests.
*
*******************************************************************************/
#ifdef TEST
/******************************************************************************/

// Unit testing. Will return 0 if all tests succeeded.
int test_FFT_Matching()
{ 

   //-------------------------------------------------------------------------//
   // Testing parameters.
   //-------------------------------------------------------------------------//
   // Number of different test cases to try.
   int repeats = 1e3;
   
   // the length of the text and pattern.
   int m       = 19;
   int n       = 278;

   //-------------------------------------------------------------------------//
 
   int i,x;


   // The text and pattern strings.
   char t[n+1];
   char p[m+1];
   p[m] = '\0';
   t[n] = '\0';
   
   // The output arrays.
   int matches_FFT  [n-m+1];
   int matches_naive[n-m+1];
   
   // Test the FFT Matching.  
   for (x=0; x < repeats; x++)
   {  
      // Generate random text and pattern of length n and m respectively.
      // These will consist only of the letters a and b.
      randomStrings(t,p, n,m);
      
      // zero the FFT matches array.
      for (i=0;i<n-m+1;i++)
         matches_FFT[i]=0; 
         
      // Add the matches for 'a', and then 'b' using the FFT.
      match_with_FFT(matches_FFT,'a', t, p, n, m);
      match_with_FFT(matches_FFT,'b', t, p, n, m);

      // Perform naive matching for testing.
      naive_matcher(t, p, n, m, matches_naive);
   
      // Check that the two outputs are the same.
      for (i=0;i<n-m+1;i++)
         assert(matches_naive[i] == matches_FFT[i]);    
   }
        
   return 0;

}

/******************************************************************************/
#endif
/******************************************************************************/
