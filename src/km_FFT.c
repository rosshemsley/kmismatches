#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <fftw3.h>
#include <assert.h>
#include <string.h>
#include <time.h>

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

double * x       = NULL;
double *_p       = NULL;
double *_t       = NULL;
double *_r       = NULL;
double *t_masked = NULL;

// FFT Plans which we can re-use.
fftw_plan plan_pattern_FFT = NULL;
fftw_plan plan_text_FFT    = NULL;
fftw_plan plan_FFT_inverse = NULL;

/******************************************************************************/

// Mask the text into a buffer of doubles.
// We assume that the text is the same length as the buffer for now.

static inline void maskText(    double*    r, 
                                const char       symbol, 
                                const char*      text, 
                                const int        n                             )
{
   for (int i= 0; i<n; i++)
      r[i] = (text[i] == symbol) ? 1.0 : 0.0;

}

/******************************************************************************/

// Mask the pattern and put it into a buffer of type double. 
// We assume that the pattern is half the length of the buffer.
// To perform matching, we need the pattern to be reversed.

static void maskPattern(               double*   r, 
                                 const char      symbol, 
                                 const char*     pattern, 
                                 const int       m, 
                                 const int       N                             )

{
   // Copy the pattern into the bottom of the buffer.
   // NB. We need to reverse the pattern.
   for (int i=0; i<m; i++)
      r[m-i-1] = (pattern[i] == symbol) ? 1.0 : 0.0;

   // Fill the top half of the buffer with zeros.
   for (int i=m; i<N; i++)
      r[i] = 0;

}

/******************************************************************************/

// This will multiply two vectors in half-complex notation, when they are of
// length N and N is even.

static inline void multiply_half_complex(       double  *r, 
                                          const double  *a, 
                                          const double  *b, 
                                          const int      N  )
{
   // TODO: CHECK THIS IS RIGHT...
   // The first entry is real.
   r[0] = a[0] * b[0];

   // TODO: make this faster using clever complex multiplcation.
   for (int i=1; i<=N/2; i++)
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
   
   int transformSize = 2*m;
   
	if(transformSize < 2048 && n > 4096){
	   transformSize = 2048;
   	}
	
	if(!fftw_import_system_wisdom()){
      printf("Failed to read system wisdom!\n");
	}

   
   
   if (N != transformSize)
   {
      N = transformSize;
      
       x       =  fftw_malloc(sizeof(double) * N);
      _p       =  fftw_malloc(sizeof(double) * N);
      _t       =  fftw_malloc(sizeof(double) * N);
      _r       =  fftw_malloc(sizeof(double) * N);
      t_masked =  fftw_malloc(sizeof(double) * (n + N ));
      
      
    
         
         

      // Any overflow over the end of the text should be set to zero.
      for (int i=n; i < n + N ; i++)
         t_masked[i] = 0.0;
    }       
    
    // Mask the text into an array which we can copy from later.
    maskText(t_masked, symbol, text, n); 
    maskPattern(x, symbol, pattern, m, N);
   
   //printf("FFT'ing pattern\n");
   // Compute the FFT of x (the pattern) and put it in _p.
   if (plan_pattern_FFT == NULL)
   {
      plan_pattern_FFT = fftw_plan_r2r_1d(N, x, _p, FFTW_R2HC, FFTW_ESTIMATE);
      fftw_execute(plan_pattern_FFT);
   } else {
      fftw_execute_r2r(plan_pattern_FFT, x , _p);
   }
   
   
   /* MATCHING */
   
   for (int i=0; i<n-m; i += N-m+1)
   {
      // Copy the text into the buffer x.
      memcpy(x, t_masked+i, sizeof(double)*N);

      // FFT x (this block of text) and store it in _t.
     // if (plan_text_FFT == NULL)
    //  {
     //    plan_text_FFT = fftw_plan_r2r_1d(N, x, _t, FFTW_R2HC, FFTW_ESTIMATE);
     //    fftw_execute(plan_text_FFT);
    //  } else {
         fftw_execute_r2r(plan_pattern_FFT, x, _t);
   //   }
      
	   // Point-wise multiply the two vectors _p and _t.
      multiply_half_complex(_r, _t, _p, N);



      // FFTW invert _r and put it into x.
      if (plan_FFT_inverse == NULL)
      {
         plan_FFT_inverse = fftw_plan_r2r_1d(N, _r,  x, FFTW_HC2R, FFTW_ESTIMATE);
         fftw_execute(plan_FFT_inverse);
	   } else {
	      fftw_execute_r2r(plan_FFT_inverse, _r, x); 
	   } 
	   



   
      // x now contains the matches.
      for (int j=0; j <= N-m; j+=1)
      {
         if (i+j > n-m+1) break;
         
         // Add to matches array.
         // The matches we are interested in start half way into the array.

        
         matches[i+j] += (int)(x[j+m-1]/N + 0.5);
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

void randomStrings(char *text, char *pattern ,int n, int m)
{	
   int i;	

   for (i=0; i<n;i++)
      // random letter from a..b
      text[i] =  (char)(rand() % 2 + 97);
   for (i=0;i<m;i++)
      pattern[i] = (char)(rand() % 2 + 97);
}

/******************************************************************************/

void naive_matcher(const char *t, const char *p, int n, int m, int *A)
{
   int i,j;
   for (i=0;i<n-m+1;i++)
   {
      int matches=0;
      for (j=0;j<m;j++)
      {
         if (i+j > n) continue;
       
         if (t[i+j] == p[j])
            matches ++;
    }
    A[i] = matches;
   }
}

/******************************************************************************/

// Unit testing. Will return 0 if all tests succeeded.
int test_FFT_Matching()
{ 

   printf("Testing FFT Matching\n");
   srand( time(NULL) );

   //-------------------------------------------------------------------------//
   // Testing parameters.
   //-------------------------------------------------------------------------//
   // Number of different test cases to try.
   int repeats = 3e2;
   
   // the length of the text and pattern.
   int m       = 1e2;
   int n       = 1e4;

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
