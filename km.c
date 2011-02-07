#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <time.h>
#include <fftw3.h>
#include <assert.h>
#include <string.h>
#include "km.h"
#include "km_FFT.h"
#include "sais.h"

/******************************************************************************/

   double *FFT_match_p = NULL;
   double *FFT_match_t = NULL;
   double *FFT_match_r = NULL;

/******************************************************************************/

// Count the frequencies of symbols in t. 
// We assume that A is the same size as the alphabet. 
void sp_km_count_symbols(const char *t, int n, int *A)
{
   int i=0;
   
   for (i=0; i<ALPHABET_SIZE; i++)
      A[i] = 0;
   
   for (i=0; i<n; i++)
      A[(int)t[i]] ++;
}

/******************************************************************************/

// Function to count the number of infrequent symbols.
int sp_km_count_infrequent_symbols(const char *A, int n, int t)
{
  int i;
  int count = 0;
  
  for (i=0; i<n; i++)
    if (A[i] >t) 
      count ++;
      
  return count;
}

/******************************************************************************/

// This assumes that the symbol occurs no more than |lookup| times.
// The idea is to store the indicies where the symbol occurs at each location
// in the lookup table.

void createLookup(int *lookup, char symbol, const char *pattern, int m, int l)
{
   int i;
   int x=0;
   
   
   for (i=0; i<m; i++)
   {
      if (pattern[i] == symbol)
      {
         lookup[x] = i;
         x++;
      }
   }
   
   for (i=x; i<l;i++)
   {
      lookup[i]=-1;
   } 
  
}

/******************************************************************************/
void markMatches                 (       int  *matches, 
                                   const char *text, 
                                         char  symbol, 
                                   const int  *lookup, 
                                         int   n, 
                                         int   l           )
{
   int i,j;
   
   for (i=0; i<n; i++)
   {  
      // Perform the marking.
      if (text[i] == symbol)
      {
         for (j=0; j<l; j++)
         {
            // TODO: can we do this better?
            if (lookup[j] == -1) break;
            if ( i - lookup[j] >= 0 )
               matches[i-lookup[j]] ++;
         }
      }      
   }
} 

/******************************************************************************/
int count_frequent_symbols(int threshold, const int *frequency_table, int n)
{
   int count = 0;
   
   for (int i=0; i<n; i++)
      if (frequency_table[i] > n) count ++;
   return count;
}

/******************************************************************************/

void sp_km_unbounded_kmismatch      ( const char   *text, 
                                     const char   *pattern, 
                                     int           n, 
                                     int           m,
                                     int           k,
		                               int          *matches,
		                               unsigned int  flags    )
{
   int i;

   // This will be a count of the occurences of symbols.
   int frequency_table[ALPHABET_SIZE];

   // Number of appearances required for a character to be classed 'frequent'.
   int FREQ_CHAR_THRESHOLD = sqrt(m);
   
   int MATCH_THRESHOLD     = 2 * sqrt(k);

   // Count the number of occurences of each symbol.
   sp_km_count_symbols(text, n, frequency_table);

   // CASE 1: pattern contains fewer than 2.sqrt{k} frequent symbols. //

   if (count_frequent_symbols(FREQ_CHAR_THRESHOLD, frequency_table, n) > MATCH_THRESHOLD)
   {
      // This will become a look up for each frequent character.
      int pattern_lookup[FREQ_CHAR_THRESHOLD];

      // Go through every symbol, looking for frequent symbols.
      // then perform FFT matching on each one.
      for (i=0; i<ALPHABET_SIZE; i++)
      {
         if (frequency_table[i] >= FREQ_CHAR_THRESHOLD)
         {
            match_with_FFT(matches,i, text, pattern,  n, m);
            printf("method 1\n");
         }
            else if (frequency_table[i] > 0)
         {
            printf("method 2\n");
            // Create lookup for this symbol.
            createLookup(pattern_lookup, i, pattern, m, FREQ_CHAR_THRESHOLD);
         
            // match this symbol in the text.
            markMatches(matches, text, i, pattern_lookup, n, FREQ_CHAR_THRESHOLD);
         }
      }
   }
      else
   {
   

   
   }  
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

// Create random pattern and text.

void randomStrings(char *text, char *pattern ,int n, int m)
{
   int i;
  
   
   for (i=0; i<n;i++)
      // random letter from a..i
      text[i] =  (char)(rand() % 26 + 97);
   
   text[n-1] = 0;   
   
   for (i=0;i<m;i++)
      pattern[i] = (char)(rand() % 26 + 97);
      
      
   pattern[m-1] = 0;

}

/******************************************************************************/


int main(int argc, char **argv)
{

   srand( time(NULL) );
   //-------------------------------------------------------------------------//
   // Testing parameters.
   //-------------------------------------------------------------------------//
   // Number of different test cases to try.
   int repeats = 1;
   
   // the length of the text and pattern.
   int m       = 10;
   int n       = 200;

   //-------------------------------------------------------------------------//
 
   int i,x;



   // The text and pattern strings.
   char t[n];
   char p[m];


   randomStrings(t,p, n, m);
   
   t[n-1] = 0;
   
   
   
   int SA[n+1];
   int LCP[n+1];
   
      printf("%s\n", t);
   
      sais(t, SA, LCP, n);
   
   exit(0);
   
   // The output arrays.
   int matches_FFT  [n-m+1];
   // int matches_naive[n-m+1];
   
   // Test the FFT Matching.  
   for (x=0; x < repeats; x++)
   {  
      // Generate random text and pattern of length n and m respectively.
      // These will consist only of the letters a and b.
      randomStrings(t, p, n,m);
      
      // zero the FFT matches array.
      for (i=0;i<n-m+1;i++)
         matches_FFT[i]=0; 
         
         
      //  printf("%s\n", t);
      //  printf("%s\n", p);
         
      sp_km_unbounded_kmismatch(t, p, n, m, 0, matches_FFT, 0);



      // Perform naive matching for testing.
      //    naive_matcher(t, p, n, m, matches_naive);
     /*
      // Check that the two outputs are the same.
      for (i=0;i<n-m+1;i++)
      {
      printf("%d ", matches_FFT[i]);
      }
      printf("\n");
      for (i=0;i<n-m+1;i++)
      {
      printf("%d ", matches_naive[i]);
      }
//         assert(matches_naive[i] == matches_FFT[i]);   

printf("\n");

*/
   }
        

   return 0;
}

/*******************************************************************************
* UNIT TESTING
*
* Compile with TEST to run unit tests.
*
*******************************************************************************/
#ifdef TEST
/******************************************************************************/

int main(int argc, char **argv)
{ 

   int status = 0;
   
   status += test_FFT_Matching();
   
   if (status==0)
   printf("All tests ran successfully.\n");
   
   return status;
  
  
  
}

/******************************************************************************/
#endif
/******************************************************************************/
