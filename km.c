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
void sp_km_count_symbols( const char *t, 
                                int   n, 
                                int  *A  )
{
   int i=0;
   
   for (i=0; i<ALPHABET_SIZE; i++)
      A[i] = 0;
   
   for (i=0; i<n; i++)
      A[(unsigned int)t[i]] ++;
}

/******************************************************************************/

// Function to count the number of infrequent symbols.
int sp_km_count_infrequent_symbols( const char *A, 
                                          int   n, 
                                          int   t  )
{
  int count = 0;
  
  for (int i=0; i<n; i++)
    if (A[i] >t) 
      count ++;
      
  return count;
}

/******************************************************************************/

// This assumes that the symbol occurs no more than |lookup| times.
// The idea is to store the indicies where the symbol occurs at each location
// in the lookup table.

void createLookup(       int   *lookup, 
                         char   symbol, 
                   const char  *pattern, 
                         int    m, 
                         int    l       )
{
   int x = 0;
   
   for (int i=0; i<m; i++)
   {
      if (pattern[i] == symbol)
      {
         lookup[x] = i;
         x++;
      }
   }

   lookup[x]=-1;
   
}

/******************************************************************************/
void markMatches(       int  *matches, 
                  const char *text, 
                        char  symbol, 
                  const int  *lookup, 
                        int   n, 
                        int   l       )
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
int count_frequent_symbols(       int  threshold, 
                            const int *frequency_table, 
                                  int  n                )
{
   int count = 0;
   
   for (int i=0; i<n; i++)
      if (frequency_table[i] > n) count ++;
   return count;
}


/******************************************************************************/
// We binary-search the SA until we find any string that starts with t.
/*
int start_p_triple(       char  t, 
                    const char *pattern, 
                    const int  *SA, 
                    const int  *LCP, 
                          int   n        )
{
   // For now we are lazy and do a linear search.
   for (int i=0; i<n;i++)
      if (pattern[SA[i]] == t) return i;
      
   return -1;
}
*/
/******************************************************************************/
/*
int extend_p_triple(       char  t, 
                     const char *pattern, 
                     const int  *SA, 
                     const int  *LCP, 
                           int   n)
{
   return 0;
}
*/

/******************************************************************************/
/*
void construct_p_representation( const char *text, 
                                 const char *pattern, 
                                 const int  *SA, 
                                 const int  *LCP, 
                                       int   n        )
{
   p_triple *p_representation = malloc(sizeof(p_triple) * n);

   int i=0;
   
   p_representation[0].i = 0;
   p_representation[0].j = start_p_triple(text[i], pattern, SA, LCP, n);
   p_representation[0].l = 1;
   
   while (i < n)
   {
      
   }  
}


*/
/******************************************************************************/

void constructSA( const char *text, 
                  const char *pattern, 
                        int   n, 
                        int   m       )
{
   // First: construct the SA/LCP for the pattern.
   
   int *SA  = malloc(sizeof(int) * n+1);
   int *LCP = malloc(sizeof(int) * n  );
   
   sais(text, SA, LCP, n);
  
 //  // Now, construct the p-representation for the text.  
   //construct_p_representation(text, pattern, SA, LCP, n);  
     

}

/******************************************************************************/

void abrahamson_kosaraju( const char *text, 
                          const char *pattern,
                                int   n,
                                int   m,
                                int  *matches  )
{
         
   // zero the FFT matches array.
   memset(matches, 0, sizeof(int)*(n-m+1));

   // This will be a count of the occurences of symbols.
   int frequency_table[ALPHABET_SIZE];
   
   printf("Counting Symbols\n");
   sp_km_count_symbols(pattern, m, frequency_table);

   printf("Performing Matching\n");
   
   // Number of appearances required for a character to be classed 'frequent'.
   int FREQ_CHAR_THRESHOLD =  sqrt(m * log(m)/log(2) );
   
   // printf("threshold: %d\n", FREQ_CHAR_THRESHOLD);
   
   // This will become a look up for each frequent character.
   // int pattern_lookup[FREQ_CHAR_THRESHOLD];

   // Go through every symbol, looking for frequent symbols.
   // then perform FFT matching on each one.
   for (int i=0; i<ALPHABET_SIZE; i++)
   {
   
      if (1) //frequency_table[i] >= FREQ_CHAR_THRESHOLD)
      {
         if (!(i >= 97 && i < 101)) continue;
         
         printf("method 1\n");      
         printf("%c\n", i);

         match_with_FFT(matches, i, text, pattern,  n, m);
      }
         else if (frequency_table[i] > 0)
      {
         printf("%c\n", i);
         printf("method 2\n");
         // Create lookup for this symbol.
         // createLookup(pattern_lookup, i, pattern, m, FREQ_CHAR_THRESHOLD);
      
         // match this symbol in the text.
         // markMatches(matches, text, i, pattern_lookup, n, FREQ_CHAR_THRESHOLD);
      }
   }
}                                


/******************************************************************************/
/*
void sp_km_unbounded_kmismatch      ( const    char *text, 
                                      const    char *pattern, 
                                               int   n, 
                                               int   m,
                                               int   k,
		                                         int  *matches,
		                                unsigned int   flags   )
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
       //     printf("method 1\n");
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
*/

/******************************************************************************/

void naive_matcher( const char *t, 
                    const char *p, 
                          int   n, 
                          int   m, 
                          int  *A  )
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

// Load the file at 'filename' into data, alloc'ing space required. 
// return the number of bytes loaded.

int loadData(       char ** data, 
              const char  *filename )
{
   FILE *fp;
   long n;

   if((fp = fopen(filename, "rb")) == NULL) {
      fprintf(stderr, "Cannot open file `%s': ", filename);
      exit(0);
   }
  
  
  fseek(fp, 0, SEEK_END);
  n = ftell(fp);
  rewind(fp);
  
  *data = malloc(sizeof(char) * n );
  
  if(fread(*data, sizeof(char), n, fp) != (size_t)n)
  {
     fprintf(stderr, "Could not load file.\n");
     exit(0);
  }
  
  return n;
}

/******************************************************************************/

int randDNA()
{
   int r = rand() % 4;
   
   if (r==1) return 'G';
   if (r==2) return 'T';
   if (r==3) return 'A';
   if (r==4) return 'C';

   return 0;
}



/******************************************************************************/

// Create random pattern and text.

void randomStrings( char *text, 
                    char *pattern,
                    int   n, 
                    int   m        )
{
   int i;
  
   
   for (i=0; i<n; i++)
      // random letter from a..z 
      text[i]    = (char)(rand() % 4 + 97);
   
   text[n] = 0;   
   
   for (i=0; i<m; i++)
      pattern[i] = (char)(rand() % 4 + 97);
   
   pattern[m] = 0;

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
   int m       = 5;
   int n       = 30;

   //-------------------------------------------------------------------------//
 
   int i,x;

   // The text and pattern strings.
   char *t = malloc(sizeof(char) * (n+1));
   char *p = malloc(sizeof(char) * (m+1));

   randomStrings(t,p, n, m);

   printf("%s\n", t);
   printf("%s\n", p);

   // n = loadData(&t, "./english.50MB");
   // printf("loaded, %d bytes\n", n);
   
   int * matches_FFT   = malloc(sizeof(int) * (n-m+1));
   int * matches_naive = malloc(sizeof(int) * (n-m+1));
   
   
   // Test the FFT Matching.  
   for (x=0; x < repeats; x++)
   {  
      // Generate random text and pattern of length n and m respectively.
      // These will consist only of the letters a and b.
      randomStrings(t, p, n, m);
         
      abrahamson_kosaraju(t, p, n, m, matches_FFT);
      // Perform naive matching for testing.
      naive_matcher(t, p, n, m, matches_naive);
  
      // Check that the two outputs are the same.
      for (i=0; i<n-m+1; i++)
         printf("%d ", matches_FFT[i]);
      
      printf("\n");
      
      for (i=0; i<n-m+1; i++)
         printf("%d ", matches_naive[i]);
      
      // for (i=0; i<n-m+1; i++)
      // assert(matches_naive[i] == matches_FFT[i]);   

      printf("\n");
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
