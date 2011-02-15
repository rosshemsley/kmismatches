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
#include "RMQ_succinct.h"

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

int count_frequent_symbols(       
                            const int *frequency_table, 
                                  int  threshold, 
                                  int  n                )
{
   int count = 0;
   
   for (int i=0; i<n; i++)
      if (frequency_table[i] > n) count ++;
   return count;
}

/******************************************************************************/

// This will take the current position in the SA, x and the current
// pTriple p and the next letter in the text. It will then loop through the 
// following sufficies until the first suffix which has the correct next letter
// If this is found, we return the suffix and extend l by one. If the LCP becomes
// to small (i.e. the preceeding chars don't match) then we return -1, 
// and must start a new pTriple.

int extend     (          char   t, 
                           int   l,
                           int  *x,
                     const char *pattern, 
                     const  ESA *esa,
                           int   n,     
                           int   m         )
{

   //printf("Extending\n");   
   // Go through the Suffix Array until LCP[i] < l or pattern[SA[i]+l] = t.
   
   // Check to see whether or not we can extend the current suffix.   
   if ( (esa->SA[*x]+l < m) && (pattern[esa->SA[*x] + l] == t) )
   {
      // printf("Succeeded with current\n");
      // printf("'%s'\n", pattern + SA[*x]); 
      return 1;
   }
   
   // Look for other possible sufficies which may match.
   for (int i = *x+1; i<n; i++)
   {
      // There are no sufficies with this prefix.
      if (esa->LCP[i] < l ) {
         //printf("Ran out of values\n");
         return 0;
      }
      
      // Check this character match.
      if ( (esa->SA[*x]+l < m) && (pattern[esa->SA[i] + l] == t) )
      {
         //printf("Changing to '%s'\n", pattern + SA[i]); 
         *x = i;
         return 1;
      }
   }

   return -1;
}

/******************************************************************************/

int findStart(char c, const char *pattern, const int *SA, int m)
{

   // For now we are lazy and do a linear search.
   for (int i=0; i<m; i++)
      if (pattern[SA[i]] == c) {
         return i;
      }
   return -1;
}

/******************************************************************************/


void construct_pRepresentation(       pTriple   *P,
                                const char      *text, 
                                const char      *pattern, 
                                const ESA       *esa,
                                      int        n,
                                      int        m              )
{

   int p=0;
   int i=0;
   int x=0;
   int l=0;

   // Now, go through keeping i as the most recent char in the text.
   while (i+1 < n)
   {
     // printf("\nStarting new Suffix: '%s'\n", text + i);
      l = 0;      
      
      // Find the first suffix which starts with the current symbol
      x = findStart(text[i], pattern, esa->SA, m);
      
      if (x<0)
      {
         P[p].j = -1;
         ++p;
         ++i;
         continue;
      }

      // printf("Starting with suffix: '%s'\n", pattern + SA[x]);
      P[p].i = i;
      
      // Extend the value as far as possible.
      while ( extend( text[++i], ++l, &x, pattern, esa, n, m ) );
            
      P[p].l = l;
      P[p].j = esa->SA[x];
   
      //  printf("Extended to length %d, using suffix:\n", P[p].l);
      //  printf("'%s'\n", pattern + SA[x]);
    
      ++p;
   }  
   
   if (p < n-1)
      P[p].j=-2;
   
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
   int FREQ_CHAR_THRESHOLD =  0.3*sqrt(m * log(m)/log(2) );
   
    printf("threshold: %d\n", FREQ_CHAR_THRESHOLD);
   
   // This will become a look up for each frequent character.
    
   int *pattern_lookup = malloc(sizeof(int)*FREQ_CHAR_THRESHOLD);

   // Go through every symbol, looking for frequent symbols.
   // then perform FFT matching on each one.
   for (int i=0; i<ALPHABET_SIZE; i++)
   {
   
      if (1) //frequency_table[i] >= FREQ_CHAR_THRESHOLD)
      {
         if (!(i >= 97 && i < 123)) continue;
         
         printf("method 1\n");      
         printf("%c\n", i);

         match_with_FFT(matches, i, text, pattern,  n, m);
      }
         else if (frequency_table[i] > 0)
      {
         printf("%c\n", i);
         printf("method 2\n");
         
         // Create lookup for this symbol.
         createLookup(pattern_lookup, i, pattern, m, FREQ_CHAR_THRESHOLD);
      
         // match this symbol in the text.
          markMatches(matches, text, i, pattern_lookup, n, FREQ_CHAR_THRESHOLD);
      }
   }
}                                


/******************************************************************************/
void kmismatches(         const char *text, 
                          const char *pattern,
                                int   k,
                                int   n,
                                int   m,
                                int  *matches  )
{
         
   // zero the matches array.
   memset(matches, 0, sizeof(int)*(n-m+1));

   // This will be a count of the occurences of symbols.
   int frequency_table[ALPHABET_SIZE];
   
   // Count the number of occurences of each syumbol.
   sp_km_count_symbols(pattern, m, frequency_table);

   // Number of appearances required for a character to be classed 'frequent'.
   int FREQ_CHAR_THRESHOLD = (int)(sqrt((double)k) + 0.5);
   
   // printf("threshold: %d\n", FREQ_CHAR_THRESHOLD);


   // Which k-mismatches case to perform.
   if ( count_frequent_symbols(frequency_table, FREQ_CHAR_THRESHOLD, n) > 2*sqrt(k) )
   {
   
      // This will become a look up for each frequent character.
      int *pattern_lookup = malloc(sizeof(int)*FREQ_CHAR_THRESHOLD);

      // Go through every symbol, looking for frequent symbols.
      // then perform FFT matching on each one.
      for (int i=0; i<ALPHABET_SIZE; i++)
      {
         if ( frequency_table[i] <= 0 ) continue;
      
         if ( frequency_table[i] > FREQ_CHAR_THRESHOLD )
         {
            printf("method 1\n");      
            printf("%c\n", i);

            match_with_FFT(matches, i, text, pattern,  n, m);
         }
            else
         {
            printf("%c\n", i);
            printf("method 2\n");
            
            // Create lookup for this symbol.
            createLookup(pattern_lookup, i, pattern, m, FREQ_CHAR_THRESHOLD);
         
            // match this symbol in the text.
            markMatches(matches, text, i, pattern_lookup, n, FREQ_CHAR_THRESHOLD);
         }
      }
   
   } else {
   
      printf("CASE 2, k = %d\n", k);
      k_mismatches_case2(text,pattern,frequency_table, k, n, m, matches);
   
   }

}                                


/******************************************************************************/

void display_pRepresentation(pTriple *P, const char *pattern, int n)
{
   printf("'");
   for (int i=0; i<n && P[i].j >=-1; i++)
   {
      // If the char wasn't in the text
      if (P[i].j == -1)
      {
         printf(" ");
         continue;
      } 
        
         
      for (int j=P[i].j; j< P[i].j + P[i].l; j++)
      {
         printf("%c", pattern[j]);
      }
   }
   printf("'");
}

/******************************************************************************/
void displaySA(int *SA, int *LCP, const char *pattern, int m)
{
   for (int i=0; i<m; i++)
   {
      printf("%d %s\n", LCP[i], pattern + SA[i]);
   } 
}

/******************************************************************************/

// j gives the location in the pattern, 
// the pair (x,t) give the location in the text (t is the index into the text
// where the x'th p-triple starts). 

int verifyMatch(  const pTriple  *pRepresentation,
                  const char *text,
                  const char *pattern,
                  const  ESA *esa,
                         int  x,
                         int  t,
                         int  i,
                         int  k,
                         int  n,
                         int  m                )
{
   // This implements the 'Kangarooing' method.
  
   // Call _i the position in the text, _j the position in the pattern.
   int _j = 0;
   
   // Get the index of the text substring.
   // note that i-t gives the 'lag' between the current text
   // position and the current p-block.
   int _i = pRepresentation[x].j + i-t;

   // The number of mismatches so far.
   int _k = 0;
   
   while (_j < m)
   {

      // jump to the next point where the text and pattern do not match.
      int temp = query(esa->SAi[_i]+1, esa->SAi[_j], esa->LCP, esa->n);
      int l    = esa->LCP[temp];
     
      // Does this run over the end of this p-triple? 
      // - if yes, we have found a mismatch, and need to 
      
      // Does this match carry on beyond the length of the pattern?
      // if yes then we are done.
      if ( l + _j >= m ) 
         return _k;

      // Does this run over the end of this p-triple?
      // If yes, then we know there is a mismatch, and 
      // we have to move to the next p-triple. (this happens at most three
      // times for any given verification).
      else if ( l > pRepresentation[x].l - (i-t))
      {
         // Increment number of mismatches found so far.
         _k ++;
         
         // Move to the next p-block.
         t += pRepresentation[x].l;
         x ++;
         
         // We are starting from the beginning of the next block, so there 
         // is no longer any lag between the text position and the start of
         // this block.
         i = t;
         
         // _i becomes the value j from the current p-block,
         // as we are starting from the beginning of it.
         _i = pRepresentation[x].j;
         
         // _j moves forwards one to go past this mismatching point.
         _j += l+1;
         
      }
      
      // We stopped because there was a mismatch. Remain in this p-block, 
      // increment k and continue.
      else 
      {  
         // Increment number of mismatches found so far.
         _k++;         
         
         // Move _j,_i along past this mismatch.
         _i += l+1;
         _j += l+1;  
      }
   }
   
   return _k;
}

/******************************************************************************/

// Construct an extended suffix array for some string of length n.s
void constructESA(const char *s, int n, ESA *esa)
{
   esa->n   = n;
   esa->SA  = malloc(sizeof(int) * (n+1));
   esa->SAi = malloc(sizeof(int) *  n   );
   esa->LCP = malloc(sizeof(int) *  n   );
     
   // Construct the SA and LCP in linear time.
   sais((unsigned char*)s, esa->SA, esa->LCP, n);
   
   // Construct SAi.
   for (int i=0; i<n; i++)
      esa->SAi[esa->SA[i]] = i;

}

/******************************************************************************/

// This is the second case in the k-mismatches algorithm.

void k_mismatches_case2(  const char *text, 
                          const char *pattern,
                          const int  *frequency_table,
                                int   k,
                                int   n,
                                int   m,
                                int  *matches          )
{

   // Find the positions where the pattern may match. 
   // We do this first for memory efficiency
   int sqrt_k = (int)(sqrt((double)k) + 0.5);
   
   printf("sqrt k: %d\n", sqrt_k);
   
   int *pattern_lookup = malloc(sizeof(int)*sqrt_k);

   // Find the first 2\sqrt{k} frequenct symbols, and mark all the positions
   // where they match.
   
   printf("Finding first %d characters and choosing first %d instances of them in the pattern\n", 2*sqrt_k, sqrt_k);
   for (int i=0, j=0; i<ALPHABET_SIZE &&  j< 2*sqrt_k; i++)
   {
      // Symbols that appear more than 
      if ( frequency_table[i] >= sqrt_k )
      {
         printf("%c is one of them\n", i);
         // Create lookup for this symbol.
         createLookup(pattern_lookup, i, pattern, m, sqrt_k);
      
         // match this symbol in the text.
         markMatches(matches, text, i, pattern_lookup, n, sqrt_k);
      }
   }

   //---------------//

   
   pTriple *pRepresentation = malloc(sizeof(pTriple) * n);

   // Construct the extended suffix array.
   ESA esa;   
   constructESA(pattern, m, &esa);
   
   construct_pRepresentation(pRepresentation, text, pattern, &esa, n, m);
     
   
   for (int i=0; i<n-m+1; i++)
   {
      printf("%d ", matches[i]);
   }  
   printf("\n");
     
   // INITIALISE THE RMQ structure so we can perform O(1) RMQ lookups.
   RMQ_succinct(esa.LCP, esa.n);  
     
   // We need to keep track of our location in the p-representation
   // AND the text.

   // The value t will keep track of the i position of the previously
   // seen p-triple.
   // The x value will keep track of the current position in the triple
   // array.
   int t=0;
   
   for (int i=0,x=0; i<n; i++)
   {
      if (x+1<n && (t + pRepresentation[x+1].l <= i))
      {
         ++x;
         t += pRepresentation[x].l;
      }
      
      // If there could be a possible match here.
      if (matches[i] >= k)
      {
      
         // Verify this location
         //for (int ]
      
      }
      
   
   
   }
     
     
   // printf("Actual: \n'%s'\n", text);  
   // display_pRepresentation(pRepresentation, pattern, n);
   
   // printf("\n\n");
   

}



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
   
   text[n-1] = 0;   
   
   for (i=0; i<m; i++)
      pattern[i] = (char)(rand() % 4 + 97);
   
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
   int n       = 20;

   //-------------------------------------------------------------------------//
 
   int x;

   // The text and pattern strings.
   char *t = malloc(sizeof(char) * (n+1));
   char *p = malloc(sizeof(char) * (m+1));


   randomStrings(t, p, n, m);
   printf("%s\n%s\n",t,p);
   
      int * matches  = malloc(sizeof(int) * (n-m+1));

   kmismatches(t,p,1,n,m,matches);

   exit(0);

   //printf("%s\n", t);
   printf("%s\n", p);

  n = loadData(&t, "./dna.50MB");
  printf("loaded, %d bytes\n", n);
   
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
      //   for (i=0; i<n-m+1; i++)
      //     printf("%d ", matches_FFT[i]);
      
      // printf("\n");
      
      //  for (i=0; i<n-m+1; i++)
      //    printf("%d ", matches_naive[i]);
      
      //  for (i=0; i<n-m+1; i++)
      //  assert(matches_naive[i] == matches_FFT[i]);   

      //printf("\n");
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
