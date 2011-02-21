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

#define DEBUG

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
   
   for (int i=0; i<l; i++)
   {
      if (pattern[i] == symbol)
      {
         lookup[x] = i;
         ++x;
      }
   }

   if (x<l)
   lookup[x]=-1;
   
}

/******************************************************************************/

void markMatches(       int  *matches, 
                  const char *text, 
                        char  symbol, 
                  const int  *lookup, 
                        int   n,
                        int   m, 
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
            if ( i - lookup[j] >= 0 && i-lookup[j] < n-m+1 )
            {

               matches[i-lookup[j]] ++;
            }
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

static inline int extend     (           char   t, 
                           int    l,
                           int   *x,
                     const char  *pattern, 
                     const  ESA  *esa,
                           int    n,     
                           int    m         )
{

   //printf("Extending\n");  
   //printf("l: %d\n", l); 
   // Go through the Suffix Array until LCP[i] < l or pattern[SA[i]+l] = t.
   
   // Check to see whether or not we can extend the current suffix.   
   if ( (esa->SA[*x]+l < m) && (pattern[esa->SA[*x] + l] == t) )
   {
       //printf("Succeeded with current\n");
       //printf("'%s'\n", pattern + esa->SA[*x]); 
      return 1;
   }
   
   // Look for other possible sufficies which may match.
   for (int i = *x+1; i<n; i++)
   {
      // If we have run out of sufficies.
      if (i>=m) return 0;
      
      //printf("Checking against: %d\n", i);
      // There are no sufficies with this prefix.
      if (esa->LCP[i] < l ) {
      
         //printf("Ran out of values\n");
         return 0;
      }
      
      // Check this character match.
      //printf("Comparing '%c' against '%c'\n", pattern[esa->SA[i]+l], t);
      //printf("%d, %d, l: %d. x: %d, SA[x]: %d", esa->SA[*x]+l, m, l, *x, esa->SA[*x]);
      if ( (esa->SA[i]+l < m) && (pattern[esa->SA[i] + l] == t) )
      {
         //printf("Changing to '%s'\n", pattern + esa->SA[i]); 
         *x = i;
         return 1;
      }
   }

   return -1;
}

/******************************************************************************/

static inline int findStart(char c, const char *pattern, const int *SA, int m)
{

   if (c!='A' && c!= 'C' && c!='T' && c!='G') return -1;
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

   int p = 0;
   int i = 0;
   int x = 0;
   int l = 0;

   // Now, go through keeping i as the most recent char in the text.
   while (i+1 < n)
   {
      if (i%1000000 ==0)
      {
         printf("Done %d%%\n", (int)(i/(float)n*100));
      }
   
     // printf("\nStarting new Suffix: '%s'\n", text + i);
      l = 0;      
      
      // Find the first suffix which starts with the current symbol
      x = findStart(text[i], pattern, esa->SA, m);
      
      if (x<0)
      {
         P[p].j = -1;
         P[p].l =  1;
         ++p;
         ++i;
         continue;
      }

      // printf("Starting with suffix: '%s'\n", pattern + esa->SA[x]);
      P[p].i = i;
        
      // Extend the value as far as possible.
      while ( (i<n-1) && extend( text[++i], ++l, &x, pattern, esa, n, m ) );
            
      P[p].l = l;
      P[p].j = esa->SA[x];
      ++p;
   }  
   
   // Mark the end of the p-representation.
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
         markMatches(matches, text, i, pattern_lookup, n, m, FREQ_CHAR_THRESHOLD);
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
   if (0) // count_frequent_symbols(frequency_table, FREQ_CHAR_THRESHOLD, n) > 2*sqrt(k) )
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
            markMatches(matches, text, i, pattern_lookup, n,m, FREQ_CHAR_THRESHOLD);
         }
      }
   
   } else {
   
   //   printf("CASE 2, k = %d\n", k);
      k_mismatches_case2(text,pattern,frequency_table, k, n, m, matches);
   
   }

}                                


/******************************************************************************/

void display_pRepresentation(const pTriple *P, const char *pattern, int n)
{

   for (int i=0; i<n && P[i].j >=-1; i++)
   {
      // If the char wasn't in the text
      if (P[i].j == -1)
      {
         printf("_");
         continue;
      } 
        
      for (int j=P[i].j; j< P[i].j + P[i].l; j++)
      {
         printf("%c", pattern[j]);
      }
      printf("|");
   }
   printf("\n");
}

/******************************************************************************/

int count_naive(const char *t, const char *p, int k, int l)
{
   int count =0;
   
   for (int i=0; i<l; i++)
   {
      if (count >=k) return count;
      if (t[i] != p[i]) count ++;
   }
   
   return count;
   
 

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
// Find the Longest Common Extension using the Extended Suffix Array.

static inline int LCE(int i, int j, const ESA *esa)
{
   // Trivial query.
   if (i == j)
      return (esa->n - j);     

   // Make sure the indicies are the right way around.
   int a = esa->SAi[i];
   int b = esa->SAi[j];
     
   if (a>b)
   {
      int c = a;
      a = b;
      b = c;
   }
   
   int temp =  /*query(a+1, b, esa->LCP, esa->n); //*/ query_naive( a+1, b, esa->LCP, esa->n );
   return     esa->LCP[temp];
               
}

/******************************************************************************/

// j gives the location in the pattern, 
// the pair (x,t) give the location in the text (t is the index into the text
// where the x'th p-triple starts). 

static inline int verifyMatch(  const pTriple  *pRepresentation,
                  const char     *text,
                  const char     *pattern,
                  const ESA      *esa,
                  
                        // The position in the p-representation.              
                        int       x,    /*  The current p-block               */
                        int       t,    /*  The i position of this block      */
                        int       i,    /*  The current position in the text  */
                        
                        // problem-specific variables.
                        int       k,
                        int       n,
                        int       m                                            )
{
   // The current position in the pattern.
   int j = 0;
   
   // The count of the mismatches found.
   int mismatches = 0;
   
   // The start and end of the p-block for this part of the text.
   int block_start = pRepresentation[x].j + (i-t);
   int block_end   = pRepresentation[x].j + pRepresentation[x].l-1;
      
   // Look through all the characters in the pattern.
   // NOTE: We assume the last char is \0.
   while (j < m-1)
   {
      // Ignore filtered characters:
      // TODO: allow blocks of length > 1 for ignored characters.
      if (block_start == -1)
      {
         // Move to the next block
         ++x;
         block_start = pRepresentation[x].j;
         block_end   = pRepresentation[x].j + pRepresentation[x].l-1;

         // Advance the position by 1.
         i += 1;
         t += 1;
         j += 1;
         
         // We consider this a mismatch.         
         mismatches ++;
         if (mismatches >= k) return k;      
         continue;
      }
      
      // Find the longest common extension between current positions.
      int l = LCE(block_start, j, esa);

      // If the length goes over the length of this block, we 
      // just give the maximum possible length.
      if (l + block_start > block_end)
         l = block_end - block_start+1;
          
      // If this takes us over the end of the pattern, return.
      // Otherwise we might end up incrementing mismatches too many times.
      // TODO: Do this in a more efficient way.
      if (j + l  == m-1) break;
      
      
      /**
      *   CASE 1: Substrings match all the way to the end of their current block
      *   We start the next block and continue.
      */
      if (block_start + l > block_end)
      {
         // Advance position by l.
         i += l;
         t += l;
         j += l;
         
         // Move to the next block.
         ++x;
         block_start = pRepresentation[x].j;
         block_end   = pRepresentation[x].j + pRepresentation[x].l -1;        
      } 
      /**
      *  CASE 2: It mismatches within the current block,
      *  We increment k and continue in this block.
      */
         else
      {
         // Advance the position by l+1.
         i           += l+1;
         j           += l+1;         
         block_start += l+1;

         // Register the mismatch we found.
         mismatches ++;   
         if (mismatches >= k) return k;      
      }   
   }

   return mismatches;
}

/******************************************************************************/

void freeESA(ESA *esa)
{
   free(esa->SA);
   free(esa->LCP);
   free(esa->SAi);
}

/******************************************************************************/

// Construct an extended suffix array for some string of length n.s
void constructESA(const char *s, int n, ESA *esa)
{
   
   esa->n   = n;
   
   // TODO: Change these to malloc's later.
   esa->SA  = calloc( (n+1), sizeof(int) );
   esa->SAi = calloc( (n+1), sizeof(int) );
   esa->LCP = calloc( (n+1), sizeof(int) );
     
   // Construct the SA and LCP in linear time.
   sais((unsigned char*)s, esa->SA, esa->LCP, n);

   // Construct SAi.
   for (int i=0; i<n; i++)
      esa->SAi[esa->SA[i]] = i;

   
   /*
   for (int i=0; i<n; i++)
   {
      printf("%d: %d   (%d)%s\n", i, esa->LCP[i], esa->SA[i], s + esa->SA[i]);
   }
   printf("\n");
   for (int i=0; i<n; i++)
      printf("%d: %d\n", i, esa->SAi[i]);
   */
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
   
 //  printf("sqrt k: %d\n", sqrt_k);
   
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
         markMatches(matches, text, i, pattern_lookup, n,m, sqrt_k);
      }
   }

   //---------------//

   
   pTriple *pRepresentation = malloc(sizeof(pTriple) * n);



   // Construct the extended suffix array.
   ESA esa;   
   constructESA(pattern, m, &esa);
   
   construct_pRepresentation(pRepresentation, text, pattern, &esa, n, m);
   
   //printf("%s\n", text);
   //display_pRepresentation(pRepresentation, pattern, n);

  // exit(0);
   for (int i=0; i<n-m+1; i++)
   {
      printf("%d ", matches[i]);
   }  
   printf("\n");
     
   // INITIALISE THE RMQ structure so we can perform O(1) RMQ lookups.
//   RMQ_succinct(esa.LCP, esa.n);  
     
   // We need to keep track of our location in the p-representation
   // AND the text.

   // The value t will keep track of the i position of the previously
   // seen p-triple.
   // The x value will keep track of the current position in the triple
   // array.
   int t=0;
   
   for (int i=0,x=0; i<n-m+1; i++)
   {
      if (t + pRepresentation[x].l <= i)
      {
         
         t += pRepresentation[x].l;
         ++x;
      }
      
      // If there could be a possible match here.
      if (matches[i] >= k)
      {
         printf("\nVerifying position %d\n", i);
            display_pRepresentation(pRepresentation, pattern, n);
         // Verify this location    

         
         int actual = count_naive(text+i, pattern, k, m-1);

         if (actual != verifyMatch(pRepresentation, text, pattern, &esa, x, t, i, k, n, m))
         {
           printf("ERROR\n");
           exit(0);
          
         }
     

     
     
      }
      
   
   
   }
     
     
   // printf("Actual: \n'%s'\n", text);  
   // display_pRepresentation(pRepresentation, pattern, n);
   
   // printf("\n\n");
   
   
   freeESA(&esa);
   free(pattern_lookup);
   free(pRepresentation);
}

/******************************************************************************/

void naive_kangaroo (     const char *text, 
                          const char *pattern,
                                int   k,
                                int   n,
                                int   m,
                                int  *matches          )
{

      
      
   memset(matches, 0, sizeof(int)*(n-m+1));
   for (int i=0; i<n-m+1; i++)
   {
      int mismatches =0;
      for (int j=0; j<m-1; j++)
      {
         if (pattern[j] != text[i+j]) 
         { 
            mismatches ++;
            if (mismatches >= k) break;  
         }
      }
   
      matches[i] = mismatches;
   
   }

}                                
                               

/******************************************************************************/
// Use kangarooing to achieve O(nk) k-mismatches.

void kangaroo(            const char *text, 
                          const char *pattern,
                                int   k,
                                int   n,
                                int   m,
                                int  *matches          )

{   
   // Zero the matches array.
   memset(matches, 0, sizeof(int)*(n-m+1));
   
   
   // Construct the extended suffix array.
   ESA esa;   
   constructESA(pattern, m, &esa);
  
  
  
 //  RMQ_succinct(esa.LCP, esa.n); 
  
   printf("Done RMQ\n");
  
   printf("Done ESA\n");
  
   pTriple *pRepresentation = malloc(sizeof(pTriple) * n);  
  
   // Construct the p-representation.
   construct_pRepresentation(pRepresentation, text, pattern, &esa, n, m);
   
   
   printf("Done P-Rep\n");
   
 //    display_pRepresentation(pRepresentation, pattern, n);
   
   int t=0;
   // Now, go through every position and look for mismatches.
   for (int i=0,x=0; i<n-m+1; i++)
   {
      if (i%10 == 0)
      printf("i: %d\n", i);
      // Advance through the p-reprsentation.
      if (t + pRepresentation[x].l <= i)
      {
         t += pRepresentation[x].l;
         ++x;
      }
           
      int v = verifyMatch(pRepresentation, text, pattern, &esa, x,t,i,k,n,m);

      // DEBUGGING.
      int actual = count_naive(text+i, pattern, k, m-1);
      if (actual != v)
      {
         printf("ERROR: %d, %d\n", actual, v);
         exit(0);    
      }
      
      matches[i] = v;
      
   }
   
   free(pRepresentation);
   freeESA(&esa);
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
   
   if (r==0) return 'G';
   if (r==1) return 'T';
   if (r==2) return 'A';
   if (r==3) return 'C';

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
      text[i]    = randDNA(); //(char)(rand() % 4 + 97);
   
   text[n-1] = '\0';   
   
   for (i=0; i<m; i++)
      pattern[i] = randDNA(); //(char)(rand() % 4 + 97);
   
   pattern[m-1] = '\0';

}


/******************************************************************************/
//
// NOTE: SAIS APPEARS TO FAIL WITH INPUT 'aacdbbcca', 'dddbbcddc'
//
//
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
   int m       = 100;
   int n       = 100000;

   //-------------------------------------------------------------------------//
 
   // The text and pattern strings.
   char *t       = malloc(sizeof(char) * n);
   char *p        = malloc(sizeof(char) * m);
   int  *matches  = malloc(sizeof(int)  * (n-m+1));

   for (int a=0; a < repeats; a++)
   {
    randomStrings(t, p, n, m);   
   
      //printf("Pattern: %s\n", p);
   
   //  n = loadData(&t, "./dna.50MB");
    // printf("loaded, %d bytes\n", n);
   

   
      int  *matches        = malloc(sizeof(int)  * (n-m+1));
      int  *matches_naive  = malloc(sizeof(int)  * (n-m+1));
   
      printf("Doing naive.\n");
      naive_kangaroo(t,p,100,n,m,matches_naive);
            kangaroo(t,p,100,n,m,matches);
      printf("Doing kangaroo.\n");
      
      
      for (int b=0;b<n-m+1; b++)
      {
         if (matches_naive[b] != matches[b])
         {
            fprintf(stderr, "%d: MATCHES NOT EQUAL: %d, %d\n",b, matches_naive[b], matches[b]);
            
            exit(0);
         }
      }
   }
   
   free(p);
   free(t);
   free(matches);
  // FreeRMQ_succinct();
   
   /*
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
        
   */
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
