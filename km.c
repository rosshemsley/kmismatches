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

/******************************************************************************
*
* Ideas for optimisations:
*
* - Improve p-representation construction: bin-search, child-tab
* - make non-relevant characters span multi blocks.
* - re-write the way that marking works - do two at once, or more?
* 
*
*******************************************************************************/

   double *FFT_match_p = NULL;
   double *FFT_match_t = NULL;
   double *FFT_match_r = NULL;

/******************************************************************************/
// Load a test input.

void load(                       const char     *filename, 
                                       int      *n, 
                                       int      *m, 
                                       int      *k, 
                                       int      *pos, 
                                       char    **text, 
                                       char    **pattern                       ) 
{

   FILE *f = fopen(filename, "r");

   if (f==NULL) {
      fprintf(stderr, "Failed to open file \"%s\"\n", filename);
      exit(1);
   }

   // The first line tells us how the rest of the file looks.  
   char buff[256];
   
   fgets(buff, 256, f);
   if (sscanf(buff, "%d %d %d %d", n, m, k, pos) != 4)
   {
      fprintf(stderr, "File header improperly formatted.\n");
      exit(1);
   }
  
   *text    = malloc(sizeof(char) * (*n+1));
   *pattern = malloc(sizeof(char) * (*m+1)); 
  
   fread(*pattern, sizeof(char), *m,   f);
   fgets( buff,    2,                  f);
   fread(*text,    sizeof(char), *n,   f);
   
   (*text)[*n]    = '\0';
   (*pattern)[*m] = '\0';
   
   if (strlen(*text) != *n || strlen(*pattern) != *m)
   {
      fprintf(stderr, "There was a problem reading the input file.\n");
      exit(1);
   }
   
   printf("===============================================================\n");
   printf("| Loaded test data.                                           |\n");
   printf("|                                                             |\n");
   printf("|       Text: %-10d bytes                                |\n", *n );
   printf("|    Pattern: %-10d bytes                                |\n", *m );
   printf("| Mismatches: %-10d                                      |\n", *k );
   printf("|   Position: %-10d                                      |\n", *pos);
   printf("|                                                             |\n");
   printf("===============================================================\n");
  
   // Take account of the end of termination of the string.
   (*text)[*n]    = '\0';
   (*pattern)[*m] = '\0';
   
   *n = *n+1;
   *m = *m+1;
  
}

/******************************************************************************/

void naive_matcher(              const char     *t, 
                                 const char     *p, 
                                       int       n, 
                                       int       m, 
                                       int      *A                             )
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

// Process has done i out of n rounds, and we want width w and resolution r.
static inline void loadBar(int x, int n, int r, int w)
{
   if ( x % (n/r) != 0 ) return;
   
   float ratio = x/(float)n;
   int   c     = ratio * w;

   printf("%3d%% [", (int)(ratio*100) );

 
   for (int x=0; x<c; x++)
      printf("=");
   
   for (int x=c; x<w; x++)
      printf(" ");   
  
   printf("]\n\033[F\033[J");
}

/******************************************************************************/

// Count the frequencies of symbols in t. 
// We assume that A is the same size as the alphabet. 
void sp_km_count_symbols(        const char     *t, 
                                       int       n, 
                                       int      *A                             )
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

void createLookup(                     int      *lookup, 
                                       char      symbol, 
                                 const char     *pattern, 
                                       int       m, 
                                       int       l                             )
{
   int x = 0;
   
   // zero the array.
   memset(lookup, 0, sizeof(int) * l);
   
   for (int i=0; i<m; i++)
   {
      if (pattern[i] == symbol)
      {
         lookup[x] = i;
         ++x;
      }
      
      if (x >=l ) break;
   }

   if (x<l)
   lookup[x]=-1;
   
}

/******************************************************************************/
// attempts to optimise this function

// This seems slower.
void markMatches2(       int  *matches, 
                  const char *text, 
                        char  symbol, 
                  const int  *lookup, 
                        int   n,
                        int   m, 
                        int   l       )
{

   int end = l;
   
   // Find the last char in the lookup. (Optimisation).
   for (int i=0; i<l;i++)
   {
      if (lookup[i] == -1) 
      {
         end = i;
         break;
      }  
   }


   // NOTICE OPTIMISED BOUNDARIES: we infer the max and min 
   // positoins in the text.
   for (int i=lookup[0]; i <  (n-lookup[l-1]+1);    i++)
   {  
      // Perform the marking.
      if (text[i] == symbol)
      {
         for (int j=0; j<end; j++)
         {
            // TODO: CHECK THIS IS RIGHT
            if ( i - lookup[j] >= n-m+1) break;
            
            if ( i - lookup[j] >= 0 )
               ++ matches[i-lookup[j]];
         }
      }      
   }
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

   int end = l;
   
   // Find the last char in the lookup. (Optimisation).
   for (int i=0; i<l;i++)
   {
      if (lookup[i] == -1) 
      {
         end = i;
         break;
      }  
   }

   for (int i=0; i<n; i++)
   {  
      // Perform the marking.
      if (text[i] == symbol)
      {
         for (int j=0; j<end; j++)
         {
            // TODO: CHECK THIS IS RIGHT
            if (i-lookup[j] >= n-m+1) break;
            
            if ( i - lookup[j] >= 0 )
               ++ matches[i-lookup[j]];
         }
      }      
   }
} 

/******************************************************************************/

int count_frequent_symbols(       
                            const int *frequency_table, 
                                  int  threshold)
{
   int count = 0;
   
   for (int i=0; i<ALPHABET_SIZE; i++)
      if (frequency_table[i] > threshold) count ++;
   return count;
}

/******************************************************************************/

// This will take the current position in the SA, x and the current
// pTriple p and the next letter in the text. It will then loop through the 
// following sufficies until the first suffix which has the correct next letter
// If this is found, we return the suffix and extend l by one. If the LCP becomes
// to small (i.e. the preceeding chars don't match) then we return -1, 
// and must start a new pTriple.

static inline int extend(              char      t, 
                                       int       l,
                                       int      *x,
                                 const char     *pattern, 
                                 const ESA      *esa,
                                       int       n,     
                                       int        m                            )
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
      //printf("%d, %d, l: %d. x: %d, SA[x]: %d",
      // esa->SA[*x]+l, m, l, *x, esa->SA[*x]);
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

static inline int findStart(           char      c, 
                                 const char     *pattern, 
                                 const int      *SA, 
                                       int       m                             )
{

   // Bin-search:
   
   int min = 0;
   int max = m-1;
   int mid; 
   
   do
   {
      mid = (max+min)/2;
   
      if (pattern[SA[mid]] == c)
         return mid;   

      else if (c > pattern[SA[mid]])
         min = mid+1;
         
      else
         max = mid-1;
         
   } 
   while (min <= max);
   

/*
   // For now we are lazy and do a linear search.
   for (int i=0; i<m; i++)
      if (pattern[SA[i]] == c) {
         if(test!=i) printf ("ERROR: found: %d, actual: %d\n", test, i);
         printf("min: %d, max: %d\n", min, max);
         exit(0);
      }
  
*/
  
   return -1;
}

/******************************************************************************/


void construct_pRepresentation(        pTriple  *P,
                                 const char     *text, 
                                 const char     *pattern, 
                                 const ESA      *esa,
                                       int       n,
                                       int       m                             )
{

   int p = 0;
   int i = 0;
   int x = 0;
   int l = 0;


   /* OPTIMSATION *************************************************************/
   
   // Create a look-up table for the first positions of each character in 
   // the SA.
   int LOOKUP[ALPHABET_SIZE];
   
   // Set every item to -1.
   for (int i=0; i<ALPHABET_SIZE; i++)
      LOOKUP[i] = -1;
   
   // Now, set those characters that appear to the correct value.
   for (int i=0; i<m; i++)
   {
      // only when the LCP is 0 do we have a new character.
      if (esa->LCP[i]==0)
         LOOKUP[(unsigned char)pattern[esa->SA[i]]] = i;
   
   }
   
   /***************************************************************************/

   // Now, go through keeping i as the most recent char in the text.
   while (i+1 < n)
   {
   
    //  loadBar(i,n,100,40);
      
     // printf("\nStarting new Suffix: '%s'\n", text + i);
      l = 0;      
      
      // Find the first suffix which starts with the current symbol
      // x = findStart(text[i], pattern, esa->SA, m);
      
      x = LOOKUP[(unsigned char)text[i]];
      
      if (x<0)
      {
         P[p].j = -1;
         P[p].l =  1;
         ++p;
         ++i;
         continue;
      }

      // printf("Starting with suffix: '%s'\n", pattern + esa->SA[x]);
        
      // Extend the value as far as possible.
      while ( (i<n-1) && extend( text[++i], ++l, &x, pattern, esa, n, m ) );
            
      P[p].l = l;
      P[p].j = esa->SA[x];
      ++p;
   }  
   
   // Mark the end of the p-representation.
   if (p < n-1)
      P[p].j=-2;
      
    printf("Length of p-representation: %d for %d symbols\n", p, n);

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
      printf("%d %s\n", LCP[i], pattern + SA[i]);
}

/******************************************************************************/
// Find the Longest Common Extension using the Extended Suffix Array.

static inline int LCE(int i, int j, const ESA *esa)
{
   // Trivial query.
   if (i == j)
      return (esa->n - j);     

   // Make sure the indicies are the right way around.
   register int a = esa->SAi[i];
   register int b = esa->SAi[j];
     
   if (a>b)
   {
      register int c = a;
      a = b;
      b = c;
   }
   
   register int temp =  query(a+1, b, esa->LCP, esa->n); 
                                   //*/ query_naive( a+1, b, esa->LCP, esa->n );
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
   while (mismatches<=k && j < m-1)
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
         ++i;
         ++t;
         ++j;
         
         // We consider this a mismatch.         
         ++mismatches;
         continue;
      }
          
      // Don't use RMQ to check if one character matches. (Optimsation).
      if (block_start == block_end)
      {
         if (pattern[j] != pattern[block_start])
            ++mismatches;
         
         ++t;
         ++i;
         ++j;
         ++x;
         
         block_start = pRepresentation[x].j;
         block_end   = pRepresentation[x].j + pRepresentation[x].l -1;   
         continue;    
      }
      
      // Find the longest common extension between current positions.
      register int l = LCE(block_start, j, esa);

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
     
         // REMOVE POINTLESS QUERIES FROM END. (Optimisation).
         if (block_start + l == block_end)
         {
            i += l+1;
            t += l+1;
            j += l+1;
            
            // Move to the next block.
            ++x;
            block_start = pRepresentation[x].j;
            block_end   = pRepresentation[x].j + pRepresentation[x].l -1;       
            
            mismatches++;
            continue;  
         }
       
         // Advance the position by l+1.
         i           += l+1;
         j           += l+1;         
         block_start += l+1;

         // Register the mismatch we found.
         ++mismatches;   
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
   int sqrt_k = (int)(sqrt((double)k) + 0.5);
   
    printf("Threshold: \\sqrt k %d.\n", sqrt_k);

   int num_freq_chars = count_frequent_symbols(frequency_table, sqrt_k);
   
   printf("There are %d frequent characters\n", num_freq_chars);

   // Which k-mismatches case to perform.
   if (num_freq_chars < 2* sqrt_k)
   {
      printf("CASE 1\n");
   
      // This will become a look up for each frequent character.
      int *pattern_lookup = malloc(sizeof(int)*sqrt_k);

      // Go through every symbol, looking for frequent symbols.
      // then perform FFT matching on each one.
      for (int i=0; i<ALPHABET_SIZE; i++)
      {
         if ( frequency_table[i] <= 0 ) continue;
      
         if ( frequency_table[i] > sqrt_k )
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
            createLookup(pattern_lookup, i, pattern, m, sqrt_k);
         
            // match this symbol in the text.
            markMatches(matches, text, i, pattern_lookup, n,m, sqrt_k);
         }
      }
   
   } else {
   
      printf("CASE 2\n");
      k_mismatches_case2(text, pattern, frequency_table, k, n, m, matches);
   
   }
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
   
   // Create a look up matrix with 2k positions.
   int *pattern_lookup = malloc(sizeof(int)*sqrt_k*2*sqrt_k);

   // Find the first 2\sqrt{k} frequenct symbols, and mark all the positions
   // where they match.
   
   printf("Finding first %d characters and choosing first %d "                
           "instances of them in the pattern\n", 2*sqrt_k, sqrt_k);
           
           
    // This will map symbols in the alphabet to look up tables. 
    // -1 means ignore.
    int LOOKUP[ALPHABET_SIZE];
    for (int i=0; i<ALPHABET_SIZE; i++)
      LOOKUP[i] = -1;
           
   for (int i=0, j=0; i<ALPHABET_SIZE &&  j<2*sqrt_k; i++)
   {
      // Symbols that appear more than 
      if ( frequency_table[i] >= sqrt_k )
      {
         printf("%c is a frequent character.\n", i);
         
         // Create lookup for this symbol.
         // we store the look-ups in columns. 
         createLookup(pattern_lookup + j*sqrt_k, i, pattern, m, sqrt_k);

         // The position of the look up table for this character.
         LOOKUP[i] = j*sqrt_k;
         
         j++;
      }
   }
   
   printf("Doing marking\n");
   for (int i=0; i<n; i++)
   {
      // If this symbol is one of our look-up characters.
      int l = LOOKUP[(unsigned char)text[i]];
     
      if (l != -1)
      {
         int *this_table = pattern_lookup + l;
         
         for (int j=0; j<sqrt_k; j++)
         {
            // TODO: CHECK THIS IS RIGHT
            if (i-this_table[j] >= n-m+1) break;            
            
            if ( i - this_table[j] >= 0 )
               ++ matches[i-this_table[j]];
         }
      }
   }
   
   // We no longer need this.
   free(pattern_lookup);   

   printf("Done Marking\n");
   //---------------//
   
   pTriple *pRepresentation = malloc(sizeof(pTriple) * n);

   // Construct the extended suffix array.
   
   printf("Constructing ESA and p-representation\n");
   ESA esa;   
   constructESA(pattern, m, &esa);
   construct_pRepresentation(pRepresentation, text, pattern, &esa, n, m);
   printf("Done\n");
   
   
   //printf("%s\n", text);
   //display_pRepresentation(pRepresentation, pattern, n);

  // exit(0);

     
   // INITIALISE THE RMQ structure so we can perform O(1) RMQ lookups.
   RMQ_succinct(esa.LCP, esa.n);  
     
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
      
      //printf("%d\n", matches[i]);
      
      // If there could be a possible match here.
      if (matches[i] >= k)
      {
      
         printf("\nVerifying position %d\n", i);

         // Verify this location    
         matches[i] = verifyMatch(pRepresentation, text, pattern, &esa, x, t, i, k, n, m);     
      } else matches[i] = k+1;
         
   }
   
   freeESA(&esa);

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
            if (mismatches > k) break;  
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
   printf("Done ESA\n");
  
   RMQ_succinct(esa.LCP, esa.n); 
   printf("Done RMQ\n");

   pTriple *pRepresentation = malloc(sizeof(pTriple) * n);  
  
   // Construct the p-representation.
   printf("Constructing p-Representation.\n");
   construct_pRepresentation(pRepresentation, text, pattern, &esa, n, m);
   
   int t=0;
   
   // Now, go through every position and look for mismatches.
   printf("Looking for k-Mismatches.\n");
   for (int i=0,x=0; i<n-m+1; i++)
   {
      loadBar(i,n,100,40);
   
      // Advance through the p-reprsentation.
      if (t + pRepresentation[x].l <= i)
      {
         t += pRepresentation[x].l;
         ++x;
      }
           
      int v = verifyMatch(pRepresentation, text, pattern, &esa, x,t,i,k,n,m);

      matches[i] = v;
      
   }
   
   free(pRepresentation);
   freeESA(&esa);
}

/******************************************************************************/
//
// NOTE: SAIS APPEARS TO FAIL WITH INPUT 'aacdbbcca', 'dddbbcddc'
//
//
/******************************************************************************/

int main(int argc, char **argv)
{

   if (argc !=2)
   {
      fprintf(stderr, "No input file provided\n");
      exit(1);
   }
   
   
   // the length of the text and pattern.
   int m;
   int n;
   int k;
   int pos;

   // The text and pattern strings.
   char *t = NULL;
   char *p = NULL;

   // Load the test file.
   load(argv[1], &n, &m, &k, &pos, &t, &p);

   // An array to put the mismatches in.
   int  *matches        = malloc(sizeof(int)  * (n-m+1));

   // Perform Kangarooing.
   kmismatches(t,p,k,n,m,matches);

      
   // Verify the output.   
   printf("\nCHECKING: \n");

   int match_pos = -1;
   int match_k   = -1;
   
   for (int b=0;b<n-m+1; b++)
   {
      if (matches[b] <=k) 
      {
         match_pos = b;
         match_k   = matches[b];
         break;
      }  
   }
   
   printf("Position: %d, mismatches: %d\n", match_pos, match_k);  
   if (match_pos != pos || match_k != k)
   {
      printf("FAILED TEST\n");
      exit(1);
   }  
   
   
   free(p);
   free(t);
   free(matches);
   
   // This needs to be fixed.
   // FreeRMQ_succinct();
   
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
