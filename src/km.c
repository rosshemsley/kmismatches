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
#include "esa.h"
#include "stack.h"
#include "breaks.h"
#include "loadTest.h"
#include "RMQ_succinct.h"
#include "sp_km_unbounded_matcher.h"

#define DEBUG

/*******************************************************************************
*
* Ideas for optimisations:
*
* X Improve p-representation construction: bin-search, child-tab
* - make non-relevant characters span multi blocks.
* X re-write the way that marking works - do two at once, or more?
* - change the way p-blocks are written, since the length of a block is
* - rarely longer than 256.
* - Same with LCP
*
*******************************************************************************/





/******************************************************************************/
// Display the status of a long-running process.

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
 
void sp_km_count_symbols(        const char*     t, 
                                       int       n, 
                                       int*      A                             )
{
   int i=0;
   
   for (i=0; i<ALPHABET_SIZE; i++)
      A[i] = 0;
   
   // NOTE, WE DO NOT COUNT LAST SYMBOL.
   for (i=0; i<n-1; i++)
      A[(unsigned char)t[i]] ++;
}

/******************************************************************************/
// Store the indicies where the symbol occurs in the pattern at each location
// in the lookup table.

void createLookup(                     int*      lookup, 
                                       char      symbol, 
                                 const char*     pattern, 
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
// First attempt at optimising this function
// Currently seems to be slower: TODO: improve the way the optimisations
// are implemented!
// - These ideas may be significantly faster for large pattern sizes.

void markMatches2(                     int*      matches, 
                                 const char*     text, 
                                 char            symbol, 
                                 const int*      lookup, 
                                       int       n,
                                       int       m, 
                                       int       l                             )
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
   // positions in the text.
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
// Count how many frequent symbols their are by going through the lookup table.

int count_frequent_symbols(      const int*      frequency_table, 
                                       int       threshold                     )
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
// If this is found, we return the suffix and extend l by one. If the LCP 
// becomes to small (i.e. the preceeding chars don't match) then we return -1, 
// and must start a new pTriple.

static inline int extend(              char      t, 
                                       int       l,
                                       int*      x,
                                 const char*     pattern, 
                                 const ESA*      esa,
                                       int       n,     
                                       int       m                             )
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

// Find the maximum possible extension of any suffix and this part of the text 
// and store it in P.
static inline int extendInterval(const int*      LOOKUP, 
                                       pTriple*  P, 
                                 const char*     text,
                                 const char*     pattern, 
                                 const int       n,  
                                 const int       m, 
                                 const ESA*      esa                           )
{
 
   // We start with the first l-interval.
   int i = LOOKUP[(unsigned char)text[0]];
   int j = esa->accross[i]-1;
   
   // This is the current length of the matching prefix of the text.
   int l = 1;
   
   // u is the start of the current l-index.
   int u = i;   
   
   // v is the start of the next l-index.
   int v = LI_FIRST_CHILD(i, j, esa);
      
   // We use this to force a stop on the next iteration.
   int STOP = 0;

   // Now, keep traversing the tree until we run out of possibilities. 
   while (1)
   {
      
      // Singleton.
      // We have reached a leaf node.
      if (i == j)
      {  
        // printf("Singleton interval\n");
         while  (text[l-1] != '\0' && (pattern + esa->SA[u])[l] == text[l])
            l++;
         break;
      }
                 
      // Find the l-value of this l-interval.
      int lcp = LI_GET_LCP(i,j,esa);
      
      // Set this to zero if the matching failed.
      int match = 1;
      
      // This will walk us down branches which are multiple symbols long.
      for (; l < lcp ; l++)
      {
         if ((pattern + esa->SA[u])[l] != text[l])
         {
            match = 0;
            break;
         }
      }
      
      if (match == 0)
         break;
         
      // Does the start of this l-interval match?
      if ( (pattern + esa->SA[u])[l] == text[l] )       
      {
         // Increment the 'depth' in the text.
         ++l;      
         
         // If we are in the last interval and find a match, then 
         // we want to continue looking down the tree and so we do not stop yet.
         STOP = 0;

         // Move down the tree.
         i = u;
         j = v-1;
         
         // find the index of the first child interval in the new l-interval.
         v = LI_FIRST_CHILD(i,j,esa);
        
         continue;
      }  
      
      // Move accross the tree to the next child interval.
      u = v;
      v = esa->accross[v];
   
      // If we have reached the end of this l-interval,
      // then v will be 0. However, there could still be a match
      // here, so we use STOP to indicate that the next round should
      // be the last.
      // Obviously, if we find a match, STOP should be reset.            
      if (STOP)
         break;
      
      // We have reached the end of this l-interval.
      // Setting v=j+1 will mean that the end of the interval will be set to 
      // j later on.
      // We want to stop on the next iteration of there are no matches.
      if (v == 0)
      {
         v = j+1;
         STOP = 1;
      }       
   }
   
   // put our 'results' into this p-block.
   P->l = l;
   P->j = esa->SA[i];

   // Return the longest match we found.
   return l;

}

/******************************************************************************/

void construct_pRepresentation(        pTriple*  P,
                                 const char*     text, 
                                 const char*     pattern, 
                                 const ESA*      esa,
                                       int       n,
                                       int       m                             )
{

   // Look up table for each character, giving first occurence in the SA.
   int LOOKUP[ALPHABET_SIZE];
   
   // Set every item to -1.
   for (int i=0; i<ALPHABET_SIZE; i++)
      LOOKUP[i] = -1;
   
   // Now, set those characters that appear to the correct value.
   for (int i=0; i<m; i++)
   {
      // only when the LCP is 0 do we have a new character.
      if (esa->LCP[i] == 0)
         LOOKUP[ (unsigned char)pattern[esa->SA[i]] ] = i;   
   }

   // Position in the text.
   int t = 0;
   // Position in the p-Representation.
   int x = 0;
   
   // Go through every value in the text.
   while (t < n-1)
   {      
      // Is the symbol in the pattern?
      if (LOOKUP[(unsigned char)text[t]] == -1)
      {
         P[x].j = -1;
         P[x].l = 1;         
         ++t;
      } 
         else
      {       
         t += extendInterval(LOOKUP, &P[x], text + t, pattern, n, m, esa);
      }
        
      ++x;           
   }
      
   // Terminate the p-representation if it is shorter than the worst case.
   if (x<n);
      P[x-1].l = 0;
      P[x-1].j = -2;

}

/******************************************************************************/
// Old, naive p-representation construction.

void construct_pRepresentation_old(    pTriple*  P,
                                 const char*     text, 
                                 const char*     pattern, 
                                 const ESA*      esa,
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

void abrahamson_kosaraju(        const char*     text, 
                                 const char*     pattern,
                                       int       n,
                                       int       m,
                                       int *     matches                       )
{

   
         
   // Zero the matches array.
   memset(matches, 0, sizeof(int)*(n-m+1));

   // This will be a count of the occurences of symbols.
   int lookup[ALPHABET_SIZE];
   
   // Count the number of occurences of each symbol and store in lookup
   sp_km_count_symbols(pattern, m, lookup);

   int threshold =  sqrt( m * log(m)/log(2) );
         
   printf("Threshold: %d\n", threshold);



   /*--- Frequent characters ---  */

   // Go through every symbol, looking for frequent symbols.
   // then perform FFT matching on each one.
   for (int i=0; i<ALPHABET_SIZE; i++)
   {
      if (lookup[i] >= threshold)
      {  
         printf("Method 1\n");      
         printf("%c\n", i);

         match_with_FFT(matches, i, text, pattern,  n, m);
      }
   }   
 
   /*--- Infrequent characters ---*/
   
   // Now, deal with the infrequent characters using a simple counting
   // algorithm.
   
   // This will become a look up for each frequent character.
   
   
   int block_size = 10;
   int block      = 0;
      
   // Lookup table to give the position of the symbol lookup row 
   // in the look up matrix.
   int LOOKUP[ALPHABET_SIZE];   
   
   int *symbol_lookup = malloc(sizeof(int)* threshold *(block_size+1));
   
   for (int i=0; i<ALPHABET_SIZE; i++)
      LOOKUP[i]=-1;
   
   for (int i=0; i<ALPHABET_SIZE; i++)
   { 
      if (lookup[i] < threshold && lookup[i] > 0)
      {
         printf("Method 2\n");
         printf("%c\n", i);

         LOOKUP[(unsigned char)i] = block*threshold;
                     
         // Create lookup for this symbol, and store it in the look up 
         // matrix.
         createLookup(symbol_lookup + block*threshold ,i,pattern, m, threshold);
         
         printf("Creating block: %d\n", block);
         block ++;
      }   
         
      // Perform marking when we have filled up the lookup matrix,
      // or when we are at the end of the alphabet)
      if (block >= block_size || i == ALPHABET_SIZE-1)
      {
         printf("REACHED END OF BLOCK SIZE: %d\n", block);
               
         markMatches(LOOKUP, symbol_lookup, threshold, matches, text, n,m);
      
         // Zero the look-up table.
         for (int i=0; i<ALPHABET_SIZE; i++)
            LOOKUP[i]=-1;
         
         // reset the block value.
         block = 0;
      }
         
   }
   
   // We have calculuated the number of matches.
   // subtract this from m to get the number of mismatches.
   for (int i=0; i<n-m+1; i++)
      matches[i] = m - matches[i]-1;
   
}                                

/******************************************************************************/

void display_pRepresentation(    const pTriple*  P, 
                                 const char*     pattern, 
                                       int       n                             )
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

int count_naive(                 const char*     t, 
                                 const char*     p, 
                                       int       k, 
                                       int       l                             )
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

// j gives the location in the pattern, 
// the pair (x,t) give the location in the text (t is the index into the text
// where the x'th p-triple starts). 

static inline int verifyMatch(   const pTriple*  pRepresentation,
                                 const char*     text,
                                 const char*     pattern,
                                 const ESA*      esa,
                  
                                       // The position in the p-representation.              
                                       int       x,    
                                       int       t,    
                                       int       i,    
                        
                                       // problem-specific variables.
                                 const int       k,
                                 const int       n,
                                 const int       m                             )
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

// This is the second case in the k-mismatches algorithm.
void k_mismatches_case2(  const char *text, 
                          const char *pattern,
                          const int  *frequency_table,
                          const int   k,
                          const int   n,
                          const int   m,
                                int  *matches                                  )
{

   // Find the positions where the pattern may match. 
   // We do this first for memory efficiency
   int sqrt_k = (int)(sqrt((double)k) + 0.5);
   
   // Create a look up matrix with 2k positions.
   int *pattern_lookup = malloc(sizeof(int)*sqrt_k*2*sqrt_k);

   // Find the first 2\sqrt{k} frequenct symbols, and mark all the positions
   // where they match.
   
   //  printf("Finding first %d characters and choosing first %d "                
   //         "instances of them in the pattern\n", 2*sqrt_k, sqrt_k);
           
           
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
       //  printf("%c is a frequent character.\n", i);
         
         // Create lookup for this symbol.
         // we store the look-ups in columns. 
         createLookup(pattern_lookup + j*sqrt_k, i, pattern, m, sqrt_k);

         // The position of the look up table for this character.
         LOOKUP[i] = j*sqrt_k;
         
         j++;
      }
   }
   
   markMatches(LOOKUP, pattern_lookup, sqrt_k, matches, text, n,m);
   
   // We no longer need this.
   free(pattern_lookup);   

   //---------------//
   
   pTriple *pRepresentation = malloc(sizeof(pTriple) * n);
   
   //pTriple *pRepresentation_old = malloc(sizeof(pTriple) * n);   

   // Construct the extended suffix array.
   
   printf("Constructing ESA and p-representation\n");
   ESA esa;   
   constructESA(pattern, m, &esa, 0);
   construct_pRepresentation(pRepresentation, text, pattern, &esa, n, m);
   
   //construct_pRepresentation_old(pRepresentation_old, text, pattern, &esa, n, m);   
   //display_pRepresentation(pRepresentation_old,     pattern, n);
   //display_pRepresentation(pRepresentation, pattern, n);
   
   /*
     printf("%d\n", n);
   for (int i=0; i<n; i++)
   {
      if (pRepresentation[i].j != pRepresentation_old[i].j)
      {
         
         printf("j FAILED HERE: %d\n",i);
         printf("found: %d, expected: %d\n", pRepresentation[i].j, pRepresentation_old[i].j);
         
         exit(0);
      }
      
      if (pRepresentation[i].l != pRepresentation_old[i].l)
      {
         printf("l FAILED HERE: %d\n",i);
         printf("found: %d, expected: %d\n", pRepresentation[i].l, pRepresentation_old[i].l);
                  printf("found: %d, expected: %d\n", pRepresentation[i].j, pRepresentation_old[i].j);
         exit(0);
      }
      
      if (pRepresentation[i].j == -2)
      {
         break;
      }
   
   }
   
  
   */
   
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
   
   int count=0;
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
      
        // printf("\nVerifying position: %d\n", i);
         
         // Verify this location    
         matches[i] = verifyMatch(pRepresentation, text, pattern, &esa, x, t, i, k, n, m);     
         ++count;
      } else matches[i] = k+1;
         
   }
   printf("Perfomed %d verifications\n", count);
   freeESA(&esa);
   free(pRepresentation);
}

/******************************************************************************/


void kmismatches(         const char *text, 
                          const char *pattern,
                          const int   k,
                          const int   n,
                          const int   m,
                                int  *matches  )
{
                           
   // zero the matches array.
   memset(matches, 0, sizeof(int)*(n-m+1));

   // This will be a count of the occurences of symbols.
   int frequency_table[ALPHABET_SIZE];
   
   // Count the number of occurences of each syumbol.
   sp_km_count_symbols(pattern, m, frequency_table);

   // Number of appearances required for a character to be classed 'frequent'.
   int sqrt_k = (int)(sqrt((double)k) );
   
   printf("Threshold: \\sqrt k %d.\n", sqrt_k);

   int num_freq_chars = count_frequent_symbols(frequency_table, sqrt_k);
   
   printf("There are %d frequent characters\n", num_freq_chars);

   // Which k-mismatches case to perform.
   if (num_freq_chars < 2* sqrt_k)
   {
      printf("CASE 1\n");
      
      /*--- Frequent characters ---*/

      // Go through every symbol, looking for frequent symbols.
      // then perform FFT matching on each one.
      for (int i=0; i<ALPHABET_SIZE; i++)
      {
         if (frequency_table[i] >= sqrt_k)
         {  
            printf("Method 1\n");      
            printf("%c\n", i);

            match_with_FFT(matches, i, text, pattern,  n, m);
         }
      }   
      
      

      /*--- Infrequent characters ---*/
      
      // Now, deal with the infrequent characters using a simple counting
      // algorithm.
          
      int block_size = 10;
      int block      = 0;
         
      // Lookup table to give the position of the symbol lookup row 
      // in the look up matrix.
      int LOOKUP[ALPHABET_SIZE];   
      
      int *symbol_lookup = malloc(sizeof(int)* sqrt_k *(block_size+1));
      
      for (int i=0; i<ALPHABET_SIZE; i++)
         LOOKUP[i]=-1;
      
      for (int i=0; i<ALPHABET_SIZE; i++)
      { 
         if (frequency_table[i] < sqrt_k && frequency_table[i] > 0)
         {
            printf("Method 2\n");
            printf("'%c'\n", i);

            LOOKUP[(unsigned char)i] = block*sqrt_k;
                        
            // Create lookup for this symbol, and store it in the look up 
            // matrix.
            createLookup(symbol_lookup + block*sqrt_k,i,pattern,m,sqrt_k);                      
            
            printf("Creating block: %d\n", block);
            block ++;
         }   
            
         // Perform marking when we have filled up the lookup matrix,
         // or when we are at the end of the alphabet)
         if (block >= block_size || (i == ALPHABET_SIZE-1 && block > 0))
         {
            printf("REACHED END OF BLOCK SIZE: %d\n", block);
                  
            markMatches(LOOKUP, symbol_lookup, sqrt_k, matches, text, n,m);
         
            // Zero the look-up table.
            for (int i=0; i<ALPHABET_SIZE; i++)
               LOOKUP[i]=-1;
            
            // reset the block value.
            block = 0;
         }           
      }
      
      printf("Matches in place: %d\n", matches[5860]);
      
      // We have calculuated the number of matches.
      // subtract this from m to get the number of mismatches.
      for (int i=0; i<n-m+1; i++)
         matches[i] = m - matches[i]-1;
         
   } else {
   
      printf("CASE 2\n");
      k_mismatches_case2(text, pattern, frequency_table, k, n, m, matches);
   
   }
}                     

/******************************************************************************/
// We give a lookup table which takes symbols in the alphabet to rows in the
// lookup_matrix. l Gives the max length of any given row (-1 terminates a row
// early). We thus go through the text marking the positions where each 
// symbol would be in relation to a start of the pattern. We do this in blocks
// like this since it reduces memory bottle necking.

void markMatches(                const int*      lookup,        
                                 const int*      lookup_matrix,
                                       int       l,
                                       int*      matches,
                                 const char*     text,
                                 const int       n,      
                                 const int       m                             )
                                       
{
   printf("Doing Marking\n");
   
    
   for (int i=0; i<n; i++)
   {
      // If this symbol is one of our look-up characters.
      // Then this will give us the offset into the lookup matrix.
      int x = lookup[(unsigned char)text[i]];
     
      if (x != -1)
      {
         // Create a pointer to the row we are interested in.
         const int *this_table = lookup_matrix + x;
         
         for (int j=0; j<l; j++)
         {
            // if this table is terminated early.
            if (this_table[j] == -1) break;
            
            // TODO: CHECK THIS IS RIGHT
            if (i-this_table[j] >= n-m+1) break;            
            
            if ( i - this_table[j] >= 0 )
               ++ matches[i-this_table[j]];
         }
      }
   }
   
   printf("Done Marking\n");
}                                        



void hamming_naive(       const char *text, 
                          const char *pattern,
                          const int   n,
                          const int   m,
                                int  *matches          )
{
   memset(matches, 0, sizeof(int)*(n-m+1));
   for (int i=0; i<n-m+1; i++)
   {
      int mismatches =0;
      for (int j=0; j<m-1; j++)
      {
         if (pattern[j] != text[i+j]) 
            mismatches ++;
      }
   
      matches[i] = mismatches;
   
   }
}

/******************************************************************************/

void kmismatches_naive(   const char *text, 
                          const char *pattern,
                          const int   k,
                          const int   n,
                          const int   m,
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
                          const int   k,
                          const int   n,
                          const int   m,
                                int  *matches          )

{   
   // Zero the matches array.
   memset(matches, 0, sizeof(int)*(n-m+1));
   
   
   // Construct the extended suffix array.
   ESA esa;   
   constructESA(pattern, m, &esa, 0);   
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
     // loadBar(i,n,100,40);
   
      // Advance through the p-reprsentation.
      if (t + pRepresentation[x].l <= i)
      {
         t += pRepresentation[x].l;
         ++x;
      }
           
      int v = verifyMatch(pRepresentation, text, pattern, &esa, x,t,i,k,n,m);

      matches[i] = m-v;
      
   }
   
   free(pRepresentation);
   freeESA(&esa);
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
