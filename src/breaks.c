#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <time.h>
#include <fftw3.h>
#include <assert.h>
#include <string.h>
#include "km.h"
#include "esa.h"
#include "stack.h"
#include "loadTest.h"
#include "RMQ_succinct.h"
#include "./sp_km_unbounded_matcher.h"

/******************************************************************************/


/******************************************************************************/
// Get the period of this block of length n.
// The period can clearly be at most n/2.

// TODO: Do this better!
// This is just a naive O(n^2) algorithm.
int getPeriod(const char *t, int n)
{

   for (int i=1; i<n/2+1;i++)
   {
      int match=1;
      for (int j=0; j<i; j++)
      {        
         // TODO: stop overflow
         if (t[j] != t[j+i])
         {
            // This cannot be the period.
            match=0; 
            break; 
         }
      }
      
      if (match)
         return i;
   }

   return 0;
}

/******************************************************************************/
// Partition into k-breaks.

int partition(const char *t, int l, int n, int *breaks)
{

   int current=0;
   
   for (int i=0; i<n; i++)
   {
     
      // This is the starting point of this periodic stretch.      
      int start = i;
      
      // Determine Period.
      int per = getPeriod(t+i, l);
      
      //printf("Found period %d\n", per);
      
      // Jump forwards by the length of the period. 
      // Since we know that the first cycle matches.    
      i += per;
      
      // If the period is non-zero, walk through this periodic
      // stretch for as long as possible.
      if (per!=0)
      {
         for(;i<n;i++)            
            // The periodic stretch has ended.                     
            if ((i+per < n) && t[i] != t[i+per])
               break;           
         
         // This is the position of the mismatch.
         i = i+per+1;         
      }
      
      // i contains the position of the msmatch.
      // now, we cannot form a break which cuts of the
      // periodic stretch too early: there must be at least i-start
      // characters in the periodic stretch.
      // We can now construct a new k-break.
      if (i<n)
      {
         if (i-l < start)
         {
            breaks[current] = start;
            i = start+l-1;
         } else {
            breaks[current] = i-l;
            i--;
         }   
         
         current++;
      }
          
   }   
   
   
   // Check that this doesn't go out of bounds etc.
   if (current < n/l)
      breaks[current] = -1;
   return current;
}



/******************************************************************************/

void displaySubStr(const char *t, int _l, int _r, int n, int l, const int *breaks, int b)
{
   int z=0;
      
   for (int i=0;i<n && i<_r; i++)
   {
      if (breaks[z] == i)
      {
         if (i>=_l && i <_r)
            printf("[");
         int x = i+l;
         for (;i<x && i<n; i++)
            if (i>=_l && i <_r)
               printf("%c", t[i]);
         if (i>=_l && i <_r)
            printf("]");
         z++;
         i--;
      } else {      
         if (i>=_l && i <_r)
            printf("%c", t[i]);      
      }
   }
}

/******************************************************************************/
// b is the maximum number of breaks.

void displayBreaks(const char *t, const int *breaks, int n, int l, int b)
{
   int z=0;
      
   for (int i=0;i<n;i++)
   {
      if (breaks[z] == i)
      {
         printf("[");
         int x = i+l;
         for (;i<x && i<n; i++)
            printf("%c", t[i]);
         printf("]");
         z++;
         i--;
      } else {      
         printf("%c", t[i]);      
      }
   }
   printf("\n");

}

/******************************************************************************/

// We seek a value of l and m such that there are at least 2k l-breaks, and 
// fewer than 2k m-breaks.
int find_l(const char *t, int n, int k, int *ln, int *mn, int *lbreaks, int *mbreaks)
{
   // Do a linear search for now. 
   for (int l=k; l>=2; l--)
   {
      *ln = partition(t, l-1, n, lbreaks);
      *mn = partition(t, l,   n, mbreaks);
   
      if (*ln >= 2*k && *mn <= 2*k) 
      {  
         printf("There are %d l-breaks\n", *ln);
         printf("There are %d m-breaks\n", *mn);         
         printf("DONE: l=%d, k=%d\n",l-1,k);
         return l-1;
      }
   }
   
   printf("No such l exists\n");   
   return -1;
}

/******************************************************************************/
// Simple Kangarooing, for when we have a full ESA for the text and
// not just a p-representation.

// Calcululate the number of mismatches between the substrings
// starting at i and j respectively. 
// If there are more than k, return.

// TODO: CHECK END CONDITIONS.

int verify(int i, int j, int m, int k, const ESA* esa)
{
   
   // The number of mismatches.
   int mismatches=0;
   
   // The position in the pattern.
   int end = j+m;
   
   int length=0;
   
   while (j < end)
   {
      
     // printf("Finding longest extension\n");
      // The longest number of shared characters.
      int l = LCE(i, j, esa);
      
      length+=l;
      
      //printf(" found %d matching chars\n", l);
      i += l+1;
      j += l+1;
      
      mismatches ++;
      
      if (mismatches > k)
      {
         return m-length-1;
      }  
   }
   
   
   return m-length-1;

}

/******************************************************************************/
// When there are at least 2k k breaks, we can do the following in O(n)

void simpleMatcher(              const char*     text,
                                 const char*     pattern,
                                 const int*      kbreaks,
                                       int*      matches,
                                       int       k,
                                       int       n,
                                       int       m,
                                       int       bn,
                                 const ESA*      esa   )
{

   // Zero the matches array.



   
   
   // Go through all of the k-breaks, and mark the starting positions.
   
   
   printf("bn is: %d, k is; %d\n", bn, k);
   
   // Loop through all of the k-breaks.
   for (int i=0; i<bn; i++)
   {
   

      // The k-break we are currently considering.
      const char *thisBreak = pattern + kbreaks[i];
      
      // Find the first location of this break  in the text (if applicable)
      // in the suffix array.
      // TODO: Add a lookup table for first level.

      //  printf("Next break \n");

      int x = findSubstringPosition(thisBreak, k, 0, esa->n, esa); 
     
      // TODO: Fix this!
      // This currently fails sometimes, there seems to be a bug in SAIS.
      assert( x>=0 );
      
      // Find all locations of this k-break and mark in the matches 
      // array the starting position.
      do
      {      
         // Text location of this match.
         int j = esa->SA[x];
         
         // Make sure this suffix comes from the text and not the pattern.
         if (j<n-1 )
         {         
         // If this does not run off the end of the matches array, then
         // mark in the matches array the possible starting position of the 
         // pattern.
            if (j-kbreaks[i] >= 0)
            {
               ++matches[j-kbreaks[i]];
             //  printf("Marking %d\n", j-kbreaks[i]);   
            }
         } //else // printf("Not in text\n"); 
         ++x;
      } while (x < esa->n  &&  esa->LCP[x] >= k);
      
   }  
   
   /**
   *  kangaroo accross all the potential matching positions.
   */
   //   printf("matches at right val: %d\n", matches[348717i80]);
   
   for (int i=0;i<n-m+1;i++)
   {
      // If there could be a match here.
      if (matches[i] >= k)
      {
         printf("Verifying: %d\n", i);
         matches[i] = verify(i, n-1, m,  k, esa);             
      } else 
         matches[i] = k+1;
   }
}


/******************************************************************************/
// Return  0 for equality.
// return -1 for less than.
// Return  1 for greater than.

int set_cmp(int a, int b, int l)
{
   if (a < b   ) return -1;  
   if (a >= b+l ) return  1;
   
   return 0;
  
}

/******************************************************************************/
// Perform a binary search for this break.
// Given a region of length from {x,...,x+l}, find any value in lookup
// which is contained within that region.

// There can be at most two matches in this special case we are considering.
// Thus we look at the previous value, and if it also satisfies the critiera, 
// we return it. In this way, we always know we have the first.
int binaryBreakSearch(int x, int l, const int *lookup, int n)
{
   int mid;
   int min = 0;
   int max = n-1;
   
   do
   {
      mid = min + (max-min)/2;     
     
      int e = set_cmp(lookup[mid], x,l);
      
      if (e > 0)
         max = mid - 1;
      else if (e < 0)
         min = mid + 1;
      else 
      {
         // We always want to return the first possible match.
         if (mid>0 && set_cmp(lookup[mid-1],x,l) == 0)
            return mid-1;
         else 
            return mid;         
      }   
      
   } while ( min <= max );

   return -1;
  
}

/******************************************************************************/


// Perform matching in O(k log k) for contiguous blocks of length l.
int algorithm_2(                       int       x,
                                       int       l,
                                       int       n,
                                       int       m,
                                       int       k,
                                       int       ln,
                                 const ESA*      esa,
                                 const int*      lbreaks,
                                 const int*      dbreaks,
                                 const int*      lookup,
                                 const int*      indicies,                      
                                       int*      matches                       )
{
   // Perform Marking //   
   
   // For each of the first 2k breaks.
   for (int i=0; i<2*k; i++)
   {
      // This is the correct lookup offset for this break (I hope).
      const int* breakArr = indicies + (n/k) * dbreaks[i];   
      int        block    = (x + lbreaks[i])/k;
      
      // Figure out which of the disjoint breaks we are looking at.
      
      // This is the index of the distinc break:
      printf("Looking for break: (lbreak) %d (dbreak) %d\n", lbreaks[i], dbreaks[i]);
      // This is the index of the current pattern-break.
      printf("Looking in block: %d\n", block);
      
            
      // This is the block index, it gives the index of the block
      // that this break would have to fall into if the pattern were lined
      // up to start at x.     
      
      // Have a peek at the block.
      for (int z=0; z<k; z++)
         printf("%d\n", breakArr[ z ]);
            
      // This is the index into the lookup array of the start of this sorted
      // array.      
      int start = breakArr[ block ];
      
      // For now, assume that we go over a break boundary.
      // This probably doesn't add any overhead.
      

      int end;
      // The following gives us the end point:      
      // TODO: CHECK THE LOGIC HERE
      if (block + 2 <= ln)
         end = breakArr[ block+2 ];
      else 
         end = breakArr[ block+1 ];
         
      printf("Start: %d\n", start);
      printf("End:   %d\n", end);
      printf("smallest: %d, largest: %d\n" ,x+lbreaks[i], x+lbreaks[i]+l-1);


      printf("Doing search on the following:\n");
      
      for (int y=start; y<end; y++)
      {
         printf("%d\n", (lookup)[y]);
      }

      int f = binaryBreakSearch(x + lbreaks[i], l, lookup + start, lookup[end]-lookup[start]) + start;
      
      
      printf("Found: %d (%d)\n", f, lookup[f]);
      matches[lookup[f]-lbreaks[i]] += 1;
      
      if ( set_cmp(lookup[f+1], x+lbreaks[i], l) == 0 )
      {
         printf("Found: %d (%d)\n", f+1, lookup[f+1]);
         matches[lookup[f+1]-lbreaks[i]] += 1;
      }
      
      // arr now contains all the pointers to the starts of each instance of this
      // break for each block of length k
      
      // We find out which block we are intersted in: it is (x+breaks[i])/k, 
      //((x+breaks[i])/k +1)
      
      
      printf("\n");               
     
      
   }
      
   // Perform Verification on the, at most 4 locations with >=k marks //
   
   
   for (int i=x; i< x+l-1; i++)
   {
      // If there could be a match here.
      if (matches[i] >= k)
      {
         printf("Verifying: %d\n", i);
         matches[i] = verify(i, n-1, m,  k, esa);             
      } else 
         matches[i] = k+1;
   }
   

   return 0;
}

/******************************************************************************/
// TODO: Think carefully about unsigned types.

// Create n/k sets of 2k pointers to sorted arrays of distinct l-breaks.
int constructLookups(            const int*      breaks, 
                                       int       bn,                                  
                                 const char*     text,
                                 const char*     pattern, 
                                 const ESA*      esa, 
                                       int       l, 
                                       int       k, 
                                       int       n, 
                                       int       m,
                                       int*      dbreaks,
                                       int*      lookup, 
                                       int*      indicies                      )
{
   int *breakPositions = malloc(sizeof(int)*n);
   int *breakCounts    = calloc(2*k, sizeof(int));
  // int *dbreaks = malloc(sizeof(int)*bn);
   
   for (int i=0; i<n; i++)
      breakPositions[i] = -1;
   
   // This counts the number of disjoint breaks.
   int count = -1;
   
   for (int i=0; i<2*k; i++)
   {   
      // printf("Looking for break: %d\n", breaks[i]); 
      // Find the first instance of this substring.    
      int j = findSubstringPosition(pattern + breaks[i], l, 0, esa->n, esa); 

      printf("Location: %d, string; %.10s\n", j, pattern + breaks[i]);

      // Assume the esa is a generalised suffix tree, and so we always find a 
      // match.      
      assert(j >= 0);
      
      int non_dup = 0;
      do {
         int x = esa->SA[j];
       
         printf("Match: %d\n", x);
         if (x<n-1)
         {

            // Duplicate break.
            // Label it in the dbreaks array and continue.
            if (breakPositions[x] != -1)
            {
               dbreaks[i] = breakPositions[x];
               break;
            }
            
            if (non_dup == 0)
            {
               count++;
               printf("last count:%d\n", breakCounts[count]);
               // This is a new break.
               dbreaks[i]            = count;
               breakPositions[x]     = count;      
               breakCounts[count]    = 1;         
               non_dup = 1;
           } else
           {
               x = esa->SA[j];
               breakPositions[x]  = count;
               breakCounts[count] ++;
            }
         }
         ++j;
      } while (j < esa->n  &&  esa->LCP[j] >= l);
   }
      
   count +=1;
  

   printf("First 2k distinct breaks:\n");
   for(int i=0; i<count; i++)
      printf("%d, %.2s\n", dbreaks[i], pattern + dbreaks[i]);
   
   printf("Break locations:\n");
   for (int i=0; i<n; i++)
      if (breakPositions[i] <0)
         printf(" ");
      else printf("%d", breakPositions[i]);
  printf("\n");
 
   printf("Counts: \n");
   for (int i=0; i<count; i++)
   {
      printf("%d\n", breakCounts[i]);
   }
 
   printf("\n");
   
   // breakPositions now contains a look up for all the breaks.
   
   // Create an array which will point to the start positions of the
   // disjoint breaks in the lookup array we will create shortly.
   // it stores the last value too, to make end caes easier.
   int *breakIndicies = malloc(sizeof(int)*(count+1));
   
   // Now, set those break indicies to point to placeholders in the
   // lookup array.
   
   // The first break starts at index 0.
   // TODO: think about end case.
   breakIndicies[0]=0;
   for (int i=1; i<count+1; i++)
      breakIndicies[i] = breakIndicies[i-1] + breakCounts[i-1];

   printf("Break indicies: \n");
   
   for (int i=0; i<count; i++)
      printf("%d\n", breakIndicies[i]);

   // Allocate just enough space in the lookup table to contain all of
   // the matches for the breaks.
   //lookup = malloc(sizeof(int)*  breakIndicies[count]);
   
   printf("Total number of disjoint breaks: %d\n", breakIndicies[count]);
   
   // We now copy the breaks into the lookup table, keeping them sorted
   // by maintaining pointers into the correct positions in the array.
   int *temp_pointers = malloc(sizeof(int) * count);
   
   // The temporary pointers point to the last location inserted in the 
   // lookup for each break.
   for (int i=0; i<count; i++)
      temp_pointers[i] = breakIndicies[i];
   
   // Go through the break positions array, and put all the instances
   // of each break in the correct place in the lookup.
   for (int i=0; i<n; i++)
   { 
      if(breakPositions[i]<0) continue;
     // printf("  break position: %d\n", breakPositions[i]);
      //printf("  pat position:   %d\n",   dbreaks[breakPositions[i]]);
      lookup[temp_pointers[breakPositions[i]]] = i;
      temp_pointers[breakPositions[i]]++;
   }

   printf("Lookup: (length: %d)\n", breakIndicies[count]);
   
   for (int i=0; i<count; i++)
   {
      printf("  Next: \n");
      for (int j=breakIndicies[i]; j<breakIndicies[i+1]; j++)
      {
         printf("%d\n", lookup[j]);
      } 
      
   }
   
   // We now have a lookup containing all the positions of the disjoint breaks,
   // sorted first by break index, and then by pattern index.
   // We now require n/k sets of 2k (>=count) pointers to index this array in 
   // order to perform "algorithm-2" efficiently
   
   // No longer needed.
   free(breakPositions);
   free(temp_pointers);
      
//   int partitions = (int)( (float)n/(k+0.5) );


      
     
   // We now set the pointers to give the locations of the start points for
   // each set of breaks within each of the n/k blocks of length k

   // For each disjoint break, go through all the locations setting pointers
   // to the start and end of the set of breaks within the blocks of length k.
   for (int i=0; i<count; i++)
   {
      // This is the sorted array of all the i'th breaks in the text.
      int *breakArr = lookup + breakIndicies[i];
      int  length   = breakCounts[i];
      
      int *indiciesArr = indicies + i * (n/k);
      
      // We go through all of the break instances, and set pointers every time
      // we go over a block of length k.
      int boundary = 1;
      
      // The first boundary is at 0.
      indiciesArr[0]  = breakIndicies[i];
      for (int j=0; j<length; j++)
      {
         //printf("Looking at: %d:%d\n",j, breakArr[j]);
         
         // This is the block within this break should be:
         int x = breakArr[j]/k;
         
         //printf("Fits into block: %d\n", x);
         
         // If we have moved to a new boundary.
         if (x>=boundary)
         {
            // It may be that we have 'skipped a few' of the i's.
            // We must therefore copy back the correct end value in order
            // for the algorithm to work properly.
            for (int k=boundary; k<=x; k++)
            {
               //printf("Setting index: %d to %d\n", k, j);
               indiciesArr[k] = j + breakIndicies[i];         
            }  
            boundary = x+1;
         }
         
      }
      

   }
      // Test output //
      for (int i=0; i<n; i++)      
        printf("%-3d Boundary: %d, k-value: %d\n", i, indicies[i], (i%(n/k))*k);
   // Return the number of disjoint breaks found.
   return count;
}


/******************************************************************************/
/*
// Find the matches in a region of length l.
int algorithm_2(int x0, char* text, char *pattern, int breaks, int bn, int l, ESA *esa, int k)
{
   
   
   return 0;
}
*/

/******************************************************************************/
// Try to use the periodicity properties of the pattern to match.
// If the pattern is not sufficiently aperiodic (or k is too large)
// then we return 0. othewrise we return 1 to indicate success.

int periodicMatching(            const char*     text, 
                                 const char*     pattern,
                                       int       k,
                                       int       n,
                                       int       m,
                                       int*      matches                       )
{

   memset(matches, 0, sizeof(int)*(n-m+1));

   // This is the largest possible value of b.
   int  pn      = m;   
   int *breaks  = calloc(pn, sizeof(int));   
   
   // Construct the ESA for the text.
   ESA esa;   
   
   //   printf("%.10s\n", text + n+m-1 );
   // We need a generalised suffix array: 
   // Do this by using an auxillary array called tp (text-pattern)

   char *tp = malloc( sizeof(char) * ( n+m-1 ) );

   // Copy the text (minus the '\0') into tp,
   // and the pattern in after the text.
   memcpy(tp,       text,    sizeof(char) * (n-1) );   
   memcpy(tp + n-1, pattern, sizeof(char) *  m    );
  
   // Construct the suffix array.
   constructESA(tp, n+m-1, &esa, NO_CHILD_TAB);  

   // Partition in the text into its l-breaks.
   pn = partition(pattern, k, m, breaks);
   
   // displayBreaks(p, breaks, m, k, pn);   
   printf("There are %d pattern breaks\n", pn);
  
   if (pn >= 2*k)
   {    
  
      printf("There are enough k-breaks\n");        
      // Only use the first 2k kbreaks for matching.
      simpleMatcher(text, pattern, breaks, matches, k, n, m, 2*k, &esa);   
 
      return 1;
   }
   else 
   {
      printf("Constructing lookup-tables for Algorithm 2\n");
      // ** Initialise the structures for algorithm 2 **
     
      // 1) Find l-boundary
      int  ln;
      int  mn;
      int *lbreaks = malloc(sizeof(n)*n);
      int *mbreaks = malloc(sizeof(n)*n);
      
      printf("Finding l-boundary\n");      
      int l = find_l(pattern, m, k, &ln, &mn, lbreaks, mbreaks);
//      int _m = l+1; 
      
      
   
      // 2) Create look-up structure.
   
      int *lookup   =  malloc(sizeof(int)* n);
      int *indicies =  malloc(sizeof(int)* 2*n);
      int *dbreaks  =  malloc(sizeof(int)* 2*k);
      
      printf("Constructing look-ups\n");

      displayBreaks(pattern, mbreaks, m, l, ln);
      //  displaySA(&esa);
      
      
      constructLookups(lbreaks, ln, text, pattern, &esa, l, k, n, m, dbreaks, lookup, indicies);
   
   
      for (int i=0; i<n; i+=l)
      {
         algorithm_2(i, l, n, m, k, ln, &esa, lbreaks, dbreaks, lookup, indicies, matches);
      }

      return 1;
   }
   
}


/******************************************************************************/
/*
void match(char *t, char *p, int *pbreaks, int *tbreaks, int k, int n, int m, int pn, int tn, int *matches)
{
   
   // This index last used in the tbreaks array, so that 
   // we can perform O(m) look-ups instead of O(n).
   int b_pos = 0;
   
   // These will give the leftmost and rightmost breaks for the set X.
   int left, right;
   
   // For each block of length 2m in the text.
   for (int i=0; i<n; i+=m)
   {
      // u and v are the start and end of this interval.
      int u = i;
      int v = i+2*m;
   
      // find the special set X in the text which must contain all the matches.
      
      // Find the middle break.
      int middle=0;
      while (middle < tn && tbreaks[middle] < u+m){ middle ++; }
            
      if (middle-3*k < 0 || tbreaks[middle-3*k] < u) 
      left = i;      
      else left = tbreaks[middle-3*k];
      
      if (middle+3*k >= tn || tbreaks[middle+3*k] > v) 
      right = i+2*m;
      else right = tbreaks[middle+3*k];
      
      // These are the start and end values in the pbreaks array.
      // They give the start and end positions of the breaks
      // which we consider in the breaks for the text.
      int X_start=0, X_end=0;
     
      // Find the starting position of X in the breaks array.
      while (tbreaks[X_start] < left) X_start ++;
      while (tbreaks[X_end]   < right) X_end ++;         
         
      // Mark all the possible matching locations: at most 24^3 marks
      // Go through every pbreak, tbreak pairing, marking the 
      // possible starting positions.
      for (int x=X_start; x<X_end; x++)
      {
         for (int y=0; y<pn; y++)
         {
           // printf("Matching block starting at %d with (p)block starting at %d\n", tbreaks[x], pbreaks[y]);
            // The start position for this break
            // given that the breaks overlap at the last character
            // of the pattern break.           
            int start = tbreaks[x] - pbreaks[y]-k+1;
            int end   = start+2*k-1;
            if (end   <= 0)  continue;

            if (start < 0)      start = 0;
            
            if (end   > n-m+1)  end   = n-m+1;
         
                 //       printf("   start: %d, end: %d\n", start, end);
         
            // Mark all the possible starting locations.
            for (int z=start; z<end; z++)
               ++matches[z];
         
         }         
      }
            

         //  break;
   } 
   
   // Verify the locations: at most 24k^2.
}

******************************************************************************

void match2(const char *t, int n, int m, int l, int k, int b, int *breaks)
{
   // This is the current position in the breaks array.
   // this means we can look up breaks in O(m) time instead of O(n).
   int pos=0;

   // Go through blocks of the text of size 2m.
   for (int i=0; i<n; i+=2*m)
   {
      printf("Working with:\n");         
      displaySubStr(t,i, i+2*m, n, l, breaks, b);
      printf("\n");      
   
      // Seek the index of the last break starting before the middle of the 
      // current block of the text.
      int middle=0;
      int left, right;
      while (middle < b && breaks[middle] < i+m){ middle++; }
      printf("middle: %d, chars: %d\n", middle, breaks[middle]);



      // Find the index 'left'
      if (middle-3*k < 0 || breaks[middle-3*k] < i) 
      left = i;
      else left = breaks[middle-3*k];
      
      printf("left: %d\n",left);
        
      printf("%d\n", breaks[middle]+3*k);

      printf("breaks: %d\n", middle+3*k);
      // Find the index 'right'
      if (middle+3*k >=b || breaks[middle+3*k] > i+2*m) 
      right = i+2*m;
      else right = breaks[middle+3*k];
      
      printf("right: %d\n",right);
      
      
      displaySubStr(t, i, left, n, l, breaks, b);
      printf("|");
      displaySubStr(t, left, right, n, l, breaks, b);
      printf("|");
      displaySubStr(t,right, i+2*m, n, l, breaks, b);
      

      printf("\n");
      exit(0);
   
   }
   
   
}
*/

