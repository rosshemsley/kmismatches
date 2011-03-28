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
         {
          //  printf("%d, %c, %c\n",i, t[i], t[i+per]);
            // The periodic stretch has ended.
            if (t[i] != t[i+per])
            {
            //   printf("Periodic stretch ends at %d\n", i);
               break;           
           }
         }
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

void displayBreaks(const char *t, int *breaks, int n, int l, int b)
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

// We seek a value of l such that there are at least 2k l-breaks, and l<k
int find_l(const char *t, int n, int k, int *bn, int *breaks)
{
   // Do a linear search for now. 
   int success=0;
   for (int l=k; l>=2; l--)
   {
      int b = partition(t, l, n, breaks);
   
      if (b> 2*k) 
      {  
         *bn = b;
         success=1;
         printf("DONE: l=%d, b=%d, k=%d\n",l,b,k);
         return l;
      }
   }
   
   if (!success)
   {
      printf("No such l exists\n");
   }
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
                                       int       bn                            )
{

   // Zero the matches array.
   memset(matches, 0, sizeof(int)*(n-m+1));

   // Construct the ESA for the text.
   // TODO: Don't bother with child values for this?
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

      int x = findSubstringPosition(thisBreak, k, 0, esa.n, &esa); 
     
      // TODO: Fix this!
      // This currently fails sometimes, there seems to be a bug in SAIS.
      assert( x>=0 );
      
      // Find all locations of this k-break and mark in the matches 
      // array the starting position.
      do
      {      
         // Text location of this match.
         int j = esa.SA[x];
         
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
      } while (x < esa.n  &&  esa.LCP[x] >= k);
      
   }  
   
   /**
   *  kangaroo accross all the potential matching positions.
   */
   //   printf("matches at right val: %d\n", matches[34871780]);
   
   for (int i=0;i<n-m+1;i++)
   {
      // If there could be a match here.
      if (matches[i] >= k-1)
      {
         printf("Verifying: %d\n", i);
         matches[i] = verify(i, n-1, m,  k, &esa);             
      } else 
         matches[i] = k+100;
   }
}

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

   // This is the largest possible value of b.
   int  pn      = m;   
   int *breaks  = calloc (pn, sizeof(int));   
      
   // Partition in the text into its l-breaks.
   pn = partition(pattern, k, m, breaks);
   
   // displayBreaks(p, breaks, m, k, pn);   
   printf("There are %d pattern breaks\n", pn);
  
   if (pn >= 2*k)
   {    
  
      printf("There are enough k-breaks\n");        
      // Only use the first 2k kbreaks for matching.
      simpleMatcher(text, pattern, breaks, matches, k, n, m, 2*k);   
 
      return 1;
   }
   else 
   {
      printf("There are insufficient k-breaks\n");
      return 0;
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
