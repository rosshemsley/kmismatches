#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <time.h>
#include <fftw3.h>
#include <assert.h>
#include <string.h>
#include "km.h"
#include "sais.h"
#include "stack.h"
#include "loadTest.h"
#include "RMQ_succinct.h"
#include "./sp_km_unbounded_matcher.h"

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
// return 1 if a > b, 0 if a==b and -1 if a < b

int str_gth(const char *a, const char *b, int n)
{
   
   
   

   for (int i=0; i<n; i++)
   {
      
      if (a[i]=='\0') return 0;
      if (b[i]=='\0') return 1;

      printf("Comparing: %c to %c\n", a[i], b[i]);
      
      if ((unsigned char)a[i] > (unsigned char)b[i]) return 1;
      if ((unsigned char)a[i] < (unsigned char)b[i]) return -1;
   }
   
   return 0;
}

/******************************************************************************/
// Find the first location of a substring in O(n\log m) [O(n + \log m) 
// expected time].
// This algorithm comes from Gusfield.

int findSubstring(const char *p, const char *t, const int *SA, int n)
{
   int min = 0;
   int max = n-1;
   int mid;
   
   do 
   {
      mid   = min+(max-min)/2;      
      int c = str_gth(p, t + SA[mid], n);    
      
      
      printf("min, max: %d, %d, %d\n", min, max,mid);
      if (c == 1)
         min = mid+1;
         
      else if (c == -1)
         max = mid-1;
         
      else if (c == 0)
         return mid;
         
   } while (min < max);

   // Could not find the substring.
   return -1;
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
/******************************************************************************/

int main(int argc, char **argv)
{
   
   char *t="helloaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaadsfasdfasdfasdfasdfasdfasdfasdgsdfhg";
   char *p="helloaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";  
   
   int k = 2;
   int n = strlen(t);
   int m = strlen(p);

   int pos=0;

   if (argc > 1)
      load(argv[1], &n, &m, &k, &pos, &t, &p);
    
    printf("Text length:%d, pat length: %d\n", n,m);

   printf("%s\n", t);

    p="accidentally in love";  

    m = strlen(p);

   // This is the largest possible value of b.
   int  pn      = m/k+1;   
   int *breaks = calloc (pn, sizeof(int));   
      
   // Partition in the text into its l-breaks.
   pn = partition(p, k, m, breaks);
   
  // displayBreaks(p, breaks, m, k, pn);
   
   printf("There are %d pattern breaks\n", pn);


   int *SA = malloc(n*sizeof(int)+2);
   int *LCP = malloc(n*sizeof(int)+2);
   sais((unsigned char*)t, SA,LCP, n);

   int x = findSubstring(p, t, SA, n);
   

   
   printf("Found substring: %d\n",x);
  
//   int * matches = calloc(n-m+1, sizeof(int));
   
//   match(t, p, pbreaks, tbreaks, k, n, m, pn, tn, matches);
    //  printf(" ");
   printf("\n");
}

/******************************************************************************/













