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
#include "RMQ_succinct.h"
#include "./sp_km_unbounded_matcher.h"

/******************************************************************************/
// Get the period of this block of length n.
// The period can clearly be at most n/2.

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

int partition(const char *t, int k, int n, int *kbreaks)
{
   // There can be at most n/k distinct k-breaks.



   int current=0;
   
   for (int i=0; i<n; i++)
   {
     
      // This is the starting point of this periodic stretch.      
      int start = i;
      
      // Determine Period.
      int per = getPeriod(t+i, k);
      
      printf("Found period %d\n", per);
      
      // Jump forwards by the length of the period. 
      // Since we know that the first cycle matches.    
      i += per;
      
      // If the period is non-zero, walk through this periodic
      // stretch for as long as possible.
      if (per!=0)
      {
         for(;i<n;i++)
         {
            printf("%d, %c, %c\n",i, t[i], t[i+per]);
            // The periodic stretch has ended.
            if (t[i] != t[i+per])
            {
               printf("Periodic stretch ends at %d\n", i);
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
         if (i-k < start)
         {
            kbreaks[current] = start;
            i = start+k-1;
         } else {
            kbreaks[current] = i-k;
            i--;
         }   
         
            current++;
      }
          
   }   
}

/******************************************************************************/

// We seek a value of l such that there are at least 2k l-breaks, and l<k
void find_l(const char *t, int n, int k, int *breaks)
{
   // Do a linear search for now. 
   
   for (int l=0; l<k; l++)
   {
//      partition(
   }
}

/******************************************************************************/
// b is the maximum number of breaks.
void displayBreaks(char *t, int *breaks, int n, int b)
{
   int z=0;
      
   for (int i=0;i<n;i++)
   {
      if (breaks[z] == i)
      {
         printf("[");
         int x = i+k;
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


int main(int argc, char **argv)
{
   
   char t[]="abcabcabcdabcabcabcabcfabcabcababababababababfabababababfababk";
//   char t[]="asdfabababababasdfsd";
  
   
   int k=7;
   int n=strlen(t);

   int  b      = n/k+1;   
   int *breaks = calloc (b, sizeof(int));   
   
   partition(t,k,n,breaks);

   for (int i=0; i<b; i++)
   {
      printf("%d\n", breaks[i]);
   } 
  printf("%s\n\n", t);
  
  displayBreaks(t, breaks, n, b);
  

}


