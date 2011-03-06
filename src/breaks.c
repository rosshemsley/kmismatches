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
     
      // Determine Period.
      int per = getPeriod(t+i, k);
      
      printf("Found period %d\n", per);
      
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
      
      // We can now construct a new k-break.
      if (i<n)
      {
         kbreaks[current] = i-k;
         current++;
      }
          
   }   
}

/******************************************************************************/


int main(int argc, char **argv)
{
   
   char t[]="abcabcabcdabcabcabcabcfabcabcababababababababfababababab";


   
   int k=7;
   int n=strlen(t);
   
      int *kbreaks = calloc (  n/k+1,sizeof(int));   
   
   partition(t,k,n,kbreaks);
   
   for (int i=0; i<n/k+1; i++)
   {
      printf("%d\n", kbreaks[i]);
   }
  
   int z=0;
   printf("%s\n\n", t);
   for (int i=0;i<n;i++)
   {
      if (kbreaks[z] ==i)
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
