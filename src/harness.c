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
//
// NOTE: SAIS APPEARS TO FAIL WITH INPUT 'aacdbbcca', 'dddbbcddc'
//
/******************************************************************************/

int main(int argc, char **argv)
{

   if (argc < 2)
   {
      fprintf(stderr, "No input file provided\n");
      exit(1);
   }
   
   int naive          = 0;
   int naive_hamming  = 0;
   int _kangaroo      = 0;
   int abrahamson     = 0;
   int bs_abrahamson  = 0;
   
   
   if (argc == 3)
   {
      for (int i=2; i<argc; i++)
      {
         if (strcmp(argv[i], "-naive") == 0)
         {
            printf("Using naive k-mismatches algorithm\n");
            naive=1;
         }
            else if(strcmp(argv[i], "-naive_nm") == 0)
         {
            printf("Using Naive Hamming distance algorithm\n");
            naive_hamming = 1;  
         }
            else if(strcmp(argv[i], "-kangaroo") == 0)
         {
            printf("Using Kangarooing\n");
            _kangaroo = 1;  
         } 
            else if (strcmp(argv[i], "-abrahamson") == 0)
         {
            printf("Using Abrahamson/Kosaraju\n");
            abrahamson=1;
         } 
            else if (strcmp(argv[i], "-bs_abrahamson") == 0)
         {
            printf("Using Ben Smither's Abrahamson/Kosaraju\n");
            bs_abrahamson = 1;
         }
            else 
         {
            printf("Invalid arguments. Exiting.\n");
            exit(1);
         } 
      }   
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


   // An array to put the matching locations in.
   int  *matches = malloc(sizeof(int)  * (n-m+1));


   // Perform matching.
   if (naive)
      kmismatches_naive(t,p,k,n,m,matches);
      
   else if(naive_hamming)
      hamming_naive(t,p,n,m,matches);
      
   else if(_kangaroo)
      kangaroo(t,p,k,n,m,matches);
      
   else if(abrahamson)
      abrahamson_kosaraju(t,p,n,m,matches);
      
   else if(bs_abrahamson)
   {
      struct SP_KM_MATCHING_POSITIONS *listOfMatches;
      listOfMatches = sp_km_create_new_list_of_matches();
      
      int numMatches = 0;

	   sp_km_unbounded_kmismatch(t,p,n,m,k, &numMatches,listOfMatches, 0);
	   
	   // Don't attempt to test this method.
      exit(0);
   }   
   
   else
      kmismatches(t,p,k,n,m,matches);
      
     
   // Verify the output.   
   int pass=0;
   
   for (int b=0;b<n-m+1; b++)
   {
      if (matches[b] <=k) 
      {
         printf("Position: %d, mismatches: %d\n", b, matches[b]);  
         if (b == pos && matches[b] == k)
         {
            printf("PASSED TEST\n");
            pass = 1;
         }
      }  
   }
   
   
   if (!pass)
   {
      fprintf(stderr, "FAILED TEST\n");
      exit(1);
   }   
  
  
   free(p);
   free(t);
   free(matches);   
   // This needs to be fixed.
   // FreeRMQ_succinct();   
   
   
   return 0;
   
}

/******************************************************************************/

