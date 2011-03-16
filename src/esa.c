#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <string.h>
#include "sais.h"
#include "stack.h"
#include "RMQ_succinct.h"
#include "esa.h"


/******************************************************************************/

void displaySA(                  const ESA*      esa, 
                                 const char*     pattern, 
                                       int       m                             )
{
   printf("|i  |SA |LCP|U  |D  |A  |\n");
   for (int i=0; i<m; i++)
      printf("|%3d|%3d|%3d|%3d|%3d|%3d|\n", i, esa->SA[i], esa->LCP[i],                                     
                                    esa->up[i], esa->down[i], esa->accross[i] );
      
   printf("\n\n");
   
   for (int i=0; i<m;i++)
      printf("%3d: %d %s\n", i, esa->LCP[i], (pattern+esa->SA[i]));
}

/******************************************************************************/
// Find the Longest Common Extension using the Extended Suffix Array.

inline int LCE(                        int       i, 
                                       int       j, 
                                 const ESA*      esa                           )
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
          //*/  query_naive( a+1, b, esa->LCP, esa->n );
                                  
   return     esa->LCP[temp];
               
}

/******************************************************************************/
// Construct the child_table for the ESA. TODO: Compress this down to one field
// using the optimised outlined in the relevant paper.
// 
// These algorithms are copied verbatim from the paper "Replacing the suffix 
// tree with the enhanced suffix array".

void constructChildValues(ESA *esa)
{

   int n = esa->n+1;
   
   stack *s = newStack();
   
   push(s, 0);
   
   // TODO: Make sure that this correctly reaches the end.
   for (int i=1; i<n; i++)
   {   
      while (esa->LCP[i] < esa->LCP[ peek(s) ])
         pop(s);
         
      if (esa->LCP[i] == esa->LCP[ peek(s) ])
         esa->accross[pop(s)] = i;
   
      push(s, i);      
   }
   
   /**   Construct Up/Down values.   ***************/
   
   // Reset the stack.   
   emptyStack(s);
   
   int lastIndex = -1;
   push(s, 0);
   for (int i=1; i<n; i++)
   {
      while (esa->LCP[i] < esa->LCP[ peek(s) ] )
      {
         lastIndex = pop(s);
         int top   = peek(s);
         
         if (    esa->LCP[i]   <= esa->LCP[top] 
              && esa->LCP[top] != esa->LCP[lastIndex] )
              
            esa->down[top] = lastIndex;
         
         if (lastIndex != -1)
         {
            esa->up[i] = lastIndex;
            lastIndex  = -1;
         }
      }     
      push(s, i);
   }  
   
   freeStack(s);
}

/******************************************************************************/

// Construct an extended suffix array for some string of length n.s
void constructESA(const char *s, int n, ESA *esa, ESA_FLAGS flags)
{
   esa->t     = s;
   esa->n     = n;
   esa->flags = flags;
         
   // TODO: Change these to malloc's later.
   esa->SA  = calloc( (n+2), sizeof(int) );
   esa->LCP = calloc( (n+2), sizeof(int) );
   
   
   printf("Constructing SA/LCP\n");
   // Construct the SA and LCP in linear time.
   sais((unsigned char*)s, esa->SA, esa->LCP, n);

   // This is needed for the child table values to be computed
   // correctly.
   esa->LCP[n] = 0;
   
   if (! (flags & NO_CHILD_TAB) )   
   {
      printf("Constructing Child table\n");
      // Child table, we attempt standard construction first,
      // then optimise it to occupy just one field.
      esa->up      = calloc( (n+2), sizeof(int) );
      esa->down    = calloc( (n+2), sizeof(int) );
      esa->accross = calloc( (n+2), sizeof(int) ); 

      // Create the child table.
      constructChildValues( esa );
   }
   
   // Construct the inverse suffix array.
   if (! (flags & NO_INV)  )
   {
      printf("Constructing SAi\n");
      esa->SAi = calloc( (n+2), sizeof(int) );
   
      for (int i=0; i<n; i++)
         esa->SAi[esa->SA[i]] = i;
   }
         
   // Initialise the RMQ structure.           
   if (! (flags & NO_RMQ) )
   {
      printf("initialising RMQ\n");
      RMQ_succinct(esa->LCP, n);  
   }
   
   printf("Done ESA\n");
}

/******************************************************************************/
// return 1 if a > b, 0 if a==b and -1 if a < b
// mle will store the location at which the match fails.
// In this way we can reduce the number of comparisons we have to do.

int str_gth(const char *a, const char *b, int n, int *i)
{
   for (; *i<n; (*i)++)
   {      
      if (a[*i]=='\0') return 0;
      if (b[*i]=='\0') return 0;
      
      if ((unsigned char)a[*i] > (unsigned char)b[*i]) return 1;
      if ((unsigned char)a[*i] < (unsigned char)b[*i]) return -1;
   }
   
   return 0;
}

/******************************************************************************/
// Find the first location of a the longest substring in O(n\log m) 
// [ O(n + \log m)  expected time ].
// This algorithm comes from Gusfield.

// l0 is the start position, r0 is the end position.
// l is the length found.
int findLongestSubstring(        const char*     p,                                
                                       int       m,   
                                       int*      l,                                                              
                                       int       l0, 
                                       int       r0, 
                                 const ESA*      esa                           )
{
   int min   = l0;
   int max   = r0;
   
   // The prefix lengths.
   int min_p = 0;
   int max_p = 0;
   
   // The number of comparisons done at any particular substring.
   int x;
   
   int longest_match = 0; 
   int longest_pos   = 0;
   
   int mid;

   do 
   {
      x = MIN(min_p, max_p);
      
      //printf("x: %d\n", x);
      
      mid   = min+(max-min)/2;      
      int c = str_gth(p, esa->t + esa->SA[mid], m, &x);    
      
      if (x>longest_match)
      {  
         longest_match = x;
         longest_pos   = esa->SA[mid];
      }
      
      //  printf("min, max: %d, %d, %d\n", min, max, mid);
      if (c == 1)
      {
         min   = mid+1;
         min_p = x;
      }  
      else if (c == -1)
      {
         max   = mid-1;
         max_p = x;
      }  
      else if (c == 0)
      {
         *l = longest_match;
         return esa->SA[mid];
       }  
   } while (min <= max);
   
   *l = longest_match;
   return longest_pos;
}

/******************************************************************************/
// This will find the first position of a given substring in the suffix
// array using the binary-search technique.
// We can then use the LCP array to find all given matches in time
// proportional to the number of occurences.

int findSubstringPosition(       const char*     p,           
                                       int       m,     
                                       int       l0, 
                                       int       r0, 
                                 const ESA*      esa                           )
{
   // This allows us to use a lookup table for the initial values. 
   int min   = l0;
   int max   = r0;
   
   // The prefix lengths.
   int min_p = 0;
   int max_p = 0;
      
   int mid;

   int x;

   do 
   {
      // This is the longest prefix of precomputed values.
      x = MIN(min_p, max_p);
      
      mid   = min+(max-min)/2;      
      int c = str_gth(p, esa->t + esa->SA[mid], m, &x);    
            
      
      //  printf("min, max: %d, %d, %d\n", min, max, mid);
      if (c == 1)
      {
         min   = mid+1;
         min_p = x;
      }  
   
      else if (c == -1)
      {
         max   = mid-1;
         max_p = x;
      }  
   
      // We found a match.
      else if (c == 0)
      {
         
         // Find the first instance of it in the suffix array.
         // This takes time proportional to the number of matches.
         while (mid > 0 && esa->LCP[mid] >= m) --mid;
         
         return mid;         
      }
   
   } while (min <= max);
   
   // Substring was not found anywhere.
   return -1;
}

/******************************************************************************/

void freeESA(ESA *esa)
{
   free(esa->SA);
   free(esa->LCP);
   free(esa->SAi);
}

/******************************************************************************/
