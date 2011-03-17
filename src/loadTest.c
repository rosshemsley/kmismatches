#include "loadTest.h"
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <string.h>
#include "loadTest.h"

/******************************************************************************/
// Load a test input.

void load(                       const char*     filename, 
                                       int*      n, 
                                       int*      m, 
                                       int*      k, 
                                       int*      pos, 
                                       char**    text, 
                                       char**    pattern                       ) 
{
   
   // Open the input file.
   FILE *f = fopen(filename, "r");

   if (f == NULL) {
      fprintf(stderr, "Failed to open file \"%s\"\n", filename);
      exit(1);
   }

   // The first line tells us how the rest of the file looks.  
   char  buff[256];   
   fgets(buff, 256, f);
   
   float R;
   long seed;
   
   if (sscanf(buff, "%d %d %d %d %f %ld", n, m, k, pos, &R, &seed) != 6)
   {
      fprintf(stderr, "File header improperly formatted.\n");
      exit(1);
   }
  
   // Allocate memory fot the text and pattern.
   // The text and pattern will be alloc'd to one contiguous block of memory.
   // This allows us to generate generalised suffix arrays.
   *text         = malloc( sizeof(char) * (*n+1) );
   *pattern      = malloc( sizeof(char) * (*m+1) ); 
  
   // Load the text and pattern into memory.
   fread(*pattern, sizeof(char), *m,   f);
   fgets( buff,    2,                  f);
   fread(*text,    sizeof(char), *n,   f);
   
   // Terminate: NOTE: We must add one char to the length to do this.
   (*text)[*n]    = '\0';
   (*pattern)[*m] = '\0';
/*   
   // We didn't load enough data for some reason.
   if (strlen(*text) != *n || strlen(*pattern) != *m)
   {
      fprintf(stderr, "There was a problem reading the input file.\n");
      exit(1);
   }
  */ 
   // Display the status.
   printf("===============================================================\n");
   printf("| Loaded test data.                                           |\n");
   printf("|                                                             |\n");
   printf("|       Text: %-10d bytes                                |\n", *n );
   printf("|    Pattern: %-10d bytes                                |\n", *m );
   printf("| Mismatches: %-10d                                      |\n", *k );
   printf("| Uniformity: %.3f                                           |\n", R   );      
   printf("|   Position: %-10d                                      |\n", *pos);
   printf("|                                                             |\n");
   printf("===============================================================\n");
  
   // Take account of the end of termination of the string.
   (*text)[*n]    = '\0';
   (*pattern)[*m] = '\0';
   *n = *n+1;
   *m = *m+1;
  
}

