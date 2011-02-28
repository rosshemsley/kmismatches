
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <limits.h>

/******************************************************************************/

#define SWAP(x,y) int t;t=x;x=y;y=t;

/******************************************************************************/

void printUsage(){
	printf("About: Takes a file and generates a text and pattern.\n\n");
	printf("Usage: [inputfile] [outputfile] [m] [k] [options]\n\n");
	printf("Options: -n value -- specify a maximum length for the text.\n");
}

/******************************************************************************/

int loadData(       char ** data, 
              const char  *filename )
{
   FILE *fp;
   long n;

   if((fp = fopen(filename, "rb")) == NULL) {
      fprintf(stderr, "Cannot open file \"%s\".\n", filename);
      exit(0);
   }
  
  
  fseek(fp, 0, SEEK_END);
  n = ftell(fp);
  rewind(fp);
  
  *data = malloc(sizeof(char) * (n+1) );
  
  if(fread(*data, sizeof(char), n, fp) != (size_t)n)
  {
     fprintf(stderr, "Could not load file.\n");
     exit(0);
  }
  
  printf("%c", (*data)[n]);
  
  (*data)[n] = '\0';
  
  return n;
}

/******************************************************************************/


// main [infile] [outfile] [n] [m] [k] 
// Construct an input.
int main(int argc, char ** argv)
{



   unsigned long seed = time(NULL);
   // Seed the generator.
   srand( seed );

   int k = -1;
   int n =-1;
   int m =-1;
  
   char *infile  = NULL;
   char *outfile = NULL;
   
   int error=0;
   
   if (argc >= 5)
   {

      infile  = argv[1];
      outfile = argv[2];
      m       = atoi(argv[3]);
      k       = atoi(argv[4]);
      
      if (argc == 7)
      {

         if(strcmp(argv[5],"-n") == 0)
         n = atoi(argv[6]);
      }
   
   } else error=1;
   
   // Check for errors. 
   if (error>0)
   {
      printUsage();
      exit(1);
   }

   // Check bounds.

   if (! (k>0 && k<=m) )
   {
      printf("Invalid k value, '%d'\n", k);
      exit(1);
   } 
   
   if (n!=-1)
   {
      if (n<m)
      {
         printf("Invalid n value, '%d'\n", n);
         exit(1);
      }
      
   }
   



   // Load the input text we will be using.
   char *data;
   int t = loadData(&data, infile);
   

   
   if (n == -1)
      n=t;


  
   // The pattern will be of length m+1 with the terminating character.
   char *pattern     = malloc(sizeof(char) * (m+1));
   char *old_pattern = malloc(sizeof(char) * (m+1));
   
   // We choose a piece of the text to use as the pattern.
   int p = rand() % n-m+1;
   
   // Copy the the text from position p into the pattern.
   for (int i=0; i<m; i++)
      pattern[i] = data[i+p];
      
   pattern[m] = '\0';

   // Reference copy of the pattern.
   memcpy(old_pattern, pattern, sizeof(char)*(m+1));
   
   // We randomise the pattern by transpositions.
   // The first transposition creates at most two mismatches, all those 
   // after that create at most one more (we make sure no transposition
   // is the inverse of the previous transposition).
   int y;
   int x  = rand() % m;
      
   // The first character to be transposed.
   char c     = pattern[x]; 
   int  first = x;
   for (int i=0; i<k; i++)
   {
      // Choose other index for the transposition.
      while( (y = rand() % m) == x || pattern[x] == pattern[y] );
      

      // Do the transposition.
      SWAP(pattern[x], pattern[y]);
      
      x=y;
   }
   // Undo the first transposition.
   pattern[first] = c;

   int mismatches=0;
   // Count the mismatches.
   for (int i=0; i<m;i++)
   {
      if (pattern[i] != old_pattern[i])
         mismatches ++;
   }
   
   // Terminate the data at position n.
   data[n] = '\0';
   
   
   
      
   FILE *output = fopen(outfile, "w");
   
   fprintf(output, "%d %d %d %d\n", n, m, mismatches, p);
   fprintf(output, "%s\n", pattern);
   fprintf(output, "%s", data);
   fclose(output);
   
   
   printf("===============================================================\n" );
   printf("| Made test data.                                             |\n" );
   printf("|                                                             |\n" );
   printf("|    Loaded: %-10d bytes                                 |\n", t   );
   printf("|                                                             |\n" );
   printf("|       Text: %-10d bytes                                |\n", n   );
   printf("|    Pattern: %-10d bytes                                |\n", m   );
   printf("| Mismatches: %-10d                                      |\n", mismatches);
   printf("|   Position: %-10d                                      |\n", p   );
   printf("|       Seed: %-10ld                                      |\n",seed);
   printf("|                                                             |\n" );
   printf("===============================================================\n" );
   
   
   
   
   
   
   
   return 0;



}

/******************************************************************************/
