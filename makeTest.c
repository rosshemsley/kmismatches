
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <limits.h>

/******************************************************************************/

void printUsage(){
	printf("About: Takes a file and generates a text and pattern. ");
	printf("Usage: <input file> <outputfile> <m> <k> [options]\n\n");
	printf("Options: -n value -- specify a maximum length on the created text");
}

/******************************************************************************/

int loadData(       char ** data, 
              const char  *filename )
{
   FILE *fp;
   long n;

   if((fp = fopen(filename, "rb")) == NULL) {
      fprintf(stderr, "Cannot open file `%s': ", filename);
      exit(0);
   }
  
  
  fseek(fp, 0, SEEK_END);
  n = ftell(fp);
  rewind(fp);
  
  *data = malloc(sizeof(char) * n );
  
  if(fread(*data, sizeof(char), n, fp) != (size_t)n)
  {
     fprintf(stderr, "Could not load file.\n");
     exit(0);
  }
  
  return n;
}

/******************************************************************************/


// main [infile] [outfile] [n] [m] [k] 
// Construct an input.
int main(int argc, char ** argv)
{

   int k = -1;
   int n =-1;
   int m =-1;
  
   char *inFile  = NULL;
   char *outFile = NULL;
   
   int error=0;
   
   if (argc >= 5)
   {
      inFile  = argv[1];
      outFile = argv[2];
      m       = atoi(argv[3]);
      k       = atoi(argv[4]);
      
      if (argc==7)
      {
         if(strcmp(argv[5],"-n") == 0)
         n = atoi(argv[6]);
      } else error = 1; 
   
   } 
    
    
   if (error)
   {
      printUsage();
      exit(1);
   }

   
   




}

/******************************************************************************/
