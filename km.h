/******************************************************************************/

#ifndef km_ross_H_
#define km_ross_H_

/*******************************************************************************
* Constants:
*******************************************************************************/

// Maximum size of the alphabet. 
#define ALPHABET_SIZE (1 << CHAR_BIT)

/*******************************************************************************
* Types:
*******************************************************************************/

// This is a container for an Extended Suffix Array.
typedef struct _ESA
{
   int   n;
   int *SA;
   int *SAi;
   int *LCP;
   int *up;
   int *down;
   int *accross;
} ESA;

// Container for a p-block in the p-representation
// TODO: it's no longer a triple..
typedef struct _pTriple
{
   int j;
   int l;
}  pTriple;


/*******************************************************************************
* Function:
*******************************************************************************/

//----------------------------------------------------------------------------//   
// Load a correctly formatted test file.
// This will return a text, pattern and the details of the position where a
// match occurs.
//----------------------------------------------------------------------------//   
void load(                       const char*     filename, 
                                       int*      n, 
                                       int*      m, 
                                       int*      k, 
                                       int*      pos, 
                                       char**    text, 
                                       char**    pattern                      );    	
//----------------------------------------------------------------------------//   
// Perform k-mismatches using the method by Amir et al.   
//----------------------------------------------------------------------------//               
                     
void kmismatches(                const char*     text, 
                                 const char*     pattern,
                                       int       k,
                                       int       n,
                                       int       m,
                                       int*      matches                      );
//----------------------------------------------------------------------------//      
// Perform k-mismatches using the naive O(nk) method.
//----------------------------------------------------------------------------//                                 
void kmismatches_naive(          const char*     text, 
                                 const char*     pattern,
                                       int       k,
                                       int       n,
                                       int       m,
                                       int*      matches                      );
//----------------------------------------------------------------------------//
// Perform matching using Abrahamson/Kosaraju's method.
//----------------------------------------------------------------------------//   
void abrahamson_kosaraju(        const char*     text, 
                                 const char*     pattern,
                                       int       n,
                                       int       m,
                                       int*      matches                      );
//----------------------------------------------------------------------------//
// Perform matching using the Kangarooing method O(nk).
//----------------------------------------------------------------------------//   
void kangaroo(                   const char*     text, 
                                 const char*     pattern,
                                       int       k,
                                       int       n,
                                       int       m,
                                       int*      matches                      );
//----------------------------------------------------------------------------//
// Calculuate the hamming distance at every location using a naive O(nm) method.
//----------------------------------------------------------------------------//
void hamming_naive(             const char*      text, 
                                const char*      pattern,
                                      int        n,
                                      int        m,
                                      int*       matches                      );   
//----------------------------------------------------------------------------//                             
// Mark all the positoins where a given symbol matches using a lookup matrix.
//----------------------------------------------------------------------------//                                
void markMatches(               const int*       lookup,        
                                const int*       lookup_matrix,
                                      int        l,
                                      int*       matches,
                                const char*      text,
                                      int        n,       
                                      int        m                            );
/******************************************************************************/
#endif
/******************************************************************************/

