/******************************************************************************/

#ifndef km_ross_H_
#define km_ross_H_

#include "esa.h"

/*******************************************************************************
* Macros
*******************************************************************************/

//----------------------------------------------------------------------------//   
// Maximum size of the alphabet. 
//----------------------------------------------------------------------------//   
#define ALPHABET_SIZE (1 << CHAR_BIT)
#define SWAP(x,y) int t;t=x;x=y;y=t;

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
//----------------------------------------------------------------------------//                             
// Verify a match in O(k) given a full generalised suffix array.
//----------------------------------------------------------------------------//
int verify(                           int        i, 
                                      int        j, 
                                      int        m, 
                                      int        k, 
                                const ESA*       esa                          );
//----------------------------------------------------------------------------//                             
// Verify a match in O(k) given a location in a p-representation
//----------------------------------------------------------------------------//                                
int verifyMatch(                 const pTriple*  pRepresentation,
                                 const char*     text,
                                 const char*     pattern,
                                 const ESA*      esa,                  
                                       // The position in the p-representation.              
                                       int       x,    
                                       int       t,    
                                       int       i,                            
                                       // problem-specific variables.
                                 const int       k,
                                 const int       n,
                                 const int       m                            );
//----------------------------------------------------------------------------//                             
// Verify a match in O(m).
//----------------------------------------------------------------------------//                                     
int verify_naive(                const char*     t, 
                                 const char*     p, 
                                       int       m, 
                                       int       k                            );
                            
/******************************************************************************/
#ifdef TEST                                    
int test_km();                                      
#endif
/******************************************************************************/
#endif
/******************************************************************************/

