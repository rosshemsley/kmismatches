
/*******************************************************************************
* Macros:
*******************************************************************************/

//------------------------------------------------------------------------------ 
// Get the first child interval for an l-Interval.
//------------------------------------------------------------------------------
#define LI_FIRST_CHILD(i,j,esa) (i < esa->up[j+1] && esa->up[j+1] <= j)        \
                                 ? esa->up[j+1] : esa->down[i]

// Get the next child along.                                 
#define LI_NEXT_CHILD(i,esa) esa->accross[i]
                                                        
//------------------------------------------------------------------------------  
// Get the LCP value of an l-Interval.
//------------------------------------------------------------------------------ 
#define LI_GET_LCP(i,j,esa) (i < esa->up[j+1] && esa->up[j+1] <= j)            \
                             ? esa->LCP[esa->up[j+1]] : esa->LCP[esa->down[i]]
//------------------------------------------------------------------------------                                                       
#define MIN(X,Y) (X) < (Y) ? (X) : (Y)
//------------------------------------------------------------------------------ 
#define MAX(X,Y) (X) > (Y) ? (X) : (Y)              


#define ALPHABET_SIZE (1 << CHAR_BIT)               
/*******************************************************************************
* Types:
*******************************************************************************/

// Control the generation of the suffix array.

typedef enum { NO_CHILD_TAB=0x1, NO_INV=0x2, NO_RMQ=0x4 } ESA_FLAGS;

// This is a container for an Extended Suffix Array.
typedef struct _ESA
{
   int   n;
   const char *t;
   int  *SA;
   int  *SAi;
   int  *LCP;
   int  *up;
   int  *down;
   int  *accross;
   ESA_FLAGS flags;
} ESA;

// Container for a p-block in the p-representation
// TODO: it's no longer a triple..
typedef struct _pTriple
{
   int j;
   int l;
}  pTriple;


/*******************************************************************************
* Functions:
*******************************************************************************/
void displaySA(                  const ESA*      esa                          );
//------------------------------------------------------------------------------                                                                               
int LCE(                               int       i, 
                                       int       j, 
                                 const ESA*      esa                          );
//------------------------------------------------------------------------------
void constructESA(               const char*     s, 
                                       int       n, 
                                       ESA*      esa,   
                                       ESA_FLAGS flags                        );
//------------------------------------------------------------------------------
void construct_pRepresentation(        pTriple*  P,
                                 const char*     text, 
                                 const char*     pattern, 
                                 const ESA*      esa,
                                       int       n,
                                       int       m                            );                                       
//------------------------------------------------------------------------------   
int str_gth(                     const char*     a, 
                                 const char*     b, 
                                       int       n, 
                                       int*      i                            );
//------------------------------------------------------------------------------   
int findLongestSubstring(        const char*     p,                                
                                       int       m,   
                                       int*      l,                                                              
                                       int       l0, 
                                       int       r0, 
                                 const ESA*      esa                          );
//------------------------------------------------------------------------------                                                                           
int findLongestSubstring_simple(const char*     p,                                
                                       int       m,   
                                       int*      l,                                                              
                                       int       l0, 
                                       int       r0, 
                                 const ESA*      esa                          );                                 
//------------------------------------------------------------------------------                                          
int findSubstringPosition(       const char*     p,           
                                       int       m,    
                                       int       l0, 
                                       int       r0, 
                                 const ESA*      esa                          );
//------------------------------------------------------------------------------                                          
void freeESA(                          ESA*      esa                          );   
//------------------------------------------------------------------------------
int getInterval(int *i, int *j, int l, char c, const ESA *esa);                                
/******************************************************************************/
#ifdef TEST
int test_ESA();
#endif
/******************************************************************************/


    
