
/*******************************************************************************
* Macros:
*******************************************************************************/

//------------------------------------------------------------------------------ 
// Get the first child interval for an l-Interval.
//------------------------------------------------------------------------------
#define LI_FIRST_CHILD(i,j,esa) (i < esa->up[j+1] && esa->up[j+1] <= j)        \
                                 ? esa->up[j+1] : esa->down[i]
//------------------------------------------------------------------------------  
// Get the LCP value of an l-Interval.
//------------------------------------------------------------------------------ 
#define LI_GET_LCP(i,j,esa) (i < esa->up[j+1] && esa->up[j+1] <= j)            \
                             ? esa->LCP[esa->up[j+1]] : esa->LCP[esa->down[i]]
//------------------------------------------------------------------------------                                                       
#define MIN(X,Y) (X) < (Y) ? (X) : (Y)
//------------------------------------------------------------------------------ 
#define MAX(X,Y) (X) > (Y) ? (X) : (Y)                             
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


/*******************************************************************************
* Functions:
*******************************************************************************/
void displaySA(                  const ESA*      esa, 
                                 const char*     pattern, 
                                       int       m                            );
//------------------------------------------------------------------------------                                                                               
int LCE(                               int       i, 
                                       int       j, 
                                 const ESA*      esa                          );
//------------------------------------------------------------------------------                                                                           
void constructESA(               const char*     s, 
                                       int       n, 
                                       ESA*      esa                          );
//------------------------------------------------------------------------------   
int str_gth(                     const char*     a, 
                                 const char*     b, 
                                       int       n, 
                                       int*      i                            );
//------------------------------------------------------------------------------   
int findSubstring(                     int       l0, 
                                       int       r0, 
                                       int*      l, 
                                 const char*     p,
                                 const char*     t, 
                                 const ESA*      esa, 
                                       int       n                            );
//------------------------------------------------------------------------------                                          
void freeESA(                          ESA*      esa                          );                                   
/******************************************************************************/
    
