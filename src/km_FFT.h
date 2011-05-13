#include <fftw3.h>

void match_with_FFT(        int  *matches, 
                            char  symbol,
                      const char *text, 
                      const char *pattern, 
                            int   n,
                            int   m        );
                    
                            
// void multiply_half_complex(       double  *r, 
//                            const double  *a, 
 //                           const double  *b, 
  //                                int      N  );
          
int test_FFT_Matching();          
          

void match_with_FFT(        int  *matches, 
                            char  symbol,
                      const char *text, 
                      const char *pattern, 
                            int   n,
                            int   m        );

