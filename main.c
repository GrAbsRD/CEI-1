#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdarg.h>
#include <math.h>

#include "DataType.h"
#include "UsualMath.h"
#include "FixedPoint.h"
#include "InfiniteSeries.h"
#include "AGM.h"
#include "OrthoPolynomial.h"
#include "IntegralSolverGC.h"
#include "IntegralSolverGL.h"

/******************************************************************************
 * Suppose you compile-link-run the code in the terminal ... $
 * Compile-Link: ...$ gcc *.c -o TestK -lm
 * Run:          ...$ ./TestK 10
 *               ...$ ./TestK 15 >> Values-K.txt
 *               ...$ ./TestK 30 
 * Note that for different points of numeric integral method, the
 * precision of the K will be different when k approaches to 1.
 *****************************************************************************/ 
int main(int argc, char* argv[])
{
    const int N = 50;
    double delta = 0.02;
    double k[N];
    double K1[N]; // Compute the K(k) with the AGM algorithm
    double K2[N]; // Compute the K(k) with the infinite series method
    double K3[N]; // Compute the K(k) with the Gauss-Legendre numeric integral method
    double K4[N]; // Compute the K(k) with the Gauss-Chebyshev numeric integral method

    double precision = 1e-9;
    unsigned npts_numeric_integral = atoi(argv[1]);

    printf(" K(k): Complete Elliptic Integral of the first kind\n");
    printf(" The value of K(k) with different computation methods\n");
    for(int i = 0; i < N; ++i)
    {
        k[i]  = i*delta;
        K1[i] = K_AGM(k[i], precision);
        K2[i] = K_Series(k[i], precision);
        K3[i] = K_IntGL(k[i], npts_numeric_integral);
        K4[i] = K_IntGC(k[i], npts_numeric_integral);
        printf(" k = %.2f, AGM K = %.7f, Series K = %.7f, IntGL K= %.7f, IntGC K = %.7f\n", 
                 k[i], K1[i], K2[i], K3[i], K4[i]);
    }

    return 0;    
}