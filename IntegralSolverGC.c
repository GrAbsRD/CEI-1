#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "DataType.h"
#include "UsualMath.h"
#include "OrthoPolynomial.h"
#include "IntegralSolverGC.h"



void CalcVectorVfGC(
        double*   vfgc, // vector v depends on f, a and root
        MapR2R   f, 
        unsigned  n, 
        double*   root, 
        double    k, 
        double    a)
{
    for(int i = 0; i < n; ++i)
    {
        vfgc[i] = f(a*root[i], k);
    }
}

double IntegralSolverGC(
            double k, 
            MapR2R f, 
            unsigned n, 
            double* root, 
            double* weight, 
            double a)
{
    double* vfgc = (double*) malloc(n*sizeof(double));
    CalcVectorVfGC(vfgc, f, n, root, k, a);
    double I = InnerProduct(weight, vfgc, n);
    free(vfgc);

    return I;
}

double KerGC(double t, double k)
{
    double u = k*t;
    return 0.5/sqrt(1-u*u);
}

double K_IntGC(double k, double n)
{
    double* root   = (double*) malloc(n*sizeof(double));
    double* weight = (double*) malloc(n*sizeof(double));
    CalcRootsChebyshev(root, n);
    CalcWeightsGC(weight, root, n);

    MapR2R f = (MapR2R) KerGC;
    double a = 1.0;
    double K = IntegralSolverGC(k, f, n, root, weight, a);
    free(root);
    free(weight);

    return K;
}