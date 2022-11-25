#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

#include "DataType.h"
#include "UsualMath.h"
#include "FixedPoint.h"
#include "OrthoPolynomial.h"
#include "IntegralSolverGL.h"


void CalcVectorVfGL(
        double*  vfgl,
        MapR2R  f, 
        unsigned n, 
        double*  root, 
        double   k, 
        double   a, 
        double   b)
{
    double u = 0;
    for(int i = 0; i < n; ++i)
    {
        u = (b + a + (b-a)*root[i])/2;
        vfgl[i] = f(u, k)*(b-a)/2;
    }
}

double KerGL(double t, double k)
{
    double u = k*sin(t);
    return 1/sqrt(1-u*u);
}

double IntegralSolverGL(
        double   k, 
        MapR2R  f, 
        unsigned n,
        double*  root,
        double*  weight,
        double   a,
        double   b)
{
    double* vfgl = (double*)malloc(n*sizeof(double));
    CalcVectorVfGL(vfgl, f, n, root, k, a, b);
    double I = InnerProduct(weight, vfgl, n);
    free(vfgl);
    return I;
}                        

double K_IntGL(double k, double n)
{
    const double PI = acos(-1.0);
    MapR2R f = (MapR2R) KerGL;
    double  a = 0;
    double  b = PI/2;   
    double* root   = (double*) malloc(n*sizeof(double));
    double* weight = (double*) malloc(n*sizeof(double));
    CalcRootsLegendre(root, n);
    CalcWeightsGL(weight, root, n);
    double K = IntegralSolverGL(k, f, n, root, weight, a, b);
    free(root);
    free(weight);
    return K;
}