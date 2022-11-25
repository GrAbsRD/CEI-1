#ifndef ORTHOPOLYNOMIAL_H_INCLUDED
#define ORTHOPOLYNOMIAL_H_INCLUDED


//#include <stdlib.h>
//#include <stdio.h>

#include "DataType.h"
#include "UsualMath.h"

double OrthoPolynom(
    double x, 
    unsigned n, 
    double c_0, 
    double c_1, 
    double c_2, 
    MapN2R A,  
    MapN2R B, 
    MapN2R C);

double DerivOrthoPolynom(
    double x, 
    unsigned n, 
    double c_0, 
    double c_1, 
    double c_2, 
    MapN2R A,  
    MapN2R B, 
    MapN2R C);

double RootSolverOrthoPolynom(
    MapR2R f, 
    MapR2R Df, 
    double x0, 
    unsigned n, 
    double precision);


void  CalcRootsChebyshev(double* root, unsigned n);
void  CalcWeightsGC(double* weight, double* root, unsigned n);

double AnLegendre(unsigned n); 
double BnLegendre(unsigned n); 
double CnLegendre(unsigned n); 
double PolynomLegendre(double x, unsigned n);
double DerivPolynomLegendre(double x, unsigned n);
void   CalcRootsLegendre(double* root, unsigned n);
void   CalcWeightsGL(double* weight, double* root, unsigned n);

#endif // ORTHOPOLYNOMIAL_H_INCLUDED