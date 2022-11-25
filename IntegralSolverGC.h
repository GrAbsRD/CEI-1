#ifndef INTEGRALSOLVERGC_H_INCLUDED
#define INTEGRALSOLVERGC_H_INCLUDED

#include "DataType.h"
#include "UsualMath.h"

double KerGC(double t, double k);

void   CalcVectorVfGC(
        double*  vfgc,
        MapR2R  f,
        unsigned n,
        double*  root,
        double   k,
        double   a);
double IntegralSolverGC(
        double   k, 
        MapR2R  f, 
        unsigned n,
        double*  root,
        double*  weight,
        double   a);
double K_IntGC(double k, double n);

#endif // INTEGRALSOLVERGC_H_INCLUDED