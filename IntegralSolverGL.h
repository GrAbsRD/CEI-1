#include <math.h>
#include <stdlib.h>

#include "DataType.h"
#include "UsualMath.h"

double KerGL(double t, double k);

void   CalcVectorVfGL(
        double*  vfgl,
        MapR2R   f,
        unsigned n,
        double*  root,
        double   k,
        double   a,
        double   b);

double IntegralSolverGL(
        double   k, 
        MapR2R   f, 
        unsigned n,
        double*  root,
        double*  weight,
        double   a,
        double   b);

double K_IntGL(double k, double n);
