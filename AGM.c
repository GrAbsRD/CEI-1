#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <stdbool.h>

#include "DataType.h"
#include "AGM.h"
#include "FixedPoint.h"

Point2d AgmUpdate(Point2d p)
{
    Point2d q;
    q.x = (p.x + p.y)/2;
    q.y = sqrt(p.x*p.y);
    return q;
}

double AGM(double a, double b, double precision)
{
    Point2d guess = {a, b};
    MapP2P    f = (MapP2P) AgmUpdate;
    Point2d agm = FixedPointSolver2d(f,guess,precision);
    return  agm.x;
}

/*************************************************************
   Compute the first complete elliptic integral K(k) with the
   Arithmatic-Geometric Mean Algorithm 
 *************************************************************/
double K_AGM(double k, double precision)
{
    // arcsin(1.0) = PI/2 = 3.1415926.../2 
    double K = asin(1.0)/AGM(1.0, sqrt(1-k*k), precision);
    return K;  
}
