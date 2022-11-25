#ifndef AGM_H_INCLUDED
#define AGM_H_INCLUDED

//#include "FixedPoint.h"

Point2d AgmUpdate(Point2d p);

double AGM(double a, double b, double precision);

double K_AGM(double k, double precision);

#endif // AGM_H_INCLUDED