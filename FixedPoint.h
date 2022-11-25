#ifndef FIXEDPOINT_H_INCLUDED
#define FIXEDPOINT_H_INCLUDED

#include <stdarg.h>
#include <math.h>
#include <stdbool.h>

#include "DataType.h"
#include "UsualMath.h"


bool   IsGoodEnough1d(double guess,	double improve, double precision);
bool   IsGoodEnough2d(Point2d guess, Point2d improve, double precision);

Point2d FixedPointSolver2d(MapP2P f, Point2d guess, double precision, ...);

void   CheckSlope(double slope);
double NewtonUpdate(MapR2R f, MapR2R Df, double x, unsigned n);
double RootSolverOrthPolynom(MapR2R f, MapR2R Df, double x0, unsigned n, double precision);
void   CalcRootsLegendre(double* root, unsigned n);

#endif // FIXEDPOINT_H_INCLUDED
