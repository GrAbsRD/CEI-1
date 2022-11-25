#ifndef INFINITESERIES_H_INCLUDED
#define INFINITESERIES_H_INCLUDED

#include "DataType.h"
#include "UsualMath.h"

bool   IsNegligible(double term, double sum, double precision);

double InfiniteSeriesSolver(MapN2R CoefFun, double x, double precision);

double K_Series(double k, double precision);
double CalcCoefCm(unsigned m);

#endif  // INFINITESERIES_H_INCLUDED





