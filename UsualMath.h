#ifndef USUALMATH_H_INCLUDED
#define USUALMATH_H_INCLUDED

#include "DataType.h" 

void   PrintRealArray(char* arrayName, double* a, unsigned size);
double InnerProduct(double* w, double* v, unsigned size);
double DistEuclid(double* u, double* v, unsigned size);
double NormEuclid(double* u, unsigned size);

double Norm1d(double p);   // special case of NormEuclid
double Norm2d(Point2d p);  // special case of NormEuclid
double Dist1d(double p, double q);     // special case of DistEuclid
double Dist2d(Point2d p, Point2d q);   // special case of DistEuclid
double RelativeError1d(double p, double q);
double RelativeError2d(Point2d p, Point2d q);

#endif // USUALMATH_H_INCLUDED