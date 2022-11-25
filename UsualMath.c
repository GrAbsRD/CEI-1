#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <stdbool.h>

#include "DataType.h"
#include "UsualMath.h"

double InnerProduct(double* w, double* v, unsigned size)
{
    double sum = 0;
    for(int i = 0; i < size; ++i)
    {
        sum += w[i] * v[i];
    }
    return sum;
}

void PrintRealArray(char* arrayName, double* a, unsigned size)
{
    printf("%s: ", arrayName);
    for(int i = 0; i < size; ++i)  printf(" %f", a[i]);
    printf("\n");
}

double NormEuclid(double* u, unsigned size)
{
    double sum = 0;
    for(int i = 0; i < size; ++i)
    {
        sum += u[i]*u[i];
    }
    double norm = sqrt(sum);
    return norm;
}

double DistEuclid(double* u, double* v, unsigned size)
{
	double sum = 0;
	for(int i = 0; i < size; ++i)
	{
		sum += (u[i] - v[i])*(u[i] - v[i]);
	}
	double dist = sqrt(sum);
	return dist;
}

double Norm1d(double x)
{
	return fabs(x);
}

double Norm2d(Point2d p)
{
	const int dim = 2;
	double v[] = {p.x, p.y};
    return NormEuclid(v,dim); 
}


double Dist1d(double guess, double improve)
{
	return Norm1d(guess - improve);
}

double Dist2d(Point2d p, Point2d q)
{
	const int dim = 2;
	double u[] = {p.x, p.y};
	double v[] = {q.x, q.y};
	return DistEuclid(u, v, dim);
}

double RelativeError1d(double p, double q)
{
	const double SMALL_NUMBER = 1e-10;
	// Attention, norm may be zero!
    double relative_error = Dist1d(p, q)/(Norm1d(p) + SMALL_NUMBER); 
	return relative_error;
}

double RelativeError2d(Point2d p, Point2d q)
{
	const double SMALL_NUMBER = 1e-10;
	// Attention, norm may be zero!
    double relative_error = Dist2d(p,q)/(Norm2d(p) + SMALL_NUMBER); 
	return relative_error;
}
