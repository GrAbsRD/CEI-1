/* This program is used to find the 1-dim fixed-point specified by the
 * formula x = f(x, ...) where ... represents some extra parameters.
 */ 

#include <stdlib.h>
#include <stdio.h>

#include "DataType.h"
#include "UsualMath.h"
#include "FixedPoint.h"

bool IsGoodEnough1d(double guess, double improve, double precision)
{
	if( RelativeError1d(guess, improve) < precision )
		return true;
	else
		return false;
}

bool IsGoodEnough2d(Point2d guess, Point2d improve, double precision)
{
	if( RelativeError2d(guess, improve) < precision)
		return true;
	else
		return false;
}

Point2d FixedPointSolver2d(
	MapP2P    f, 
	Point2d   x0, 
	double    precision, 
	...)
{
	va_list ap;
	va_start(ap, precision);
	double para1 = va_arg(ap, double);
	double para2 = va_arg(ap, double);
	double para3 = va_arg(ap, double);
	va_end(ap);
	
	Point2d guess   = x0;
	Point2d improve = f(guess, para1, para2, para3);
	while( !IsGoodEnough2d(guess, improve, precision))
	{
		guess   = improve;
		improve = f(guess, para1, para2, para3);
	}
	return guess;	
}

void CheckSlope(double slope)
{   
    const double SMALL_NUMBER = 1e-10;
    if(fabs(slope) < SMALL_NUMBER)
    {
        printf("Exception! The slope/derivative is near zero!\n");
        exit(0);
    }
}

double NewtonUpdate(MapR2R f, MapR2R Df, double x, unsigned n)
{
	double slope = Df(x,n);
    CheckSlope(slope);
    double x_new =  x - f(x, n)/slope;
	return x_new;
}




