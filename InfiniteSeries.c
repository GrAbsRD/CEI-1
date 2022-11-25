#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

#include "DataType.h"
#include "UsualMath.h"
#include "InfiniteSeries.h"


bool IsNegligible(double term, double sum, double precision)
{
	const double SMALL_NUMBER = 1e-10;
	if (fabs(term)/(fabs(sum) + SMALL_NUMBER) < precision){ /* relative ratio */
		return true;  // the contribution of the term is negligible
	}else{
		return false; // the contribution of the term is not negligible
	}
}

double InfiniteSeriesSolver(MapN2R CoefFun, double x, double precision)
{
	double   c_m  = 0;
	double   term = 0;
	double   sum  = 0;  
	unsigned  m   = 0;
	do{
		c_m  = CoefFun(m);
		term = c_m * pow(x, m);	// generate the general term
		sum  += term;           // sum = sum + term
		++m;                    // m = m + 1;
	}while( !IsNegligible(term, sum, precision) );

	return sum;
}

/*************************************************************
   Compute the first complete elliptic integral K(k) with the
   Infinite series method
 *************************************************************/
double K_Series(double k, double precision)
{   
	const double PI = acos(-1.0); // PI = 3.1415926...
	MapN2R CoefFun = CalcCoefCm;
	double x = k*k;  // Attention, please!
	double sum = InfiniteSeriesSolver(CoefFun, x, precision);
 	double K = sum*PI/2;
    return K;
}

double CalcCoefCm(unsigned m)
{
	if ( m == 0) return 1;

	double prod = 1.0;
	for(int i = 1; i < m+1; ++i)
	{
		prod  *= 1.0 - 0.5/i; 
	}
	double c_m = prod*prod;
	return c_m;
}


