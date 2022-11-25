#include <stdlib.h>
#include <stdio.h>

#include "DataType.h"
#include "UsualMath.h"
#include "FixedPoint.h"
#include "OrthoPolynomial.h"

double OrthoPolynom(
    double   x, 
    unsigned n, 
    double   c_0, 
    double   c_1, 
    double   c_2, 
    MapN2R A, 
    MapN2R B, 
    MapN2R C)
{
    double Y0 = c_0;
    double Y1 = c_1*x + c_2;
    double Yn = 0;
    if(n == 0) Yn = Y0;
    if(n == 1) Yn = Y1;
    for(int i = 2; i <= n; ++i)
    {
        Yn = (A(i)*x + B(i))*Y1 - C(i)*Y0;
        Y0 = Y1;
        Y1 = Yn;
    }
    return Yn;
}

double DerivOrthoPolynom(
    double x, 
    unsigned n, 
    double c_0, 
    double c_1, 
    double c_2, 
    MapN2R A, 
    MapN2R B, 
    MapN2R C)
{
    double Y0 = c_0;
    double Y1 = c_1*x + c_2;
    double DY0 = 0;
    double DY1 = c_1;
    double Yn  = 0;
    double DYn = 0;
    if(n == 0) { Yn = Y0; DYn = DY0; }
    if(n == 1) { Yn = Y1; DYn = DY1; }
    for(int i = 2; i <= n; ++i)
    {
        Yn  = (A(i)*x + B(i))*Y1 - C(i)*Y0;
        DYn = A(i)*Y1 + (A(i)*x + B(i))*DY1 - C(i)*DY0; 
        Y0  = Y1;
        Y1  = Yn;
        DY0 = DY1;
        DY1 = DYn;
    }
    return DYn;
}

double RootSolverOrthoPolynom(
    MapR2R   f, 
    MapR2R   Df, 
    double   x0, 
    unsigned n, 
    double precision)
{
    double guess   = x0;
	double improve = x0;
    do{
        guess   = improve;
        improve = NewtonUpdate(f, Df, guess, n);
    }while(!IsGoodEnough1d(guess, improve, precision));
    return improve;
}


void CalcRootsChebyshev(double* root, unsigned n)
{
	const double PI = acos(-1.0);
	for(int i = 0; i < n; ++i)
	{
        root[i] =  cos((2*i+1)*PI/(2*n));
	}
}

void CalcWeightsGC(double* weight, double* root, unsigned n)
{
    const double PI = acos(-1.0);
    for(int i = 0; i < n; ++i)
    {
        weight[i] = PI/n;
    }
}

/*****************************************************************
 *  Coefficients for computing the Legendre polynomial P_n(x) in 
 *  the iterative formula.
 *****************************************************************
 */
double AnLegendre(unsigned n) { return 2.0 -1.0/n; }
double BnLegendre(unsigned n) { return 0; }
double CnLegendre(unsigned n) { return 1.0 - 1.0/n; }

/*****************************************************************
 *  Legendre polynomial P_n(x) 
 *****************************************************************
 */
double PolynomLegendre(double x, unsigned n)
{
    double c_0 = 1;
    double c_1 = 1;
    double c_2 = 0;

    MapN2R A = AnLegendre;
    MapN2R B = BnLegendre;
    MapN2R C = CnLegendre;

    double Yn = OrthoPolynom(x, n, c_0, c_1, c_2, A, B, C);
    return Yn;
}

/*****************************************************************
 *  Derivative of the Legendre polynomial,  P'_n(x)
  *****************************************************************
 */

double DerivPolynomLegendre(double x, unsigned n)
{
    double c_0 = 1;
    double c_1 = 1;
    double c_2 = 0;

    MapN2R A = AnLegendre;
    MapN2R B = BnLegendre;
    MapN2R C = CnLegendre;

    double DYn = DerivOrthoPolynom(x, n, c_0, c_1, c_2, A, B, C);
    return DYn;
}

/*****************************************************************
  *  The n roots of Legendre polynomial P_n(x)
  *****************************************************************
 */
void CalcRootsLegendre(double* root, unsigned n)
{
	
	MapR2R   f  = (MapR2R  ) PolynomLegendre;
	MapR2R   Df = (MapR2R  ) DerivPolynomLegendre;
	double precision = 1e-9;
    
	double x0; // initial guess for the root;
	const double PI = acos(-1.0);
	for(int i = 0; i < n; ++i)
	{
		x0 = cos(PI*(i+0.75)/(n+0.25));
        root[i] =  RootSolverOrthoPolynom(f, Df, x0, n, precision);
	}
}

/*****************************************************************
 *  The n weights for the numeric integral method of Gauss-Legendre 
 *****************************************************************
 */
void CalcWeightsGL(double* weight, double* root, unsigned n)
{
    double x, u;
    for(int i = 0; i < n; ++i)
    {
        x = root[i];
        u = DerivPolynomLegendre(x, n);
        weight[i] = 2/((1-x*x)*(u*u));
    }
}