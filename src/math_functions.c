#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/*
 *	Purpose:
 *		Provide several standard math functions:
 *			int    factorial(const int n)  : Integer factorial function (n <= 32)
 *			int    gamma_i(const int n)    : Integer gamma function
 *			double gamma_d(const double x) : Double gamma function
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 *		Press(1992-2nd)_Numerical recipes in C- the art of scientific computing (Ch. 6.1)
 */

int factorial_i(const int n)
{
	static int ntop = 0,
	           a[33] = { 1 };
	int i;

	if (n < 0)
		printf("Error: Input to factorial_i must be greater than 0.\n"), exit(1);

	// Note: large values overflow
	if (n > 32)
		printf("Write a long long factorial function to avoid overflow.\n"), exit(1);

	// As ntop and a are static variables, multiple calls do not result in recomputation.
	while (ntop < n) {
		i = ntop++;
		a[ntop] = a[i]*ntop;
	}
	return a[n];
}

int gamma_i(const int n)
{
	if (n < 1)
		printf("Error: Input to gamma_i must be greater than 0.\n"), exit(1);

	return factorial_i(n-1);
}

double gamma_d(const double x)
{
	if (x <= 0.0)
		printf("Error: Input to gamma_d must be greater than 0.0.\n"), exit(1);

	if (floor(x) == x) {
		// Integer case
		return (double) gamma_i((int) x);
	} else {
		static double cof[6] = { 76.18009172947146,
		                        -86.50532032941677,
		                         24.01409824083091,
		                        -1.231739572450155,
		                         0.1208650973866179e-2,
		                        -0.5395239384953e-5     };
		int j;
		double z = x, y = x, tmp = x + 5.5,
		       ser = 1.000000000190015;
		tmp -= (z+0.5)*log(tmp);
		for (j = 0; j <=5; j++) {
			ser += cof[j]/++y;
		}

		return exp(-tmp+log(2.5066282746310005*ser/z));
	}

}
