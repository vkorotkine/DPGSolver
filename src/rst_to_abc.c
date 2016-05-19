// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "parameters.h"
#include "functions.h"

/*
 *	Purpose:
 *		Convert from rst to abc coordinates using Duffy-type transforms.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

void rst_to_abc_SI(const unsigned int Nn, const unsigned int d, const double *rst, double *a, double *b, double *c)
{
	unsigned int n;
	double r_n, s_n, t_n;
	double *r, *s, *t;

	if (d < 2)
		printf("Error: Conversion from rst to abc SI coordinates should only be used for d >= 2.\n"), exit(1);

	r = malloc(Nn * sizeof *r); // free
	s = malloc(Nn * sizeof *s); // free
	t = malloc(Nn * sizeof *t); // free

	for (n = 0; n < Nn; n++) {
		r[n] = rst[0*Nn+n];
		s[n] = rst[1*Nn+n];
		if (d > 2)
			t[n] = rst[2*Nn+n];
		else
			t[n] = -1.0/sqrt(6.0);
	}

	for (n = 0; n < Nn; n++) {
		r_n = r[n];
		s_n = s[n];
		t_n = t[n];

		if (fabs(2*sqrt(3.0)*s_n+sqrt(6.0)*t_n-3.0) > 100*EPS)
			a[n] = 6.0*r_n/(3.0-2.0*sqrt(3.0)*s_n-sqrt(6.0)*t_n);
		else // On top line of the regular TET / At the top of the regular TRI
			a[n] = 0.0;

		if (fabs(sqrt(6.0)*t_n-3.0) > 100*EPS)
			b[n] = 1.0/3.0*(8.0*sqrt(3.0)*s_n/(3.0-sqrt(6.0)*t_n)-1.0);
		else // At the top of the regular TET
			b[n] = 0.0;

		c[n] = 0.5*(sqrt(6.0)*t_n-1.0);
	}

	free(r);
	free(s);
	free(t);
}

void rst_to_abc_PYR(const unsigned int Nn, const unsigned int d, const double *rst, double *a, double *b, double *c)
{
	unsigned int n;
	double r_n, s_n, t_n;
	double *r, *s, *t;

	if (d != 3)
		printf("Error: Conversion from rst to abc PYR coordinates should only be used for d = 3.\n"), exit(1);

	r = malloc(Nn * sizeof *r); // free
	s = malloc(Nn * sizeof *s); // free
	t = malloc(Nn * sizeof *t); // free

	for (n = 0; n < Nn; n++) {
		r[n] = rst[0*Nn+n];
		s[n] = rst[1*Nn+n];
		t[n] = rst[2*Nn+n];
	}

	for (n = 0; n < Nn; n++) {
		r_n = r[n];
		s_n = s[n];
		t_n = t[n];

		if (fabs(0.8*sqrt(2.0)-t_n) > 100*EPS) {
			a[n] = r_n/(0.8-1.0/sqrt(2.0)*t_n);
			b[n] = s_n/(0.8-1.0/sqrt(2.0)*t_n);
		} else { // At the top of the pyramid
			a[n] = 0.0;
			b[n] = 0.0;
		}

		c[n] = sqrt(2.0)*t_n-0.6;
	}

	free(r);
	free(s);
	free(t);
}
