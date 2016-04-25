#include <stdlib.h>
#include <stdio.h>
#include <math.h>
//#include <string.h>

//#include "database.h"
#include "parameters.h"
#include "functions.h"

//#include "petscsys.h"

/*
 *	Purpose:
 *		Convert from rst to abc coordinates using a Duffy transform.
 *
 *	Comments:
 *		a, b, c must be filled with zeros upon entry to this function.
 *
 *	Notation:
 *
 *	References:
 */

void rst_to_abc(const unsigned int Nn, const unsigned int d, const double *rst, const double *a, const double *b,
                const double *c)
{
	unsigned int i;
	double r_i, s_i, t_i;
	double *r, *s, *t;

	if (d < 2)
		printf("Error: Conversion from rst_to_abc coordinates should only be used for d >= 2.\n"), exit(1);

	r = malloc(Nn * sizeof *r); // free
	s = malloc(Nn * sizeof *s); // free
	t = malloc(Nn * sizeof *t); // free

	for (i = 0; i < Nn; i++) {
		r[i] = rst[0*Nn+i];
		s[i] = rst[1*Nn+i];
		if (d > 2)
			t[i] = rst[2*Nn+i];
		else
			t[i] = -1.0/sqrt(6.0);
	}

	for (i = 0; i < Nn; i++) {
		r_i = r[i];
		s_i = s[i];
		t_i = t[i];

		if (array_norm_d(1,sqrt(6.0)*t_i-3.0,"Inf") > 100*EPS) {
			if (array_norm_d(1,2*sqrt(3.0)*s_i+sqrt(6.0)*t_i-3,"Inf") > 100*EPS)
				a[i] = -6.0*r_i/(2.0*sqrt(3.0)*s_i+sqrt(6.0)*t_i-3.0);
			b[i] = 1.0/3.0*(8.0*sqrt(3.0)*s_i/(3.0-sqrt(6.0)*t_i)-1.0);
		}
		c[i] = 0.5*(sqrt(6.0)*t_i-1.0);
	}

	free(r);
	free(s);
	free(t);
}
