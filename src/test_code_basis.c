#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "parameters.h"

/*
 *	Purpose:
 *		Provide functions for testing basis_* and grad_basis_*.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

void get_scaling_basis_TRI(const unsigned int i, const unsigned int j, const double b, double *con_i, double *con_j,
                           double *con_b)
{
	*con_i = sqrt((2.0*i+1.0)/2.0);
	*con_j = sqrt((i+j+1.0)/pow(2.0,2.0*i+1.0));
	*con_b = pow(1.0-b,(double) i);
}

void get_scaling_basis_TET(const unsigned int i, const unsigned int j, const unsigned int k, const double b,
                           const double c, double *con_i, double *con_j, double *con_k, double *con_b, double *con_c)
{
	get_scaling_basis_TRI(i,j,b,con_i,con_j,con_b);
	*con_k = sqrt((2.0*(i+j+k)+3.0)/pow(2.0,2.0*(i+j)+3.0));
	*con_c = pow(1.0-c,(double) (i+j));
}

void get_scaling_basis_PYR(const unsigned int i, const unsigned int j, const unsigned int k, const double c,
                           double *con_i, double *con_j, double *con_k, double *con_c)
{
	double mu_ij = max(i,j);

	*con_i = sqrt((2.0*i+1.0)/2.0);
	*con_j = sqrt((2.0*j+1.0)/2.0);
	*con_k = sqrt((2.0*(mu_ij+k)+3.0)/pow(2.0,2.0*mu_ij+3.0));
	*con_c = pow(1.0-c,(double) mu_ij);
}

void poly2(const double *r, const double *s, const double *t, const unsigned int Nn, double **f, double **f_r,
           double **f_s, double **f_t)
{
	unsigned int i;
	double *f_rst, *f_r_rst, *f_s_rst, *f_t_rst;
	double r_i, s_i, t_i;

	f_rst   = malloc(Nn * sizeof *f_rst);   // keep (requires external free)
	f_r_rst = malloc(Nn * sizeof *f_r_rst); // keep (requires external free)
	f_s_rst = malloc(Nn * sizeof *f_s_rst); // keep (requires external free)
	f_t_rst = malloc(Nn * sizeof *f_t_rst); // keep (requires external free)

	for (i = 0; i < Nn; i++) {
		r_i = r[i];
		s_i = s[i];
		t_i = t[i];

		f_rst[i] = 1.0
		         + 2.0*r_i + 3.0*s_i + 4.0*t_i
		         + 5.0*r_i*r_i + 6.0*r_i*s_i + 7.0*r_i*t_i + 8.0*s_i*s_i + 9.0*s_i*t_i + 10.0*t_i*t_i;
		f_r_rst[i] = 2.0
		           + 10.0*r_i + 6.0*s_i + 7.0*t_i;
		f_s_rst[i] = 3.0
		           + 6.0*r_i + 16.0*s_i + 9.0*t_i;
		f_t_rst[i] = 4.0
		           + 7.0*r_i + 9.0*s_i + 20.0*t_i;
	}

	*f   = f_rst;
	*f_r = f_r_rst;
	*f_s = f_s_rst;
	*f_t = f_t_rst;
}

