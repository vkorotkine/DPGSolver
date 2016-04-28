#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/*
 *	Purpose:
 *		Provide functions for testing basis_SI and grad_basis_SI.
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
