// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/// \file

#include "bases.h"

#include <math.h>
#include "gsl/gsl_sf_gamma.h"

#include "macros.h"
#include "definitions_tol.h"
#include "definitions_elements.h"

#include "matrix.h"

#include "math_functions.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

const struct const_Matrix_d* constructor_basis_tp_orthonormal
	(const int p_b, const struct const_Matrix_d*const rst)
{
	if (rst->layout != 'C')
		EXIT_UNSUPPORTED;

	const ptrdiff_t d   = rst->ext_1,
	                n_n = rst->ext_0,
	                pp1 = p_b+1,
	                n_b = compute_n_basis(d,p_b,ST_TP);

	struct Matrix_d* phi_rst = constructor_empty_Matrix_d('R',n_n,n_b); // returned
	double* phi_data = phi_rst->data;

	const double*const r = get_col_const_Matrix_d(0,rst),
	            *const s = ( d > 1 ? get_col_const_Matrix_d(1,rst) : NULL),
	            *const t = ( d > 2 ? get_col_const_Matrix_d(2,rst) : NULL);

	for (int n = 0; n < n_n; ++n) {
		for (int k = 0, k_max = GSL_MIN(GSL_MAX((d-2)*pp1,1),pp1); k < k_max; ++k) {
		for (int j = 0, j_max = GSL_MIN(GSL_MAX((d-1)*pp1,1),pp1); j < j_max; ++j) {
		for (int i = 0, i_max = GSL_MIN(GSL_MAX((d-0)*pp1,1),pp1); i < i_max; ++i) {
			           *phi_data  = jac_jacobi_normalized(r[n],i,0.0,0.0);
			if (d > 1) *phi_data *= jac_jacobi_normalized(s[n],j,0.0,0.0);
			if (d > 2) *phi_data *= jac_jacobi_normalized(t[n],k,0.0,0.0);
			++phi_data;
		}}}
	}

	return (const struct const_Matrix_d*) phi_rst;
}

const struct const_Matrix_d* constructor_basis_si_orthonormal
	(const int p_b, const struct const_Matrix_d*const rst)
{
	if (rst->layout != 'C')
		EXIT_UNSUPPORTED;

	const ptrdiff_t d   = rst->ext_1,
	                n_n = rst->ext_0,
	                n_b = compute_n_basis(d,p_b,ST_SI);

	if (d < 2 || d > 3)
		EXIT_UNSUPPORTED;

	const struct const_Matrix_d*const abc = constructor_abc_from_rst_si(rst); // destructed
	const double*const a = get_col_const_Matrix_d(0,abc),
	            *const b = get_col_const_Matrix_d(1,abc),
	            *const c = ( d > 2 ? get_col_const_Matrix_d(2,abc) : NULL);

	struct Matrix_d* phi_rst = constructor_empty_Matrix_d('R',n_n,n_b); // returned
	double* phi_data = phi_rst->data;

	for (int n = 0; n < n_n; ++n) {
		for (int i = 0, i_max = p_b;             i <= i_max; i++) {
		for (int j = 0, j_max = p_b-i;           j <= j_max; j++) {
		for (int k = 0, k_max = (d-2)*(p_b-i-j); k <= k_max; k++) {
			if (d == 2)
				*phi_data = 2.0/pow(3.0,0.25)*pow(1.0-b[n],i);
			else
				*phi_data = 4.0/pow(2.0,0.25)*pow(1.0-b[n],i)*pow(1.0-c[n],i+j);

			*phi_data *= jac_jacobi_normalized(a[n],i,0.0,0.0);
			*phi_data *= jac_jacobi_normalized(b[n],j,2.0*i+1.0,0.0);
			if (d == 3)
				*phi_data *= jac_jacobi_normalized(c[n],k,2.0*(i+j+1),0.0);
			++phi_data;
		}}}
	}
	destructor_const_Matrix_d(abc);

	return (const struct const_Matrix_d*) phi_rst;
}

ptrdiff_t compute_n_basis (const int d, const int p_b, const int super_type)
{
	if (super_type == ST_TP)
		return pow(p_b+1,d);
	else if (super_type == ST_SI)
		return round(gsl_sf_fact(d+p_b)/(gsl_sf_fact(d)*gsl_sf_fact(p_b)));
	else if (super_type == ST_PYR)
		EXIT_ADD_SUPPORT;
	else
		EXIT_UNSUPPORTED;
}

const struct const_Matrix_d* constructor_abc_from_rst_si (const struct const_Matrix_d*const rst)
{
	if (rst->layout != 'C')
		EXIT_UNSUPPORTED;

	const ptrdiff_t d   = rst->ext_1,
	                n_n = rst->ext_0;

	if (d < 2 || d > 3)
		EXIT_UNSUPPORTED;

	double* abc_d = malloc(n_n*d * sizeof *abc_d); // moved
	struct Matrix_d* abc = constructor_move_Matrix_d_d('C',n_n,d,true,abc_d); // returned

	const double*const r = get_col_const_Matrix_d(0,rst),
	            *const s = ( d > 1 ? get_col_const_Matrix_d(1,rst) : NULL),
	            *const t = ( d > 2 ? get_col_const_Matrix_d(2,rst) : NULL);

	double* a = get_col_Matrix_d(0,abc),
	      * b = get_col_Matrix_d(1,abc),
	      * c = ( d > 2 ? get_col_Matrix_d(2,abc) : NULL);

	for (int n = 0; n < n_n; ++n) {
		const double r_n = r[n],
		             s_n = s[n],
		             t_n = ( d == 2 ? -1.0/sqrt(6.0) : t[n] );

		if (fabs(2.0*sqrt(3.0)*s_n+sqrt(6.0)*t_n-3.0) > 1e2*EPS)
			a[n] = 6.0*r_n/(3.0-2.0*sqrt(3.0)*s_n-sqrt(6.0)*t_n);
		else // On top line of the regular TET / At the top of the regular TRI
			a[n] = 0.0;

		if (fabs(sqrt(6.0)*t_n-3.0) > 1e2*EPS)
			b[n] = 1.0/3.0*(8.0*sqrt(3.0)*s_n/(3.0-sqrt(6.0)*t_n)-1.0);
		else // At the top of the regular TET
			b[n] = 0.0;

		if (d == 3)
			c[n] = 0.5*(sqrt(6.0)*t_n-1.0);
	}
	return (const struct const_Matrix_d*) abc;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
