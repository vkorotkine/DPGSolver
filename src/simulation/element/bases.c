// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/// \file

#include "bases.h"

#include <math.h>
#include "gsl/gsl_math.h"
#include "gsl/gsl_sf_gamma.h"

#include "macros.h"
#include "definitions_tol.h"
#include "definitions_elements.h"

#include "multiarray.h"
#include "matrix.h"

#include "math_functions.h"

// Static function declarations ************************************************************************************* //

/** \brief Evaluate the Bernstein polynomial.
 *  \todo Delete this after the std interval basis is verified.
 *  \return See brief. */
static double bernstein
	(const int p,      ///< The order.
	 const int i,      ///< The index.
	 const double r_01 ///< The reference coordinate in [0,1].
	);

/** \brief Evaluate the derivative of the Bernstein polynomial.
 *  \todo Delete this after the std interval basis is verified.
 *  \return See brief. */
//static double grad_bernstein
double grad_bernstein
	(const int p,      ///< The order.
	 const int i,      ///< The index.
	 const double r_01 ///< The reference coordinate in [0,1].
	);

/** \brief Evaluate the Bernstein polynomial on the standard interval.
 *  \return See brief. */
static double bernstein_std
	(const int p,   ///< The order.
	 const int i,   ///< The index.
	 const double r ///< The reference coordinate in [-1,1].
	);

// Interface functions ********************************************************************************************** //
// Constructor functions ******************************************************************************************** //

const struct const_Matrix_d* constructor_basis_tp_orthonormal (const int p_b, const struct const_Matrix_d*const rst)
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
	}}}}

	return (const struct const_Matrix_d*) phi_rst;
}

const struct const_Multiarray_Matrix_d* constructor_grad_basis_tp_orthonormal
	(const int p_b, const struct const_Matrix_d*const rst)
{
	if (rst->layout != 'C')
		EXIT_UNSUPPORTED;

	const ptrdiff_t d   = rst->ext_1,
	                n_n = rst->ext_0,
	                pp1 = p_b+1,
	                n_b = compute_n_basis(d,p_b,ST_TP);

	struct Multiarray_Matrix_d* grad_phi_rst =
		constructor_empty_Multiarray_Matrix_d(false,1,(ptrdiff_t[]){d}); // returned

	double* grad_phi_data[d];
	for (int dim = 0; dim < d; ++dim) {
		grad_phi_rst->data[dim] = constructor_empty_Matrix_d('R',n_n,n_b); // keep
		grad_phi_data[dim] = grad_phi_rst->data[dim]->data;
	}

	const double*const r = get_col_const_Matrix_d(0,rst),
	            *const s = ( d > 1 ? get_col_const_Matrix_d(1,rst) : NULL),
	            *const t = ( d > 2 ? get_col_const_Matrix_d(2,rst) : NULL);

	if (d == 1) {
		for (int n = 0; n < n_n; ++n) {
		for (int i = 0; i < pp1; ++i) {
			*grad_phi_data[0]++  = jac_djacobi_normalized (r[n],i,0.0,0.0);
		}}
	} else if (d == 2) {
		for (int n = 0; n < n_n; ++n) {
		for (int j = 0; j < pp1; ++j) {
		for (int i = 0; i < pp1; ++i) {
			*grad_phi_data[0]++ = jac_djacobi_normalized(r[n],i,0.0,0.0)
			                     *jac_jacobi_normalized (s[n],j,0.0,0.0);
			*grad_phi_data[1]++ = jac_jacobi_normalized (r[n],i,0.0,0.0)
			                     *jac_djacobi_normalized(s[n],j,0.0,0.0);
		}}}
	} else if (d == 3) {
		for (int n = 0; n < n_n; ++n) {
		for (int k = 0; k < pp1; ++k) {
		for (int j = 0; j < pp1; ++j) {
		for (int i = 0; i < pp1; ++i) {
			*grad_phi_data[0]++ = jac_djacobi_normalized(r[n],i,0.0,0.0)
			                     *jac_jacobi_normalized (s[n],j,0.0,0.0)
			                     *jac_jacobi_normalized (t[n],k,0.0,0.0);
			*grad_phi_data[1]++ = jac_jacobi_normalized (r[n],i,0.0,0.0)
			                     *jac_djacobi_normalized(s[n],j,0.0,0.0)
			                     *jac_jacobi_normalized (t[n],k,0.0,0.0);
			*grad_phi_data[2]++ = jac_jacobi_normalized (r[n],i,0.0,0.0)
			                     *jac_jacobi_normalized (s[n],j,0.0,0.0)
			                     *jac_djacobi_normalized(t[n],k,0.0,0.0);
		}}}}
	}

	return (const struct const_Multiarray_Matrix_d*) grad_phi_rst;
}

const struct const_Matrix_d* constructor_basis_si_orthonormal (const int p_b, const struct const_Matrix_d*const rst)
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
	}}}}
	destructor_const_Matrix_d(abc);

	return (const struct const_Matrix_d*) phi_rst;
}

const struct const_Matrix_d* constructor_basis_pyr_orthonormal (const int p_b, const struct const_Matrix_d*const rst)
{
	if (rst->layout != 'C')
		EXIT_UNSUPPORTED;

	const ptrdiff_t d   = rst->ext_1,
	                n_n = rst->ext_0,
	                n_b = compute_n_basis(d,p_b,ST_PYR);

	if (d != 3)
		EXIT_UNSUPPORTED;

	const struct const_Matrix_d*const abc = constructor_abc_from_rst_pyr(rst); // destructed
	const double*const a = get_col_const_Matrix_d(0,abc),
	            *const b = get_col_const_Matrix_d(1,abc),
	            *const c = get_col_const_Matrix_d(2,abc);

	struct Matrix_d* phi_rst = constructor_empty_Matrix_d('R',n_n,n_b); // returned
	double* phi_data = phi_rst->data;

	for (int n = 0; n < n_n; ++n) {
	for (int i = 0, i_max = p_b; i <= i_max; i++) {
	for (int j = 0, j_max = p_b; j <= j_max; j++) {
		const int mu_ij = GSL_MAX(i,j);
		for (int k = 0, k_max = p_b-mu_ij; k <= k_max; k++) {
			*phi_data  = pow(2.0,1.25)*pow(1.0-c[n],mu_ij);
			*phi_data *= jac_jacobi_normalized(a[n],i,0.0,0.0);
			*phi_data *= jac_jacobi_normalized(b[n],j,0.0,0.0);
			*phi_data *= jac_jacobi_normalized(c[n],k,2.0*(mu_ij+1),0.0);

			++phi_data;
		}
	}}}
	destructor_const_Matrix_d(abc);

	return (const struct const_Matrix_d*) phi_rst;
}

const struct const_Matrix_d* constructor_basis_tp_bezier (const int p_b, const struct const_Matrix_d*const rst)
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
if (0) {
	for (int n = 0; n < n_n; ++n) {
	for (int k = 0, k_max = GSL_MIN(GSL_MAX((d-2)*pp1,1),pp1); k < k_max; ++k) {
		const double t_01 = ( d > 2 ? 0.5*(1.0+t[n]) : 0.5 );
		for (int j = 0, j_max = GSL_MIN(GSL_MAX((d-1)*pp1,1),pp1); j < j_max; ++j) {
			const double s_01 = ( d > 1 ? 0.5*(1.0+s[n]) : 0.5 );
			for (int i = 0, i_max = GSL_MIN(GSL_MAX((d-0)*pp1,1),pp1); i < i_max; ++i) {
				const double r_01 = 0.5*(1.0+r[n]);
				           *phi_data  = bernstein(p_b,i,r_01);
				if (d > 1) *phi_data *= bernstein(p_b,j,s_01);
				if (d > 2) *phi_data *= bernstein(p_b,k,t_01);
				++phi_data;
			}
		}
	}}
} else {
	for (int n = 0; n < n_n; ++n) {
	for (int k = 0, k_max = GSL_MIN(GSL_MAX((d-2)*pp1,1),pp1); k < k_max; ++k) {
	for (int j = 0, j_max = GSL_MIN(GSL_MAX((d-1)*pp1,1),pp1); j < j_max; ++j) {
	for (int i = 0, i_max = GSL_MIN(GSL_MAX((d-0)*pp1,1),pp1); i < i_max; ++i) {
		           *phi_data  = bernstein_std(p_b,i,r[n]);
		if (d > 1) *phi_data *= bernstein_std(p_b,j,s[n]);
		if (d > 2) *phi_data *= bernstein_std(p_b,k,t[n]);
		++phi_data;
	}}}}
}

	return (const struct const_Matrix_d*) phi_rst;
}

// Helper functions ************************************************************************************************* //

ptrdiff_t compute_n_basis (const int d, const int p_b, const int super_type)
{
	if (super_type == ST_TP)
		return pow(p_b+1,d);
	else if (super_type == ST_SI)
		return round(gsl_sf_fact(d+p_b)/(gsl_sf_fact(d)*gsl_sf_fact(p_b)));
	else if (super_type == ST_PYR)
		return round(1.0/6.0*((p_b+1)*(p_b+2)*(2*p_b+3)));
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
	            *const s = get_col_const_Matrix_d(1,rst),
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

const struct const_Matrix_d* constructor_abc_from_rst_pyr (const struct const_Matrix_d*const rst)
{
	if (rst->layout != 'C')
		EXIT_UNSUPPORTED;

	const ptrdiff_t d   = rst->ext_1,
	                n_n = rst->ext_0;

	if (d != 3)
		EXIT_UNSUPPORTED;

	double* abc_d = malloc(n_n*d * sizeof *abc_d); // moved
	struct Matrix_d* abc = constructor_move_Matrix_d_d('C',n_n,d,true,abc_d); // returned

	const double*const r = get_col_const_Matrix_d(0,rst),
	            *const s = get_col_const_Matrix_d(1,rst),
	            *const t = get_col_const_Matrix_d(2,rst);

	double* a = get_col_Matrix_d(0,abc),
	      * b = get_col_Matrix_d(1,abc),
	      * c = get_col_Matrix_d(2,abc);

	for (int n = 0; n < n_n; ++n) {
		const double r_n = r[n],
		             s_n = s[n],
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
	return (const struct const_Matrix_d*) abc;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static double bernstein (const int p, const int i, const double r_01)
{
	if ((i == -1) || (p-i == -1))
		return 0.0;

	return binomial_coef(p,i)*pow(r_01,i)*pow(1.0-r_01,p-i);
}

//static double grad_bernstein (const int p, const int i, const double r_01)
double grad_bernstein (const int p, const int i, const double r_01)
{
	return p*(bernstein(p-1,i-1,r_01) - bernstein(p-1,i,r_01));
}

static double bernstein_std (const int p, const int i, const double r)
{
	if ((i == -1) || (p-i == -1))
		return 0.0;

	return binomial_coef(p,i)*pow(0.5*(1.0-r),p-i)*pow(0.5*(1.0+r),i);
}
