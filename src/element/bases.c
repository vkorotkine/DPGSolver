/* {{{
This file is part of DPGSolver.

DPGSolver is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or any later version.

DPGSolver is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
Public License for more details.

You should have received a copy of the GNU General Public License along with DPGSolver.  If not, see
<http://www.gnu.org/licenses/>.
}}} */
/// \file

#include "bases.h"

#include <assert.h>
#include <math.h>
#include <string.h>
#include "gsl/gsl_math.h"
#include "gsl/gsl_sf_gamma.h"

#include "macros.h"
#include "definitions_bases.h"
#include "definitions_core.h"
#include "definitions_elements.h"
#include "definitions_math.h"
#include "definitions_tol.h"

#include "multiarray.h"
#include "matrix.h"

#include "math_functions.h"

// Static function declarations ************************************************************************************* //

/** \brief Evaluate the Bernstein polynomial on the standard interval.
 *  \return See brief. */
static double bernstein_std
	(const int p,   ///< The order.
	 const int i,   ///< The index.
	 const double r ///< The reference coordinate in [-1,1].
	);

/** \brief Evaluate the derivative of the Bernstein polynomial on the standard interval.
 *  \return See brief. */
static double grad_bernstein_std
	(const int p,   ///< The order.
	 const int i,   ///< The index.
	 const double r ///< The reference coordinate in [-1,1].
	);

// Interface functions ********************************************************************************************** //
// Constructor functions ******************************************************************************************** //
// Tensor-Product Orthonormal *************************************************************************************** //

const struct const_Matrix_d* constructor_basis_tp_orthonormal (const int p_b, const struct const_Matrix_d*const rst)
{
	assert(rst->layout == 'C');

	const int d   = (int)rst->ext_1,
	          pp1 = p_b+1;
	const ptrdiff_t n_n = rst->ext_0,
	                n_b = compute_n_basis(d,p_b,ST_TP);

	struct Matrix_d* phi_rst = constructor_empty_Matrix_d('R',n_n,n_b); // returned
	double* phi_data = phi_rst->data;

	if (d == 0) {
		assert(n_n == 1);
		assert(n_b == 1);
		*phi_data = 1.0;
	} else {
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
	}

	return (const struct const_Matrix_d*) phi_rst;
}

const struct const_Multiarray_Matrix_d* constructor_grad_basis_tp_orthonormal
	(const int p_b, const struct const_Matrix_d*const rst)
{
	assert(rst->layout == 'C');

	const int d   = (int)rst->ext_1,
	          pp1 = p_b+1;
	const ptrdiff_t n_n = rst->ext_0,
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
			*grad_phi_data[0]++ = jac_djacobi_normalized (r[n],i,0.0,0.0);
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

// Simplex Orthonormal ********************************************************************************************** //

const struct const_Matrix_d* constructor_basis_si_orthonormal (const int p_b, const struct const_Matrix_d*const rst)
{
	assert(rst->layout == 'C');

	const int d   = (int)rst->ext_1;
	const ptrdiff_t n_n = rst->ext_0,
	                n_b = compute_n_basis(d,p_b,ST_SI);

	assert(!(d < 2 || d > 3));

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

const struct const_Multiarray_Matrix_d* constructor_grad_basis_si_orthonormal
	(const int p_b, const struct const_Matrix_d*const rst)
{
	assert(rst->layout == 'C');

	const int d   = (int)rst->ext_1;
	const ptrdiff_t n_n = rst->ext_0,
	                n_b = compute_n_basis(d,p_b,ST_SI);

	assert(!(d < 2 || d > 3));

	const struct const_Matrix_d*const abc = constructor_abc_from_rst_si(rst); // destructed
	const double*const a = get_col_const_Matrix_d(0,abc),
	            *const b = get_col_const_Matrix_d(1,abc),
	            *const c = ( d > 2 ? get_col_const_Matrix_d(2,abc) : NULL);

	struct Multiarray_Matrix_d* grad_phi_rst =
		constructor_empty_Multiarray_Matrix_d(false,1,(ptrdiff_t[]){d}); // returned

	double* grad_phi_data[d];
	for (int dim = 0; dim < d; ++dim) {
		grad_phi_rst->data[dim] = constructor_empty_Matrix_d('R',n_n,n_b); // keep
		set_to_value_Matrix_d(grad_phi_rst->data[dim],0.0);

		grad_phi_data[dim] = grad_phi_rst->data[dim]->data;
	}

	const double con = ( d == 2 ? 2.0/pow(3.0,0.25) : 4.0/pow(2.0,0.25) );

	const int n_p_der = 2*d-1;
	double grad_phi_rst_part[d][n_p_der];
	for (int n = 0; n < n_n; n++) {
	for (int i = 0, i_max = p_b;             i <= i_max; i++) {
	for (int j = 0, j_max = p_b-i;           j <= j_max; j++) {
	for (int k = 0, k_max = (d-2)*(p_b-i-j); k <= k_max; k++) {
		const double a_n = a[n],
		             b_n = b[n],
		             c_n = ( d == 3 ? c[n] : -1.0 );

		for (int dim = 0; dim < d; dim++) {
		for (int pder = 0; pder < n_p_der; pder++ ) {
			grad_phi_rst_part[dim][pder] = 0.0;
		}}

		const double jPa  = jac_jacobi_normalized (a_n,i,0.0,          0.0),
		             jPb  = jac_jacobi_normalized (b_n,j,2.0*i+1.0,    0.0),
		             jPc  = jac_jacobi_normalized (c_n,k,2.0*(i+j+1.0),0.0),
		             djPa = jac_djacobi_normalized(a_n,i,0.0,          0.0),
		             djPb = jac_djacobi_normalized(b_n,j,2.0*i+1.0,    0.0),
		             djPc = jac_djacobi_normalized(c_n,k,2.0*(i+j+1.0),0.0);

		// Obtain contributions from each partial derivative
		grad_phi_rst_part[0][0] = djPa* jPb;
		grad_phi_rst_part[1][0] = djPa* jPb;
		grad_phi_rst_part[1][1] =  jPa*djPb;
		grad_phi_rst_part[1][2] =  jPa* jPb;

		if (d == 3) {
			// djPa = 0 when i = 0
			// djPb = 0 when j = 0
			// djPc = 0 when k = 0
			grad_phi_rst_part[0][0] *= jPc;
			grad_phi_rst_part[1][0] *= jPc;
			grad_phi_rst_part[1][1] *= jPc;
			grad_phi_rst_part[1][2] *= jPc;
			grad_phi_rst_part[2][0]  = djPa* jPb* jPc;
			grad_phi_rst_part[2][1]  =  jPa*djPb* jPc;
			grad_phi_rst_part[2][2]  =  jPa* jPb* jPc;
			grad_phi_rst_part[2][3]  =  jPa* jPb*djPc;
			grad_phi_rst_part[2][4]  =  jPa* jPb* jPc;
		}

		grad_phi_rst_part[0][0] *=  2.0                  ;
		grad_phi_rst_part[1][0] *=  2.0/3.0*SQRT3*a_n;
		grad_phi_rst_part[1][1] *=  2.0/3.0*SQRT3    ;
		grad_phi_rst_part[1][2] *= -2.0/3.0*SQRT3*i  ;
		if (i > 0) {
			grad_phi_rst_part[0][0] *= pow(1.0-b_n,i-1.0);
			grad_phi_rst_part[1][0] *= pow(1.0-b_n,i-1.0);
			grad_phi_rst_part[1][2] *= pow(1.0-b_n,i-1.0);
		} else {
			grad_phi_rst_part[0][0] *= 0.0; // redundant (i = 0 -> djPa = 0)
			grad_phi_rst_part[1][0] *= 0.0; // redundant (i = 0 -> djPa = 0)
			grad_phi_rst_part[1][2] *= 0.0;
		}
		grad_phi_rst_part[1][1] *= pow(1.0-b_n,i);

		if (d == 3) {
			grad_phi_rst_part[0][0] *=  2.0;
			grad_phi_rst_part[1][0] *=  2.0;
			grad_phi_rst_part[1][1] *=  2.0;
			grad_phi_rst_part[1][2] *=  2.0;
			grad_phi_rst_part[2][0] *=  2.0/3.0*SQRT6*a_n            ;
			grad_phi_rst_part[2][1] *=  1.0/6.0*SQRT6*(3.0*b_n+1.0)  ;
			grad_phi_rst_part[2][2] *= -1.0/6.0*SQRT6*(3.0*b_n+1.0)*i;
			grad_phi_rst_part[2][3] *=  1.0/2.0*SQRT6                ;
			grad_phi_rst_part[2][4] *= -1.0/2.0*SQRT6*(i+j)          ;
			if (i > 0) {
				grad_phi_rst_part[2][0] *= pow(1.0-b_n,i-1.0);
				grad_phi_rst_part[2][2] *= pow(1.0-b_n,i-1.0);
			} else {
				grad_phi_rst_part[2][0] *= 0.0;
				grad_phi_rst_part[2][2] *= 0.0;
			}
			grad_phi_rst_part[2][1] *= pow(1.0-b_n,i);
			grad_phi_rst_part[2][3] *= pow(1.0-b_n,i);
			grad_phi_rst_part[2][4] *= pow(1.0-b_n,i);

			if (i+j > 0) {
				grad_phi_rst_part[0][0] *= pow(1.0-c_n,i+j-1.0);
				grad_phi_rst_part[1][0] *= pow(1.0-c_n,i+j-1.0);
				grad_phi_rst_part[1][1] *= pow(1.0-c_n,i+j-1.0);
				grad_phi_rst_part[1][2] *= pow(1.0-c_n,i+j-1.0);
				grad_phi_rst_part[2][0] *= pow(1.0-c_n,i+j-1.0);
				grad_phi_rst_part[2][1] *= pow(1.0-c_n,i+j-1.0);
				grad_phi_rst_part[2][2] *= pow(1.0-c_n,i+j-1.0);
				grad_phi_rst_part[2][4] *= pow(1.0-c_n,i+j-1.0);
			} else {
				grad_phi_rst_part[0][0] *= 0.0;
				grad_phi_rst_part[1][0] *= 0.0;
				grad_phi_rst_part[1][1] *= 0.0;
				grad_phi_rst_part[1][2] *= 0.0;
				grad_phi_rst_part[2][0] *= 0.0;
				grad_phi_rst_part[2][1] *= 0.0;
				grad_phi_rst_part[2][2] *= 0.0;
				grad_phi_rst_part[2][4] *= 0.0;
			}
			grad_phi_rst_part[2][3] *= pow(1.0-c_n,i+j);
		}

		// Sum contributions from all partial derivatives
		for (int dim = 0; dim < d; dim++) {
		for (int pder = 0; pder < n_p_der; pder++) {
			*grad_phi_data[dim] += grad_phi_rst_part[dim][pder];
		}}

		// Add scaling constant and increment
		for (int dim = 0; dim < d; dim++) {
			*grad_phi_data[dim] *= con;
			grad_phi_data[dim]++;
		}
	}}}}
	destructor_const_Matrix_d(abc);

	return (const struct const_Multiarray_Matrix_d*) grad_phi_rst;
}

// Pyramid Orthonormal ********************************************************************************************** //

const struct const_Matrix_d* constructor_basis_pyr_orthonormal (const int p_b, const struct const_Matrix_d*const rst)
{
	assert(rst->layout == 'C');

	const int d = (int)rst->ext_1;
	const ptrdiff_t n_n = rst->ext_0,
	                n_b = compute_n_basis(d,p_b,ST_PYR);

	assert(d == 3);

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

const struct const_Multiarray_Matrix_d* constructor_grad_basis_pyr_orthonormal
	(const int p_b, const struct const_Matrix_d*const rst)
{
	assert(rst->layout == 'C');

	const int d = (int)rst->ext_1;
	const ptrdiff_t n_n = rst->ext_0,
	                n_b = compute_n_basis(d,p_b,ST_PYR);

	assert(d == 3);

	const struct const_Matrix_d*const abc = constructor_abc_from_rst_pyr(rst); // destructed
	const double*const a = get_col_const_Matrix_d(0,abc),
	            *const b = get_col_const_Matrix_d(1,abc),
	            *const c = get_col_const_Matrix_d(2,abc);

	struct Multiarray_Matrix_d* grad_phi_rst =
		constructor_empty_Multiarray_Matrix_d(false,1,(ptrdiff_t[]){d}); // returned

	double* grad_phi_data[d];
	for (int dim = 0; dim < d; ++dim) {
		grad_phi_rst->data[dim] = constructor_empty_Matrix_d('R',n_n,n_b); // keep
		set_to_value_Matrix_d(grad_phi_rst->data[dim],0.0);

		grad_phi_data[dim] = grad_phi_rst->data[dim]->data;
	}

	const double con = pow(2.0,1.25);

	const int n_p_der = 4;
	double grad_phi_rst_part[d][n_p_der];
	for (int n = 0; n < n_n; n++) {
	for (int i = 0, i_max = p_b; i <= i_max; i++) {
	for (int j = 0, j_max = p_b; j <= j_max; j++) {
		const int mu_ij = GSL_MAX(i,j);
		for (int k = 0, k_max = p_b-mu_ij; k <= k_max; k++) {
			const double a_n = a[n],
			             b_n = b[n],
			             c_n = c[n];

			for (int dim = 0; dim < d; dim++) {
			for (int pder = 0; pder < n_p_der; pder++ ) {
				grad_phi_rst_part[dim][pder] = 0.0;
			}}

			const double jPa  = jac_jacobi_normalized (a_n,i,0.0,            0.0),
			             jPb  = jac_jacobi_normalized (b_n,j,0.0,            0.0),
			             jPc  = jac_jacobi_normalized (c_n,k,2.0*(mu_ij+1.0),0.0),
			             djPa = jac_djacobi_normalized(a_n,i,0.0,            0.0),
			             djPb = jac_djacobi_normalized(b_n,j,0.0,            0.0),
			             djPc = jac_djacobi_normalized(c_n,k,2.0*(mu_ij+1.0),0.0);

			// Obtain contributions from each partial derivative
			grad_phi_rst_part[0][0] = djPa* jPb* jPc;
			grad_phi_rst_part[1][1] =  jPa*djPb* jPc;
			grad_phi_rst_part[2][0] = djPa* jPb* jPc;
			grad_phi_rst_part[2][1] =  jPa*djPb* jPc;
			grad_phi_rst_part[2][2] =  jPa* jPb*djPc;
			grad_phi_rst_part[2][3] =  jPa* jPb* jPc;

			grad_phi_rst_part[0][0] *=  2.0          ;
			grad_phi_rst_part[1][1] *=  2.0          ;
			grad_phi_rst_part[2][0] *=  sqrt(2.0)*a_n;
			grad_phi_rst_part[2][1] *=  sqrt(2.0)*b_n;
			grad_phi_rst_part[2][2] *=  sqrt(2.0)    ;
			grad_phi_rst_part[2][3] *= -sqrt(2.0)*mu_ij    ;

			if (mu_ij > 0) {
				grad_phi_rst_part[0][0] *= pow(1.0-c_n,mu_ij-1.0);
				grad_phi_rst_part[1][1] *= pow(1.0-c_n,mu_ij-1.0);
				grad_phi_rst_part[2][0] *= pow(1.0-c_n,mu_ij-1.0);
				grad_phi_rst_part[2][1] *= pow(1.0-c_n,mu_ij-1.0);
				grad_phi_rst_part[2][3] *= pow(1.0-c_n,mu_ij-1.0);
			} else {
				grad_phi_rst_part[0][0] *= 0.0;
				grad_phi_rst_part[1][1] *= 0.0;
				grad_phi_rst_part[2][0] *= 0.0;
				grad_phi_rst_part[2][1] *= 0.0;
				grad_phi_rst_part[2][3] *= 0.0;
			}
			grad_phi_rst_part[2][2] *= pow(1.0-c_n,(double) mu_ij);

			// Sum contributions from all partial derivatives
			for (int dim = 0; dim < d; dim++) {
			for (int pder = 0; pder < n_p_der; pder++) {
				*grad_phi_data[dim] += grad_phi_rst_part[dim][pder];
			}}

			// Add scaling constant
			for (int dim = 0; dim < d; dim++) {
				*grad_phi_data[dim] *= con;
				grad_phi_data[dim]++;
			}
		}
	}}}
	destructor_const_Matrix_d(abc);

	return (const struct const_Multiarray_Matrix_d*) grad_phi_rst;
}

// Tensor-Product Bezier ******************************************************************************************** //

const struct const_Matrix_d* constructor_basis_tp_bezier (const int p_b, const struct const_Matrix_d*const rst)
{
	assert(rst->layout == 'C');

	const int d   = (int)rst->ext_1,
	          pp1 = p_b+1;
	const ptrdiff_t n_n = rst->ext_0,
	                n_b = compute_n_basis(d,p_b,ST_TP);

	struct Matrix_d* phi_rst = constructor_empty_Matrix_d('R',n_n,n_b); // returned
	double* phi_data = phi_rst->data;

	const double*const r = get_col_const_Matrix_d(0,rst),
	            *const s = ( d > 1 ? get_col_const_Matrix_d(1,rst) : NULL),
	            *const t = ( d > 2 ? get_col_const_Matrix_d(2,rst) : NULL);

	if (d == 0) {
		assert(n_n == 1);
		assert(n_b == 1);
		*phi_data = 1.0;
	} else {
		for (int n = 0; n < n_n; ++n) {
		for (int k = 0, k_max = GSL_MIN(GSL_MAX((d-2)*pp1,1),pp1); k < k_max; ++k) {
		for (int j = 0, j_max = GSL_MIN(GSL_MAX((d-1)*pp1,1),pp1); j < j_max; ++j) {
		for (int i = 0; i < pp1; ++i) {
			           *phi_data  = bernstein_std(p_b,i,r[n]);
			if (d > 1) *phi_data *= bernstein_std(p_b,j,s[n]);
			if (d > 2) *phi_data *= bernstein_std(p_b,k,t[n]);
			++phi_data;
		}}}}
	}

	return (const struct const_Matrix_d*) phi_rst;
}

const struct const_Multiarray_Matrix_d* constructor_grad_basis_tp_bezier
	(const int p_b, const struct const_Matrix_d*const rst)
{
	assert(rst->layout == 'C');

	const int d   = (int)rst->ext_1,
	          pp1 = p_b+1;
	const ptrdiff_t n_n = rst->ext_0,
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
			*grad_phi_data[0]++ = grad_bernstein_std(p_b,i,r[n]);
		}}
	} else if (d == 2) {
		for (int n = 0; n < n_n; ++n) {
		for (int j = 0; j < pp1; ++j) {
		for (int i = 0; i < pp1; ++i) {
			*grad_phi_data[0]++ = grad_bernstein_std(p_b,i,r[n])
			                     *bernstein_std     (p_b,j,s[n]);
			*grad_phi_data[1]++ = bernstein_std     (p_b,i,r[n])
			                     *grad_bernstein_std(p_b,j,s[n]);
		}}}
	} else if (d == 3) {
		for (int n = 0; n < n_n; ++n) {
		for (int k = 0; k < pp1; ++k) {
		for (int j = 0; j < pp1; ++j) {
		for (int i = 0; i < pp1; ++i) {
			*grad_phi_data[0]++ = grad_bernstein_std(p_b,i,r[n])
			                     *bernstein_std     (p_b,j,s[n])
			                     *bernstein_std     (p_b,k,t[n]);
			*grad_phi_data[1]++ = bernstein_std     (p_b,i,r[n])
			                     *grad_bernstein_std(p_b,j,s[n])
			                     *bernstein_std     (p_b,k,t[n]);
			*grad_phi_data[2]++ = bernstein_std     (p_b,i,r[n])
			                     *bernstein_std     (p_b,j,s[n])
			                     *grad_bernstein_std(p_b,k,t[n]);
		}}}}
	}

	return (const struct const_Multiarray_Matrix_d*) grad_phi_rst;
}

// Simplex Bezier *************************************************************************************************** //

const struct const_Matrix_d* constructor_basis_si_bezier (const int p_b, const struct const_Matrix_d*const rst)
{
	assert(rst->layout == 'C');

	const int d   = (int)rst->ext_1;
	const ptrdiff_t n_n = rst->ext_0,
	                n_b = compute_n_basis(d,p_b,ST_SI);

	assert(!(d < 2 || d > 3));

	struct Matrix_d* phi_rst = constructor_empty_Matrix_d('C',n_n,n_b); // returned
	double* phi_data = phi_rst->data;

	const struct const_Matrix_d*const abc = constructor_abc_from_rst_si(rst); // destructed
	const double*const a = get_col_const_Matrix_d(0,abc),
	            *const b = get_col_const_Matrix_d(1,abc),
	            *const c = ( d > 2 ? get_col_const_Matrix_d(2,abc) : NULL);

	for (int i = 0, i_max = p_b;               i <= i_max; i++) {
	for (int j = 0, j_max = p_b-i;             j <= j_max; j++) {
	for (int k = 0, k_max = p_b-i-j;           k <= k_max; k++) {
	for (int l = 0, l_max = (d-2)*(p_b-i-j-k); l <= l_max; l++) {
		if (i+j+k+l != p_b)
			continue;

		for (int n = 0; n < n_n; ++n) {
			           *phi_data =  bernstein_std(i+j,    j,a[n])
			                       *bernstein_std(i+j+k,  k,b[n]);
			if (d > 2) *phi_data *= bernstein_std(i+j+k+l,l,c[n]);
			++phi_data;
		}
	}}}}
	destructor_const_Matrix_d(abc);
	transpose_Matrix_d(phi_rst,true);

	return (const struct const_Matrix_d*) phi_rst;
}

const struct const_Multiarray_Matrix_d* constructor_grad_basis_si_bezier
	(const int p_b, const struct const_Matrix_d*const rst)
{
	assert(rst->layout == 'C');

	const int d   = (int)rst->ext_1;
	const ptrdiff_t n_n = rst->ext_0,
	                n_b = compute_n_basis(d,p_b,ST_SI);

	assert(!(d < 2 || d > 3));

	struct Multiarray_Matrix_d* grad_phi_rst =
		constructor_empty_Multiarray_Matrix_d(false,1,(ptrdiff_t[]){d}); // returned

	double* grad_phi_data[d];
	for (int dim = 0; dim < d; ++dim) {
		grad_phi_rst->data[dim] = constructor_empty_Matrix_d('C',n_n,n_b); // keep
		grad_phi_data[dim] = grad_phi_rst->data[dim]->data;
	}


	const struct const_Matrix_d*const abc = constructor_abc_from_rst_si(rst); // destructed
	const double*const a = get_col_const_Matrix_d(0,abc),
	            *const b = get_col_const_Matrix_d(1,abc);

	/** \note Application of the chain-rule to the transformation from the reference to Duffy-type coordinates
	 *        results in singular terms which cancel with terms from the Bezier basis in a manner which is not
	 *        completely intuitive. Please consult [bezier_bases.pdf] for details of the derivation of the
	 *        gradient expressions below.
	 *
	 *  <!-- References: -->
	 *  [bezier_bases.pdf]: bases/bezier_bases.pdf
	 */
	if (d == 2) {
		for (int i = 0, i_max = p_b;     i <= i_max; i++) {
		for (int j = 0, j_max = p_b-i;   j <= j_max; j++) {
		for (int k = 0, k_max = p_b-i-j; k <= k_max; k++) {
			if (i+j+k != p_b)
				continue;

			for (int dim = 0; dim < d; ++dim) {
			for (int n = 0; n < n_n; ++n) {
/// \todo Return constants from a function.
				const double da_scale[] = { ( i+j > 0 ? (p_b)/((double)(i+j))*1.0            : 0.0),
				                            ( i+j > 0 ? (p_b)/((double)(i+j))*SQRT3/3.0*a[n] : 0.0), };
				const double db_scale[] = { 0.0,
				                            2.0*SQRT3/3.0, };
				*grad_phi_data[dim]++ =  grad_bernstein_std(i+j,    j,a[n])*da_scale[dim]
				                        *     bernstein_std(i+j+k-1,k,b[n])
				                      +       bernstein_std(i+j,    j,a[n])
				                        *grad_bernstein_std(i+j+k,  k,b[n])*db_scale[dim];
			}}
		}}}
	} else if (d == 3) {
		const double*const c = get_col_const_Matrix_d(2,abc);
		for (int i = 0, i_max = p_b;               i <= i_max; i++) {
		for (int j = 0, j_max = p_b-i;             j <= j_max; j++) {
		for (int k = 0, k_max = p_b-i-j;           k <= k_max; k++) {
		for (int l = 0, l_max = (d-2)*(p_b-i-j-k); l <= l_max; l++) {
			if (i+j+k+l != p_b)
				continue;

			for (int dim = 0; dim < d; ++dim) {
			for (int n = 0; n < n_n; ++n) {
				const double da_scale[] = { ( i+j > 0 ? p_b/((double)(i+j))*1.0            : 0.0),
				                            ( i+j > 0 ? p_b/((double)(i+j))*SQRT3/3.0*a[n] : 0.0),
				                            ( i+j > 0 ? p_b/((double)(i+j))*SQRT6/6.0*a[n] : 0.0), };
				const double db_scale[] = { 0.0,
				                            ( i+j+k > 0 ? (p_b)/((double)(i+j+k))*2.0*SQRT3/3.0             : 0.0),
				                            ( i+j+k > 0 ? (p_b)/((double)(i+j+k))*SQRT6/12.0*(3.0*b[n]+1.0) : 0.0), };
				const double dc_scale[] = { 0.0,
				                            0.0,
				                            SQRT6/2.0, };
				*grad_phi_data[dim]++ =  grad_bernstein_std(i+j,      j,a[n])*da_scale[dim]
				                        *     bernstein_std(i+j+k-1,  k,b[n])
				                        *     bernstein_std(i+j+k+l-1,l,c[n])
				                      +       bernstein_std(i+j,      j,a[n])
				                        *grad_bernstein_std(i+j+k,    k,b[n])*db_scale[dim]
				                        *     bernstein_std(i+j+k+l-1,l,c[n])
				                      +       bernstein_std(i+j,      j,a[n])
				                        *     bernstein_std(i+j+k,    k,b[n])
				                        *grad_bernstein_std(i+j+k+l,  l,c[n])*dc_scale[dim];
			}}
		}}}}
	}
	destructor_const_Matrix_d(abc);

	for (int dim = 0; dim < d; ++dim)
		transpose_Matrix_d(grad_phi_rst->data[dim],true);

	return (const struct const_Multiarray_Matrix_d*) grad_phi_rst;
}

// Pyramid Bezier *************************************************************************************************** //

const struct const_Matrix_d* constructor_basis_pyr_bezier (const int p_b, const struct const_Matrix_d*const rst)
{
	assert(rst->layout == 'C');

	const int d = (int)rst->ext_1;
	const ptrdiff_t n_n = rst->ext_0,
	                n_b = compute_n_basis(d,p_b,ST_PYR);

	assert(d == 3);

	struct Matrix_d* phi_rst = constructor_empty_Matrix_d('C',n_n,n_b); // returned
	double* phi_data = phi_rst->data;

	const struct const_Matrix_d*const abc = constructor_abc_from_rst_pyr(rst); // destructed
	const double*const a = get_col_const_Matrix_d(0,abc),
	            *const b = get_col_const_Matrix_d(1,abc),
	            *const c = get_col_const_Matrix_d(2,abc);

	for (int i = 0, i_max = p_b; i <= i_max; i++) {
	for (int j = 0, j_max = p_b; j <= j_max; j++) {
		const int mu_ij = GSL_MAX(i,j);
		for (int k = 0, k_max = p_b-mu_ij; k <= k_max; k++) {
		for (int n = 0; n < n_n; ++n) {
			*phi_data  = bernstein_std(p_b-k,i,a[n])
			            *bernstein_std(p_b-k,j,b[n])
			            *bernstein_std(p_b,  k,c[n]);
			++phi_data;
		}}
	}}
	destructor_const_Matrix_d(abc);
	transpose_Matrix_d(phi_rst,true);

	return (const struct const_Matrix_d*) phi_rst;
}

const struct const_Multiarray_Matrix_d* constructor_grad_basis_pyr_bezier
	(const int p_b, const struct const_Matrix_d*const rst)
{
	assert(rst->layout == 'C');

	const int d = (int)rst->ext_1;
	const ptrdiff_t n_n = rst->ext_0,
	                n_b = compute_n_basis(d,p_b,ST_PYR);

	assert(d == 3);

	struct Multiarray_Matrix_d* grad_phi_rst =
		constructor_empty_Multiarray_Matrix_d(false,1,(ptrdiff_t[]){d}); // returned

	double* grad_phi_data[d];
	for (int dim = 0; dim < d; ++dim) {
		grad_phi_rst->data[dim] = constructor_empty_Matrix_d('C',n_n,n_b); // keep
		grad_phi_data[dim] = grad_phi_rst->data[dim]->data;
	}

	const struct const_Matrix_d*const abc = constructor_abc_from_rst_pyr(rst); // destructed
	const double*const a = get_col_const_Matrix_d(0,abc),
	            *const b = get_col_const_Matrix_d(1,abc),
	            *const c = get_col_const_Matrix_d(2,abc);

	/// See \ref constructor_grad_basis_si_bezier for relevant comments.
	for (int i = 0, i_max = p_b; i <= i_max; i++) {
	for (int j = 0, j_max = p_b; j <= j_max; j++) {
		const int mu_ij = GSL_MAX(i,j);
		for (int k = 0, k_max = p_b-mu_ij; k <= k_max; k++) {
		for (int dim = 0; dim < d; ++dim) {
		for (int n = 0; n < n_n; n++) {
			const double da_scale[] = { ( p_b-k > 0 ? p_b/((double)(p_b-k))*1.0            : 0.0),
			                            0.0,
			                            ( p_b-k > 0 ? p_b/((double)(p_b-k))*SQRT2/2.0*a[n] : 0.0), };
			const double db_scale[] = { 0.0,
			                            ( p_b-k > 0 ? p_b/((double)(p_b-k))*1.0            : 0.0),
			                            ( p_b-k > 0 ? p_b/((double)(p_b-k))*SQRT2/2.0*b[n] : 0.0), };
			const double dc_scale[] = { 0.0,
			                            0.0,
			                            SQRT2, };
			*grad_phi_data[dim]++ =  grad_bernstein_std(p_b-k,i,a[n])*da_scale[dim]
			                        *     bernstein_std(p_b-k,j,b[n])
			                        *     bernstein_std(p_b-1,k,c[n])
			                      +       bernstein_std(p_b-k,i,a[n])
			                        *grad_bernstein_std(p_b-k,j,b[n])*db_scale[dim]
			                        *     bernstein_std(p_b-1,k,c[n])
			                      +       bernstein_std(p_b-k,i,a[n])
			                        *     bernstein_std(p_b-k,j,b[n])
			                        *grad_bernstein_std(p_b,  k,c[n])*dc_scale[dim];
		}}}
	}}
	destructor_const_Matrix_d(abc);

	for (int dim = 0; dim < d; ++dim)
		transpose_Matrix_d(grad_phi_rst->data[dim],true);

	return (const struct const_Multiarray_Matrix_d*) grad_phi_rst;
}

// Helper functions ************************************************************************************************* //

ptrdiff_t compute_n_basis (const int d, const int p_b, const int super_type)
{
	if (super_type == ST_TP)
		return (ptrdiff_t)pow(p_b+1,d);
	else if (super_type == ST_SI)
		return (ptrdiff_t)round(gsl_sf_fact((unsigned)(d+p_b))/(gsl_sf_fact((unsigned)d)*gsl_sf_fact((unsigned)p_b)));
	else if (super_type == ST_PYR)
		return (ptrdiff_t)round(1.0/6.0*((p_b+1)*(p_b+2)*(2*p_b+3)));
	else
		EXIT_UNSUPPORTED;
}

const struct const_Matrix_d* constructor_abc_from_rst_si (const struct const_Matrix_d*const rst)
{
	assert(rst->layout == 'C');

	const ptrdiff_t d   = rst->ext_1,
	                n_n = rst->ext_0;

	assert(!(d < 2 || d > 3));

	double* abc_d = malloc((size_t)(n_n*d) * sizeof *abc_d); // moved
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
		             t_n = ( d == 2 ? -1.0/SQRT6 : t[n] );

		if (fabs(2.0*SQRT3*s_n+SQRT6*t_n-3.0) > 1e2*EPS)
			a[n] = 6.0*r_n/(3.0-2.0*SQRT3*s_n-SQRT6*t_n);
		else // On top line of the regular TET / At the top of the regular TRI
			a[n] = 0.0;

		if (fabs(SQRT6*t_n-3.0) > 1e2*EPS)
			b[n] = 1.0/3.0*(8.0*SQRT3*s_n/(3.0-SQRT6*t_n)-1.0);
		else // At the top of the regular TET
			b[n] = 0.0;

		if (d == 3)
			c[n] = 0.5*(SQRT6*t_n-1.0);
	}
	return (const struct const_Matrix_d*) abc;
}

const struct const_Matrix_d* constructor_bcoords_from_rst_si (const struct const_Matrix_d*const rst)
{
	assert(rst->layout == 'C');

	const ptrdiff_t d   = rst->ext_1,
	                n_n = rst->ext_0;

	assert(!(d < 2 || d > 3));

	double* bcoords_d = malloc((size_t)(n_n*(d+1)) * sizeof *bcoords_d); // moved
	struct Matrix_d* bcoords = constructor_move_Matrix_d_d('C',n_n,d,true,bcoords_d); // returned

	const double*const r = get_col_const_Matrix_d(0,rst),
	            *const s = get_col_const_Matrix_d(1,rst),
	            *const t = ( d > 2 ? get_col_const_Matrix_d(2,rst) : NULL);

	double* u = get_col_Matrix_d(0,bcoords),
	      * v = get_col_Matrix_d(1,bcoords),
	      * w = get_col_Matrix_d(2,bcoords),
	      * x = ( d > 2 ? get_col_Matrix_d(3,bcoords) : NULL);

	if (d == 2) {
		for (int n = 0; n < n_n; ++n) {
			const double r_n = r[n],
			             s_n = s[n];

			u[n] = 1.0/3.0 - 1.0/2.0*r_n - 1.0/(2.0*SQRT3)*s_n;
			v[n] = 1.0/3.0 + 1.0/2.0*r_n - 1.0/(2.0*SQRT3)*s_n;
			w[n] = 1.0/3.0               + 1.0/(SQRT3)    *s_n;
		}
	} else if (d == 3) {
		for (int n = 0; n < n_n; ++n) {
			const double r_n = r[n],
			             s_n = s[n],
			             t_n = t[n];

			u[n] = 1.0/4.0 - 1.0/2.0*r_n - 1.0/(2.0*SQRT3)*s_n - 1.0/(2.0*SQRT6)*t_n;
			v[n] = 1.0/4.0 + 1.0/2.0*r_n - 1.0/(2.0*SQRT3)*s_n - 1.0/(2.0*SQRT6)*t_n;
			w[n] = 1.0/4.0               + 1.0/(SQRT3)    *s_n - 1.0/(2.0*SQRT6)*t_n;
			x[n] = 1.0/4.0                                         + 3.0/(2.0*SQRT6)*t_n;
		}
	}
	return (const struct const_Matrix_d*) bcoords;
}

const struct const_Matrix_d* constructor_abc_from_rst_pyr (const struct const_Matrix_d*const rst)
{
	assert(rst->layout == 'C');

	const ptrdiff_t d   = rst->ext_1,
	                n_n = rst->ext_0;

	assert(d == 3);

	double* abc_d = malloc((size_t)(n_n*d) * sizeof *abc_d); // moved
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

constructor_basis_fptr get_constructor_basis_by_super_type (const int s_type, const char*const ref_basis_name)
{
	if (strcmp(ref_basis_name,"orthonormal") == 0) {
		switch (s_type) {
			case ST_TP:  return constructor_basis_tp_orthonormal;  break;
			case ST_SI:  return constructor_basis_si_orthonormal;  break;
			case ST_PYR: return constructor_basis_pyr_orthonormal; break;
			default:     EXIT_ERROR("Unsupported: %d\n",s_type);   break;
		}
	} else if (strcmp(ref_basis_name,"bezier") == 0) {
		return get_constructor_basis_bezier_by_super_type(s_type);
	}
	EXIT_ERROR("Did not find the basis with the specified inputs: (%d, %s)\n",s_type,ref_basis_name);
}

constructor_basis_fptr get_constructor_basis_bezier_by_super_type (const int s_type)
{
	switch (s_type) {
		case ST_TP:  return constructor_basis_tp_bezier;  break;
		case ST_SI:  return constructor_basis_si_bezier;  break;
		case ST_PYR: return constructor_basis_pyr_bezier; break;
		default:     EXIT_ERROR("Unsupported: %d\n",s_type); break;
	}
}

constructor_basis_fptr get_constructor_basis_by_super_type_i (const int s_type, const int ind_basis)
{
	switch (ind_basis) {
	case BASIS_ORTHO:
		return get_constructor_basis_by_super_type(s_type,"orthonormal");
		break;
	case BASIS_BEZIER:
		return get_constructor_basis_by_super_type(s_type,"bezier");
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",ind_basis);
		break;
	}
//	EXIT_ERROR("Did not find the basis with the specified inputs: (%d, %s)\n",s_type,ref_basis_name);
}

constructor_grad_basis_fptr get_constructor_grad_basis_by_super_type (const int s_type, const char*const ref_basis_name)
{
	if (s_type == ST_TP) {
		if (strcmp(ref_basis_name,"orthonormal") == 0)
			return  constructor_grad_basis_tp_orthonormal;
		else if (strcmp(ref_basis_name,"bezier") == 0)
			return  constructor_grad_basis_tp_bezier;
	} else if (s_type == ST_SI) {
		return  constructor_grad_basis_si_orthonormal;
	} else if (s_type == ST_PYR) {
		return  constructor_grad_basis_pyr_orthonormal;
	}

	EXIT_ERROR("Did not find the basis with the specified inputs: (%d, %s)\n",s_type,ref_basis_name);
}

int get_basis_i_from_s (const char*const basis_name_s)
{
	int basis_name_i = -1;
	if (strcmp(basis_name_s,"orthonormal") == 0)
		basis_name_i = BASIS_ORTHO;
	else if (strcmp(basis_name_s,"lagrange") == 0)
		basis_name_i = BASIS_LAGRANGE;
	else if (strcmp(basis_name_s,"bezier") == 0)
		basis_name_i = BASIS_BEZIER;
	else
		EXIT_ERROR("Unsupported: %s\n",basis_name_s);

	return basis_name_i;
}

double B_Spline_Basis_ip(int i, int p, double xi, const struct const_Multiarray_d* knots){

	/*
	NOTE: Special treatment is made for xi values that are at the edge of the knot
		domain, in order to ensure that the value lies within the domain. This
		is done because otherwise 0/0 errors would occur in the NURBS basis function
		computation. Look into this further to perhaps find a more elegant solution.
	*/

	const double*const knots_i = get_col_const_Multiarray_d(0, knots);

	// Handle the case where xi is close to the knot domain edges. Ensure that 
	// xi is in the domain (by adding or subtracting EPS)

	if (fabs(xi - knots_i[0]) < EPS)
		xi = knots_i[0] + EPS;  // left boundary
	else if(fabs(xi-knots_i[knots->extents[0]-1]) < EPS)
		xi = knots_i[knots->extents[0]-1] - EPS;  // right boundary
		
	
	if (p == 0){

		// ====================
		// 		Base Case
		// ====================

		if (xi < knots_i[i+1] && xi >= knots_i[i])
			return 1.0;
		else
			return 0.0;
		
	} else {

		// ====================
		// 	  Recursive Case
		// ====================

		double xi_i 			= knots_i[i],
			   xi_iPlus1 		= knots_i[i+1],
			   xi_iPlusP 		= knots_i[i+p],
			   xi_iPlusPPlus1	= knots_i[i+p+1];

		// The first term
		double num1 	= (xi - xi_i) * B_Spline_Basis_ip(i, p-1, xi, knots), 
			   denom1 	= xi_iPlusP - xi_i;

		double term1;
		if (fabs(num1) < EPS && fabs(denom1) < EPS)
			term1 = 0.0;
		else
			term1 = num1/denom1;

		// The second term
		double num2 	= (xi_iPlusPPlus1 - xi) * B_Spline_Basis_ip(i+1, p-1, xi, knots),
			   denom2 	= xi_iPlusPPlus1 - xi_iPlus1;

		double term2;
		if (fabs(num2) < EPS && fabs(denom2) < EPS)
			term2 = 0.0;
		else
			term2 = num2/denom2;
		
		double return_value = term1 + term2;

		if (isnan(return_value)){

			printf("\n********* ERROR **********\n");
			printf("i = %d \n", i);
			printf("p = %d \n", p);
			printf("xi = %e \n", xi);
			printf("knots: \n");
			print_const_Multiarray_d(knots);

			EXIT_ERROR("NAN in B_Spline_Basis_ip\n");
		} else if(isinf(return_value)){

			printf("\n********* ERROR **********\n");
			printf("i = %d \n", i);
			printf("p = %d \n", p);
			printf("xi = %e \n", xi);
			printf("knots: \n");
			print_const_Multiarray_d(knots);

			EXIT_ERROR("INF in B_Spline_Basis_ip\n");
		}

		return return_value;
	}
}

const struct const_Multiarray_d *B_Spline_Basis_p(int p, 
	const struct const_Multiarray_d* xi_vals, const struct const_Multiarray_d* knots){

	const double*const xi_vals_i = get_col_const_Multiarray_d(0, xi_vals);

	// Compute the number of basis functions (is a function of the length
	// of the knot vector and order)
	int num_basis, num_xi_vals;
	num_basis = (int)knots->extents[0] - p - 1;
	num_xi_vals = (int)xi_vals->extents[0];

	struct Multiarray_d* basis_values = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){num_basis,num_xi_vals}); // returned

	// Loop through all the basis functions and xi values and evaluate the 
	// basis functions.
	int i,j;
	double basis_ip_value;

	for (j = 0; j < num_xi_vals; j++){
		for (i = 0; i < num_basis; i++){

			basis_ip_value = B_Spline_Basis_ip(i, p, xi_vals_i[j], knots);
			get_col_Multiarray_d(j, basis_values)[i] = basis_ip_value;   

		}
	}

	return (struct const_Multiarray_d*) basis_values;
}


double derivative_B_Spline_Basis_ip(int i, int p, double xi, const struct const_Multiarray_d* knots){

	const double*const knots_i = get_col_const_Multiarray_d(0, knots);

	// First Term:
	double num1, denom1, first_term;

	num1 = p*B_Spline_Basis_ip(i, p-1, xi, knots);
	denom1 = knots_i[i+p] - knots_i[i];

	if (fabs(num1) < EPS && fabs(denom1) < EPS){
		first_term = 0.0;
	} else{
		first_term = num1/denom1;
	}

	// Second Term:
	double num2, denom2, second_term;

	num2 = p*B_Spline_Basis_ip(i+1, p-1, xi, knots);
	denom2 = knots_i[i+p+1] - knots_i[i+1];

	if (fabs(num2) < EPS && fabs(denom2) < EPS){
		second_term = 0.0;
	} else{
		second_term = num2/denom2;
	}

	double return_value = first_term - second_term;

	if (isnan(return_value)){

		printf("\n********* ERROR **********\n");
		printf("i = %d \n", i);
		printf("p = %d \n", p);
		printf("xi = %e \n", xi);
		printf("knots: \n");
		print_const_Multiarray_d(knots);

		EXIT_ERROR("NAN in derivative_B_Spline_Basis_ip\n");
	} else if(isinf(return_value)){

		printf("\n********* ERROR **********\n");
		printf("i = %d \n", i);
		printf("p = %d \n", p);
		printf("xi = %e \n", xi);
		printf("knots: \n");
		print_const_Multiarray_d(knots);
		
		EXIT_ERROR("INF in derivative_B_Spline_Basis_ip\n");
	}

	return return_value;
}


const struct const_Multiarray_d *derivative_B_Spline_Basis_p(int p, 
	const struct const_Multiarray_d* xi_vals, const struct const_Multiarray_d* knots){


	const double*const xi_vals_i = get_col_const_Multiarray_d(0, xi_vals);

	// Compute the number of basis functions (is a function of the length
	// of the knot vector and order)
	int num_basis, num_xi_vals;
	num_basis = (int)knots->extents[0] - p - 1;
	num_xi_vals = (int)xi_vals->extents[0];

	struct Multiarray_d* del_basis_values = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){num_basis,num_xi_vals}); // returned

	// Loop through all the basis functions and xi values and evaluate the 
	// basis functions. Store the results in the multiarray
	int i,j;
	double del_basis_ip_value;

	for (j = 0; j < num_xi_vals; j++){
		for (i = 0; i < num_basis; i++){

			del_basis_ip_value = derivative_B_Spline_Basis_ip(i, p, xi_vals_i[j], knots);
			get_col_Multiarray_d(j, del_basis_values)[i] = del_basis_ip_value;   

		}
	}

	return (struct const_Multiarray_d*) del_basis_values;
}


double NURBS_Basis_ip(int i, const struct const_Multiarray_d* B_Spline_Basis_values,
 const struct const_Multiarray_d* weights) {

	const double*const B_Spline_Basis_values_i = get_col_const_Multiarray_d(0, B_Spline_Basis_values);
	const double*const weights_i = get_col_const_Multiarray_d(0, weights);

	// The numerator and denominator for the rational basis function
	double num, denom; 

	// Numerator:
	num = weights_i[i] * B_Spline_Basis_values_i[i];

	// Denominator: This is the "weight function" evaluated at xi
	denom = 0.0;
	int j;
	for (j = 0; j < B_Spline_Basis_values->extents[0]; j++){
		denom += weights_i[j]*B_Spline_Basis_values_i[j];
	}

	return num / denom;
}

const struct const_Multiarray_d *NURBS_Basis_p(const struct const_Multiarray_d* B_Spline_Basis_values, 
	const struct const_Multiarray_d* weights){

	int num_xi_vals = (int)B_Spline_Basis_values->extents[1];

	//Number of NURBS basis functions is equivalent to the number of B spline ones
	int num_basis = (int)B_Spline_Basis_values->extents[0];  

	// Construct the multiarray to hold the values to be returned
	struct Multiarray_d* NURBS_basis_values = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){num_basis,num_xi_vals}); // returned

	int i,j;

	for(j = 0; j < num_xi_vals; j++){
		// Loop over the number of xi values

		// Get the column of B Spline basis function values evaluated at the given 
		// xi value and move the data into a Multiarray
		const double *const B_Spline_basis_values_xi_val_i = get_col_const_Multiarray_d(j, B_Spline_Basis_values);
		
		const struct const_Multiarray_d* B_Spline_basis_values_xi_val = constructor_move_const_Multiarray_d_d('C',
			2, (ptrdiff_t[]){B_Spline_Basis_values->extents[0],1}, false, B_Spline_basis_values_xi_val_i);

		for(i = 0; i < num_basis; i++){
			// Loop over all the basis functions

			get_col_Multiarray_d(j, NURBS_basis_values)[i] = NURBS_Basis_ip(i, B_Spline_basis_values_xi_val, weights);

		}

		destructor_const_Multiarray_d(B_Spline_basis_values_xi_val);
	}

	return (struct const_Multiarray_d*) NURBS_basis_values;
}


double derivative_NURBS_Basis_ip(int i, const struct const_Multiarray_d* B_Spline_Basis_values,
	const struct const_Multiarray_d* derivative_B_Spline_Basis_values, const struct const_Multiarray_d* weights) {

	const double*const B_Spline_Basis_values_i = get_col_const_Multiarray_d(0, B_Spline_Basis_values);
	const double*const derivative_B_Spline_Basis_values_i = get_col_const_Multiarray_d(0, derivative_B_Spline_Basis_values);
	const double*const weights_i = get_col_const_Multiarray_d(0, weights);

	double 	w_i, // The ith basis function weight
			weight_func,  // The weight function value at the given xi point
			del_weight_func;  // The weight function derivative's value at the given xi point

	// Compute the weight function terms
	weight_func = 0;
	del_weight_func = 0;

	int j;
	for (j = 0; j < B_Spline_Basis_values->extents[0]; j++){
		
		weight_func += weights_i[j]*B_Spline_Basis_values_i[j];
		del_weight_func += weights_i[j]*derivative_B_Spline_Basis_values_i[j];

	}

	w_i = weights_i[i];

	return w_i*((weight_func*derivative_B_Spline_Basis_values_i[i] - del_weight_func*B_Spline_Basis_values_i[i])/(pow(weight_func, 2.0)));

}


const struct const_Multiarray_d *derivative_NURBS_Basis_p(const struct const_Multiarray_d* B_Spline_Basis_values, 
	const struct const_Multiarray_d* derivative_B_Spline_Basis_values, const struct const_Multiarray_d* weights){

	int num_xi_vals = (int)B_Spline_Basis_values->extents[1];

	//Number of NURBS basis functions is equivalent to the number of B spline ones
	int num_basis = (int)B_Spline_Basis_values->extents[0];  

	// Construct the multiarray to hold the values to be returned
	struct Multiarray_d* del_NURBS_basis_values = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){num_basis,num_xi_vals}); // returned

	int i,j;

	for(j = 0; j < num_xi_vals; j++){
		// Loop over the number of xi values

		// Get the column of B Spline basis function values evaluated at the given 
		// xi value and move the data into a Multiarray. Do the same with the derivative data
		const double *const B_Spline_basis_values_xi_val_i = get_col_const_Multiarray_d(j, B_Spline_Basis_values);
		const double *const derivative_B_Spline_basis_values_xi_val_i = get_col_const_Multiarray_d(j, derivative_B_Spline_Basis_values);
		
		const struct const_Multiarray_d* B_Spline_basis_values_xi_val = constructor_move_const_Multiarray_d_d('C',
			2, (ptrdiff_t[]){B_Spline_Basis_values->extents[0],1}, false, B_Spline_basis_values_xi_val_i);
		const struct const_Multiarray_d* derivative_B_Spline_basis_values_xi_val = constructor_move_const_Multiarray_d_d('C',
			2, (ptrdiff_t[]){B_Spline_Basis_values->extents[0],1}, false, derivative_B_Spline_basis_values_xi_val_i);
		
		for(i = 0; i < num_basis; i++){
			// Loop over all the basis functions
			get_col_Multiarray_d(j, del_NURBS_basis_values)[i] = derivative_NURBS_Basis_ip(i, B_Spline_basis_values_xi_val, 
					derivative_B_Spline_basis_values_xi_val, weights);

		}

		destructor_const_Multiarray_d(B_Spline_basis_values_xi_val);
		destructor_const_Multiarray_d(derivative_B_Spline_basis_values_xi_val);

	}

	return (struct const_Multiarray_d*) del_NURBS_basis_values;
}


const struct const_Multiarray_d *NURBS_basis_pq(const struct const_Multiarray_d* N_p_xi, 
	const struct const_Multiarray_d* N_q_eta, const struct const_Multiarray_d* weights){

	// The denominator of the rational basis function
	int num_i_basis, num_j_basis, num_pts;

	// The number of basis functions in each direction
	num_i_basis = (int)N_p_xi->extents[0];
	num_j_basis = (int)N_q_eta->extents[0];

	num_pts = (int)N_p_xi->extents[1];
	assert(num_pts == (int)N_q_eta->extents[1]);

	// Place the basis values into a multiarray to be returned. 
	struct Multiarray_d *basis_vals = constructor_empty_Multiarray_d('C',3,(ptrdiff_t[]){num_i_basis, num_j_basis, num_pts});  // returned

	double *basis_vals_i = basis_vals->data;
	double basis_val, weight_function;

	for (int k = 0; k < num_pts; k++){
		// Loop over the points

		const double *const N_p_xi_vals  = get_col_const_Multiarray_d(k, N_p_xi);
		const double *const N_q_eta_vals = get_col_const_Multiarray_d(k, N_q_eta);

		weight_function = 0;
		// Compute the weight function at the given point
		for (int i = 0; i < num_i_basis; i++){
			for (int j = 0; j < num_j_basis; j++){
				weight_function += N_p_xi_vals[i]*N_q_eta_vals[j]*get_col_const_Multiarray_d(j, weights)[i];
			}
		}

		// Compute the NURBS basis function values
		for (int i = 0; i < num_i_basis; i++){
			for (int j = 0; j < num_j_basis; j++){

				basis_val = (get_col_const_Multiarray_d(j, weights)[i]*N_p_xi_vals[i]*N_q_eta_vals[j])/weight_function;
				basis_vals_i[i + j*num_i_basis + k*num_i_basis*num_j_basis] = basis_val;

			}
		}
	}

	return (const struct const_Multiarray_d*)basis_vals;

}


const struct const_Multiarray_d *grad_NURBS_basis_pq(const struct const_Multiarray_d* N_p_xi,
	const struct const_Multiarray_d* dN_p_xi, const struct const_Multiarray_d* N_q_eta,
	const struct const_Multiarray_d* dN_q_eta,const struct const_Multiarray_d* weights){

	int num_i_basis, num_j_basis, num_pts;

	// The number of basis functions in each direction
	num_i_basis = (int)N_p_xi->extents[0];
	num_j_basis = (int)N_q_eta->extents[0];

	assert(num_i_basis == (int)dN_p_xi->extents[0]);
	assert(num_j_basis == (int)dN_q_eta->extents[0]);

	// The number of xi,eta points to consider
	num_pts = (int)N_p_xi->extents[1];

	assert(num_pts == (int)N_q_eta->extents[1]);
	assert(num_pts == (int)dN_p_xi->extents[1]);
	assert(num_pts == (int)dN_q_eta->extents[1]);


	// Place the partial values into a multiarray to be returned. 
	struct Multiarray_d *basis_grad_vals = constructor_empty_Multiarray_d('C',4,
		(ptrdiff_t[]){num_i_basis, num_j_basis, num_pts, 2});  // returned

	double *basis_grad_vals_i = basis_grad_vals->data;

	// Compute the weight function and its gradient
	double 	w_func,	w_func_del_xi, w_func_del_eta;
	double w_ij;
	double R_ij_pq_del_xi, R_ij_pq_del_eta;

	for (int k = 0; k < num_pts; k++){
		// Loop over all the points

		const double *const N_p_xi_vals  = get_col_const_Multiarray_d(k, N_p_xi);
		const double *const dN_p_xi_vals = get_col_const_Multiarray_d(k, dN_p_xi);

		const double *const N_q_eta_vals  = get_col_const_Multiarray_d(k, N_q_eta);
		const double *const dN_q_eta_vals = get_col_const_Multiarray_d(k, dN_q_eta);

		// Get the weight function values at the kth (xi,eta) point
		w_func 			= 0.0;
		w_func_del_xi 	= 0.0;
		w_func_del_eta 	= 0.0;

		for (int i = 0; i < num_i_basis; i++){
			for (int j = 0; j < num_j_basis; j++){

				// The weight value for the given i,j index
				w_ij = get_col_const_Multiarray_d(j, weights)[i];

				w_func 			+=  N_p_xi_vals[i] *  N_q_eta_vals[j]  * w_ij;
				w_func_del_xi 	+= dN_p_xi_vals[i] *  N_q_eta_vals[j]  * w_ij;
				w_func_del_eta 	+=  N_p_xi_vals[i] * dN_q_eta_vals[j]  * w_ij;

			}
		}

		// Compute the gradient values for each i,j basis function at the given kth point
		for (int i = 0; i < num_i_basis; i++){
			for (int j = 0; j < num_j_basis; j++){

				w_ij = get_col_const_Multiarray_d(j, weights)[i];

				R_ij_pq_del_xi = (w_ij*N_q_eta_vals[j]) * 
					( (dN_p_xi_vals[i]*w_func - N_p_xi_vals[i]*w_func_del_xi) / 
						(pow(w_func, 2.0)) );

				R_ij_pq_del_eta = (w_ij*N_p_xi_vals[i]) * 
					( (dN_q_eta_vals[j]*w_func - N_q_eta_vals[j]*w_func_del_eta) / 
						(pow(w_func, 2.0)) );

				// Place the value at the correct point in the 4 dimension multiarray
				basis_grad_vals_i[i + j*num_i_basis + k*num_i_basis*num_j_basis + 0*num_i_basis*num_j_basis*num_pts] = R_ij_pq_del_xi;
				basis_grad_vals_i[i + j*num_i_basis + k*num_i_basis*num_j_basis + 1*num_i_basis*num_j_basis*num_pts] = R_ij_pq_del_eta;
			
			}
		}
	}

	return (const struct const_Multiarray_d*)basis_grad_vals;

}


// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static double bernstein_std (const int p, const int i, const double r)
{
	if ((i == -1) || (p-i == -1))
		return 0.0;

	return binomial_coef(p,i)*pow(0.5*(1.0+r),i)*pow(0.5*(1.0-r),p-i);
}

static double grad_bernstein_std (const int p, const int i, const double r)
{
	return 0.5*p*(bernstein_std(p-1,i-1,r) - bernstein_std(p-1,i,r));
}
