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
#include "gsl/gsl_math.h"
#include "gsl/gsl_sf_gamma.h"

#include "macros.h"
#include "definitions_tol.h"
#include "definitions_elements.h"

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
	assert(rst->layout == 'C');

	const ptrdiff_t d   = rst->ext_1,
	                n_n = rst->ext_0,
	                pp1 = p_b+1,
	                n_b = compute_n_basis(d,p_b,ST_TP);

	struct Multiarray_Matrix_d* grad_phi_rst = constructor_empty_Multiarray_Matrix_d(false,1,&d); // returned

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

	const ptrdiff_t d   = rst->ext_1,
	                n_n = rst->ext_0,
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

	const ptrdiff_t d   = rst->ext_1,
	                n_n = rst->ext_0,
	                n_b = compute_n_basis(d,p_b,ST_SI);

	assert(!(d < 2 || d > 3));

	const struct const_Matrix_d*const abc = constructor_abc_from_rst_si(rst); // destructed
	const double*const a = get_col_const_Matrix_d(0,abc),
	            *const b = get_col_const_Matrix_d(1,abc),
	            *const c = ( d > 2 ? get_col_const_Matrix_d(2,abc) : NULL);

	struct Multiarray_Matrix_d* grad_phi_rst = constructor_empty_Multiarray_Matrix_d(false,1,&d); // returned

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
		grad_phi_rst_part[1][0] *=  2.0/3.0*sqrt(3.0)*a_n;
		grad_phi_rst_part[1][1] *=  2.0/3.0*sqrt(3.0)    ;
		grad_phi_rst_part[1][2] *= -2.0/3.0*sqrt(3.0)*i  ;
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
			grad_phi_rst_part[2][0] *=  2.0/3.0*sqrt(6.0)*a_n            ;
			grad_phi_rst_part[2][1] *=  1.0/6.0*sqrt(6.0)*(3.0*b_n+1.0)  ;
			grad_phi_rst_part[2][2] *= -1.0/6.0*sqrt(6.0)*(3.0*b_n+1.0)*i;
			grad_phi_rst_part[2][3] *=  1.0/2.0*sqrt(6.0)                ;
			grad_phi_rst_part[2][4] *= -1.0/2.0*sqrt(6.0)*(i+j)          ;
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

	const ptrdiff_t d   = rst->ext_1,
	                n_n = rst->ext_0,
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

	const ptrdiff_t d   = rst->ext_1,
	                n_n = rst->ext_0,
	                n_b = compute_n_basis(d,p_b,ST_PYR);

	assert(d == 3);

	const struct const_Matrix_d*const abc = constructor_abc_from_rst_pyr(rst); // destructed
	const double*const a = get_col_const_Matrix_d(0,abc),
	            *const b = get_col_const_Matrix_d(1,abc),
	            *const c = get_col_const_Matrix_d(2,abc);

	struct Multiarray_Matrix_d* grad_phi_rst = constructor_empty_Multiarray_Matrix_d(false,1,&d); // returned

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
	for (int i = 0; i < pp1; ++i) {
		           *phi_data  = bernstein_std(p_b,i,r[n]);
		if (d > 1) *phi_data *= bernstein_std(p_b,j,s[n]);
		if (d > 2) *phi_data *= bernstein_std(p_b,k,t[n]);
		++phi_data;
	}}}}

	return (const struct const_Matrix_d*) phi_rst;
}

const struct const_Multiarray_Matrix_d* constructor_grad_basis_tp_bezier
	(const int p_b, const struct const_Matrix_d*const rst)
{
	assert(rst->layout == 'C');

	const ptrdiff_t d   = rst->ext_1,
	                n_n = rst->ext_0,
	                pp1 = p_b+1,
	                n_b = compute_n_basis(d,p_b,ST_TP);

	struct Multiarray_Matrix_d* grad_phi_rst = constructor_empty_Multiarray_Matrix_d(false,1,&d); // returned

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
	assert(rst->layout == 'C');

	const ptrdiff_t d   = rst->ext_1,
	                n_n = rst->ext_0;

	assert(!(d < 2 || d > 3));

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
	assert(rst->layout == 'C');

	const ptrdiff_t d   = rst->ext_1,
	                n_n = rst->ext_0;

	assert(d == 3);

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

basis_fptr get_basis_by_super_type (const int s_type, const char*const ref_basis_name)
{
	if (s_type == ST_TP) {
		if (strcmp(ref_basis_name,"ortho") == 0)
			return  constructor_basis_tp_orthonormal;
		else if (strcmp(ref_basis_name,"bezier") == 0)
			return  constructor_basis_tp_bezier;
	} else if (s_type == ST_SI) {
		return  constructor_basis_si_orthonormal;
	} else if (s_type == ST_PYR) {
		return  constructor_basis_pyr_orthonormal;
	}

	EXIT_ERROR("Did not find the basis with the specified inputs: (%d, %s)\n",s_type,ref_basis_name);
}

grad_basis_fptr get_grad_basis_by_super_type (const int s_type, const char*const ref_grad_basis_name)
{
	if (s_type == ST_TP) {
		if (strcmp(ref_grad_basis_name,"ortho") == 0)
			return  constructor_grad_basis_tp_orthonormal;
		else if (strcmp(ref_grad_basis_name,"bezier") == 0)
			return  constructor_grad_basis_tp_bezier;
	} else if (s_type == ST_SI) {
		return  constructor_grad_basis_si_orthonormal;
	} else if (s_type == ST_PYR) {
		return  constructor_grad_basis_pyr_orthonormal;
	}

	EXIT_ERROR("Did not find the basis with the specified inputs: (%d, %s)\n",s_type,ref_grad_basis_name);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static double bernstein_std (const int p, const int i, const double r)
{
	if ((i == -1) || (p-i == -1))
		return 0.0;

	return binomial_coef(p,i)*pow(0.5*(1.0-r),p-i)*pow(0.5*(1.0+r),i);
}

static double grad_bernstein_std (const int p, const int i, const double r)
{
	return 0.5*p*(bernstein_std(p-1,i-1,r) - bernstein_std(p-1,i,r));
}
