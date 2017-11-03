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
/** \file
 */

#include "test_support_bases.h"

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include "gsl/gsl_math.h"
#include "gsl/gsl_poly.h"

#include "test_support_matrix.h"

#include "macros.h"
#include "definitions_core.h"
#include "definitions_elements.h"
#include "definitions_nodes.h"

#include "multiarray.h"
#include "matrix.h"
#include "vector.h"

#include "bases.h"
#include "nodes.h"

// Static function declarations ************************************************************************************* //

/// \brief Set the scaling parameters for the tri orthonormal basis functions.
static void set_scaling_basis_tri
	(const int i,    ///< Index of the `a` polynomial component.
	 const int j,    ///< Index of the `b` polynomial component.
	 const double b, ///< Value of `b` (from abc coordinates).
	 double* con_i,  ///< Constant multiplier for the `a` polynomial component.
	 double* con_j,  ///< Constant multiplier for the `b` polynomial component.
	 double* con_b   ///< Constant multiplier for the additional `b` component.
	);

/// \brief Set the scaling parameters for the tet orthonormal basis functions.
static void set_scaling_basis_tet
	(const int i,    ///< Index of the `a` polynomial component.
	 const int j,    ///< Index of the `b` polynomial component.
	 const int k,    ///< Index of the `c` polynomial component.
	 const double b, ///< Value of `b` (from abc coordinates).
	 const double c, ///< Value of `c` (from abc coordinates).
         double* con_i,  ///< Constant multiplier for the `a` polynomial component.
	 double* con_j,  ///< Constant multiplier for the `b` polynomial component.
	 double* con_k,  ///< Constant multiplier for the `c` polynomial component.
	 double* con_b,  ///< Constant multiplier for the additional `b` component.
	 double* con_c   ///< Constant multiplier for the additional `c` component.
	);

/// \brief Set the scaling parameters for the pyramid orthonormal basis functions.
static void set_scaling_basis_pyr
	(const int i,    ///< Index of the `a` polynomial component.
	 const int j,    ///< Index of the `b` polynomial component.
	 const int k,    ///< Index of the `c` polynomial component.
	 const double c, ///< Value of `c` (from abc coordinates).
         double* con_i,  ///< Constant multiplier for the `a` polynomial component.
	 double* con_j,  ///< Constant multiplier for the `b` polynomial component.
	 double* con_k,  ///< Constant multiplier for the `c` polynomial component.
	 double* con_c   ///< Constant multiplier for the additional `c` component.
	);

/** \brief Evaluate the hard-coded tensor-product polynomial of order 2 in each direction at the input coordinates.
 *  \return The polynomial value for `derivative_index = 0` and the associated derivative otherwise.
 *
 *  The order of the polynomial is such that the following order bases are required for its exact representation:
 *  - tensor-product: line (2), quad (2), hex (2);
 *  - simplex:         tri (4), tet (6);
 *  - pyramid:         pyr (6);
 *
 *  The test for the gradient evaluation will fail if lower orders than those specified above are used.
 */
static double poly_rst
	(const double r,            ///< Reference coordinate in the r-direction.
	 const double s,            ///< Reference coordinate in the s-direction.
	 const double t,            ///< Reference coordinate in the t-direction.
	 const int derivative_index ///< The index for the "derivative" to return. See return.
	);

/** \brief Constructor for the polynomial function evaluated at the input coordinates.
 *  \return See brief. */
static const struct const_Multiarray_d* constructor_f_rst
	(const struct const_Matrix_d* rst ///< Input coordinates.
	);

/** \brief Constructor for the gradient of the polynomial function evaluated at the input coordinates.
 *  \return See brief. */
static const struct const_Multiarray_d* constructor_grad_f_rst
	(const struct const_Matrix_d* rst ///< Input coordinates.
	);

// Interface functions ********************************************************************************************** //
// Constructor functions ******************************************************************************************** //
// Tensor-Product Orthonormal *************************************************************************************** //

const struct const_Matrix_d* constructor_basis_tp_orthonormal_def (const int p_b, const struct const_Matrix_d*const rst)
{
	assert(rst->layout == 'C');

	const ptrdiff_t d   = rst->ext_1,
	                n_n = rst->ext_0,
	                n_b = compute_n_basis(d,p_b,ST_TP);

	struct Matrix_d* phi_rst = constructor_empty_Matrix_d('R',n_n,n_b); // returned
	double* phi_data = phi_rst->data;

	const double*const r = get_col_const_Matrix_d(0,rst),
	            *const s = ( d > 1 ? get_col_const_Matrix_d(1,rst) : NULL),
	            *const t = ( d > 2 ? get_col_const_Matrix_d(2,rst) : NULL);

	for (ptrdiff_t i = 0; i < n_n; ++i) {
		const double r_i = r[i],
		             s_i = ( d > 1 ? s[i] : 0.0),
		             t_i = ( d > 2 ? t[i] : 0.0);

		if (d == 1 && p_b == 3) {
			*phi_data++ = sqrt(1.0/2.0)*1.0;
			*phi_data++ = sqrt(3.0/2.0)*r_i;
			*phi_data++ = sqrt(5.0/2.0)*1.0/2.0*(3.0*pow(r_i,2.0)-1);
			*phi_data++ = sqrt(7.0/2.0)*1.0/2.0*(5.0*pow(r_i,3.0)-3.0*r_i);
		} else if (d == 2 && p_b == 2) {
			*phi_data++ = sqrt(1.0/2.0)*1.0*sqrt(1.0/2.0)*1.0;
			*phi_data++ = sqrt(3.0/2.0)*r_i*sqrt(1.0/2.0)*1.0;
			*phi_data++ = sqrt(5.0/2.0)*1.0/2.0*(3.0*pow(r_i,2.0)-1.0)*sqrt(1.0/2.0)*1.0;
			*phi_data++ = sqrt(1.0/2.0)*1.0*sqrt(3.0/2.0)*s_i;
			*phi_data++ = sqrt(3.0/2.0)*r_i*sqrt(3.0/2.0)*s_i;
			*phi_data++ = sqrt(5.0/2.0)*1.0/2.0*(3.0*pow(r_i,2.0)-1.0)*sqrt(3.0/2.0)*s_i;
			*phi_data++ = sqrt(1.0/2.0)*1.0*sqrt(5.0/2.0)*1.0/2.0*(3.0*pow(s_i,2.0)-1.0);
			*phi_data++ = sqrt(3.0/2.0)*r_i*sqrt(5.0/2.0)*1.0/2.0*(3.0*pow(s_i,2.0)-1.0);
			*phi_data++ = sqrt(5.0/2.0)*1.0/2.0*(3.0*pow(r_i,2.0)-1.0)
			             *sqrt(5.0/2.0)*1.0/2.0*(3.0*pow(s_i,2.0)-1.0);
		} else if (d == 3 && p_b == 1) {
			*phi_data++ = sqrt(1.0/2.0)*1.0*sqrt(1.0/2.0)*1.0*sqrt(1.0/2.0)*1.0;
			*phi_data++ = sqrt(3.0/2.0)*r_i*sqrt(1.0/2.0)*1.0*sqrt(1.0/2.0)*1.0;
			*phi_data++ = sqrt(1.0/2.0)*1.0*sqrt(3.0/2.0)*s_i*sqrt(1.0/2.0)*1.0;
			*phi_data++ = sqrt(3.0/2.0)*r_i*sqrt(3.0/2.0)*s_i*sqrt(1.0/2.0)*1.0;
			*phi_data++ = sqrt(1.0/2.0)*1.0*sqrt(1.0/2.0)*1.0*sqrt(3.0/2.0)*t_i;
			*phi_data++ = sqrt(3.0/2.0)*r_i*sqrt(1.0/2.0)*1.0*sqrt(3.0/2.0)*t_i;
			*phi_data++ = sqrt(1.0/2.0)*1.0*sqrt(3.0/2.0)*s_i*sqrt(3.0/2.0)*t_i;
			*phi_data++ = sqrt(3.0/2.0)*r_i*sqrt(3.0/2.0)*s_i*sqrt(3.0/2.0)*t_i;
		} else {
			EXIT_UNSUPPORTED;
		}
	}

	return (const struct const_Matrix_d*) phi_rst;
}

const struct const_Multiarray_Matrix_d* constructor_grad_basis_tp_orthonormal_def
	(const int p_b, const struct const_Matrix_d*const rst)
{
	assert(rst->layout == 'C');

	const ptrdiff_t d   = rst->ext_1,
	                n_n = rst->ext_0,
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

	for (ptrdiff_t i = 0; i < n_n; ++i) {
		const double r_i = r[i],
		             s_i = ( d > 1 ? s[i] : 0.0),
		             t_i = ( d > 2 ? t[i] : 0.0);

		if (d == 1 && p_b == 3) {
			*grad_phi_data[0]++ = 0.0;
			*grad_phi_data[0]++ = sqrt(3.0/2.0);
			*grad_phi_data[0]++ = sqrt(5.0/2.0)*3.0*r_i;
			*grad_phi_data[0]++ = sqrt(7.0/2.0)*1.0/2.0*(15.0*pow(r_i,2.0)-3.0);
		} else if (d == 2 && p_b == 2) {
			*grad_phi_data[0]++ = 0.0;
			*grad_phi_data[0]++ = sqrt(3.0/2.0)*sqrt(1.0/2.0)*1.0;
			*grad_phi_data[0]++ = sqrt(5.0/2.0)*3.0*r_i*sqrt(1.0/2.0)*1.0;
			*grad_phi_data[0]++ = 0.0;
			*grad_phi_data[0]++ = sqrt(3.0/2.0)*sqrt(3.0/2.0)*s_i;
			*grad_phi_data[0]++ = sqrt(5.0/2.0)*3.0*r_i*sqrt(3.0/2.0)*s_i;
			*grad_phi_data[0]++ = 0.0;
			*grad_phi_data[0]++ = sqrt(3.0/2.0)*sqrt(5.0/2.0)*1.0/2.0*(3.0*pow(s_i,2.0)-1.0);
			*grad_phi_data[0]++ = sqrt(5.0/2.0)*3.0*r_i*sqrt(5.0/2.0)*1.0/2.0*(3.0*pow(s_i,2.0)-1.0);

			*grad_phi_data[1]++ = 0.0;
			*grad_phi_data[1]++ = 0.0;
			*grad_phi_data[1]++ = 0.0;
			*grad_phi_data[1]++ = sqrt(1.0/2.0)*1.0*sqrt(3.0/2.0);
			*grad_phi_data[1]++ = sqrt(3.0/2.0)*r_i*sqrt(3.0/2.0);
			*grad_phi_data[1]++ = sqrt(5.0/2.0)*1.0/2.0*(3.0*pow(r_i,2.0)-1.0)*sqrt(3.0/2.0);
			*grad_phi_data[1]++ = sqrt(1.0/2.0)*1.0*sqrt(5.0/2.0)*3.0*s_i;
			*grad_phi_data[1]++ = sqrt(3.0/2.0)*r_i*sqrt(5.0/2.0)*3.0*s_i;
			*grad_phi_data[1]++ = sqrt(5.0/2.0)*1.0/2.0*(3.0*pow(r_i,2.0)-1.0)*sqrt(5.0/2.0)*3.0*s_i;
		} else if (d == 3 && p_b == 1) {
			*grad_phi_data[0]++ = sqrt(1.0/2.0)*0.0*sqrt(1.0/2.0)*1.0*sqrt(1.0/2.0)*1.0;
			*grad_phi_data[0]++ = sqrt(3.0/2.0)*1.0*sqrt(1.0/2.0)*1.0*sqrt(1.0/2.0)*1.0;
			*grad_phi_data[0]++ = sqrt(1.0/2.0)*0.0*sqrt(3.0/2.0)*s_i*sqrt(1.0/2.0)*1.0;
			*grad_phi_data[0]++ = sqrt(3.0/2.0)*1.0*sqrt(3.0/2.0)*s_i*sqrt(1.0/2.0)*1.0;
			*grad_phi_data[0]++ = sqrt(1.0/2.0)*0.0*sqrt(1.0/2.0)*1.0*sqrt(3.0/2.0)*t_i;
			*grad_phi_data[0]++ = sqrt(3.0/2.0)*1.0*sqrt(1.0/2.0)*1.0*sqrt(3.0/2.0)*t_i;
			*grad_phi_data[0]++ = sqrt(1.0/2.0)*0.0*sqrt(3.0/2.0)*s_i*sqrt(3.0/2.0)*t_i;
			*grad_phi_data[0]++ = sqrt(3.0/2.0)*1.0*sqrt(3.0/2.0)*s_i*sqrt(3.0/2.0)*t_i;

			*grad_phi_data[1]++ = sqrt(1.0/2.0)*1.0*sqrt(1.0/2.0)*0.0*sqrt(1.0/2.0)*1.0;
			*grad_phi_data[1]++ = sqrt(3.0/2.0)*r_i*sqrt(1.0/2.0)*0.0*sqrt(1.0/2.0)*1.0;
			*grad_phi_data[1]++ = sqrt(1.0/2.0)*1.0*sqrt(3.0/2.0)*1.0*sqrt(1.0/2.0)*1.0;
			*grad_phi_data[1]++ = sqrt(3.0/2.0)*r_i*sqrt(3.0/2.0)*1.0*sqrt(1.0/2.0)*1.0;
			*grad_phi_data[1]++ = sqrt(1.0/2.0)*1.0*sqrt(1.0/2.0)*0.0*sqrt(3.0/2.0)*t_i;
			*grad_phi_data[1]++ = sqrt(3.0/2.0)*r_i*sqrt(1.0/2.0)*0.0*sqrt(3.0/2.0)*t_i;
			*grad_phi_data[1]++ = sqrt(1.0/2.0)*1.0*sqrt(3.0/2.0)*1.0*sqrt(3.0/2.0)*t_i;
			*grad_phi_data[1]++ = sqrt(3.0/2.0)*r_i*sqrt(3.0/2.0)*1.0*sqrt(3.0/2.0)*t_i;

			*grad_phi_data[2]++ = sqrt(1.0/2.0)*1.0*sqrt(1.0/2.0)*1.0*sqrt(1.0/2.0)*0.0;
			*grad_phi_data[2]++ = sqrt(3.0/2.0)*r_i*sqrt(1.0/2.0)*1.0*sqrt(1.0/2.0)*0.0;
			*grad_phi_data[2]++ = sqrt(1.0/2.0)*1.0*sqrt(3.0/2.0)*s_i*sqrt(1.0/2.0)*0.0;
			*grad_phi_data[2]++ = sqrt(3.0/2.0)*r_i*sqrt(3.0/2.0)*s_i*sqrt(1.0/2.0)*0.0;
			*grad_phi_data[2]++ = sqrt(1.0/2.0)*1.0*sqrt(1.0/2.0)*1.0*sqrt(3.0/2.0)*1.0;
			*grad_phi_data[2]++ = sqrt(3.0/2.0)*r_i*sqrt(1.0/2.0)*1.0*sqrt(3.0/2.0)*1.0;
			*grad_phi_data[2]++ = sqrt(1.0/2.0)*1.0*sqrt(3.0/2.0)*s_i*sqrt(3.0/2.0)*1.0;
			*grad_phi_data[2]++ = sqrt(3.0/2.0)*r_i*sqrt(3.0/2.0)*s_i*sqrt(3.0/2.0)*1.0;
		} else {
			EXIT_UNSUPPORTED;
		}
	}

	return (const struct const_Multiarray_Matrix_d*) grad_phi_rst;
}

// Simplex Orthonormal ********************************************************************************************** //

const struct const_Matrix_d* constructor_basis_si_orthonormal_def (const int p_b, const struct const_Matrix_d*const rst)
{
	assert(rst->layout == 'C');

	const ptrdiff_t d   = rst->ext_1,
	                n_n = rst->ext_0,
	                n_b = compute_n_basis(d,p_b,ST_SI);

	const struct const_Matrix_d*const abc = constructor_abc_from_rst_si(rst); // destructed
	const double*const a = get_col_const_Matrix_d(0,abc),
	            *const b = get_col_const_Matrix_d(1,abc),
	            *const c = ( d > 2 ? get_col_const_Matrix_d(2,abc) : NULL);

	struct Matrix_d* phi_rst = constructor_empty_Matrix_d('R',n_n,n_b); // returned
	double* phi_data = phi_rst->data;

	double con_i = 0.0,
	       con_j = 0.0,
	       con_k = 0.0,
	       con_b = 0.0,
	       con_c = 0.0;
	for (ptrdiff_t n = 0; n < n_n; ++n) {
		const double a_n = a[n],
		             b_n = b[n],
		             c_n = ( d == 3 ? c[n] : -1.0 );

		if (d == 2 && p_b == 2) {
			double con = 2.0/pow(3.0,0.25);

			set_scaling_basis_tri(0,0,b_n,&con_i,&con_j,&con_b);
			*phi_data++ = con*con_b
			            * con_i*1.0
			            * con_j*1.0;

			set_scaling_basis_tri(0,1,b_n,&con_i,&con_j,&con_b);
			*phi_data++ = con*con_b
			            * con_i*1.0
			            * con_j*1.0/2.0*(3.0*b_n+1.0);

			set_scaling_basis_tri(0,2,b_n,&con_i,&con_j,&con_b);
			*phi_data++ = con*con_b
			            * con_i*1.0
			            * con_j*(5.0/2.0*pow(b_n-1.0,2.0)+6.0*(b_n-1.0)+3.0);

			set_scaling_basis_tri(1,0,b_n,&con_i,&con_j,&con_b);
			*phi_data++ = con*con_b
			            * con_i*a_n
			            * con_j*1.0;

			set_scaling_basis_tri(1,1,b_n,&con_i,&con_j,&con_b);
			*phi_data++ = con*con_b
			            * con_i*a_n
			            * con_j*1.0/2.0*(5.0*b_n+3.0);

			set_scaling_basis_tri(2,0,b_n,&con_i,&con_j,&con_b);
			*phi_data++ = con*con_b
			            * con_i*1.0/2.0*(3.0*pow(a_n,2.0)-1.0)
			            * con_j*1.0;
		} else if (d == 2 && p_b == 3) {
			double con = 2.0/pow(3.0,0.25);

			set_scaling_basis_tri(0,0,b_n,&con_i,&con_j,&con_b);
			*phi_data++ = con*con_b
			            * con_i*1.0
			            * con_j*1.0;

			set_scaling_basis_tri(0,1,b_n,&con_i,&con_j,&con_b);
			*phi_data++ = con*con_b
			            * con_i*1.0
			            * con_j*1.0/2.0*(3.0*b_n+1.0);

			set_scaling_basis_tri(0,2,b_n,&con_i,&con_j,&con_b);
			*phi_data++ = con*con_b
			            * con_i*1.0
			            * con_j*(5.0/2.0*pow(b_n-1.0,2.0)+6.0*(b_n-1.0)+3.0);

			set_scaling_basis_tri(0,3,b_n,&con_i,&con_j,&con_b);
			*phi_data++ = con*con_b
			            * con_i*1.0
			            * con_j*(35.0/8.0*pow(b_n-1.0,3.0)+15.0*pow(b_n-1.0,2.0)+15.0*(b_n-1)+4.0);

			set_scaling_basis_tri(1,0,b_n,&con_i,&con_j,&con_b);
			*phi_data++ = con*con_b
			            * con_i*a_n
			            * con_j*1.0;

			set_scaling_basis_tri(1,1,b_n,&con_i,&con_j,&con_b);
			*phi_data++ = con*con_b
			            * con_i*a_n
			            * con_j*1.0/2.0*(5.0*b_n+3.0);

			set_scaling_basis_tri(1,2,b_n,&con_i,&con_j,&con_b);
			*phi_data++ = con*con_b
			            * con_i*a_n
			            * con_j*(21.0/4.0*pow(b_n-1.0,2.0)+15.0*(b_n-1.0)+10.0);

			set_scaling_basis_tri(2,0,b_n,&con_i,&con_j,&con_b);
			*phi_data++ = con*con_b
			            * con_i*1.0/2.0*(3.0*pow(a_n,2.0)-1.0)
			            * con_j*1.0;

			set_scaling_basis_tri(2,1,b_n,&con_i,&con_j,&con_b);
			*phi_data++ = con*con_b
			            * con_i*1.0/2.0*(3.0*pow(a_n,2.0)-1.0)
			            * con_j*1.0/2.0*(7.0*b_n+5.0);

			set_scaling_basis_tri(3,0,b_n,&con_i,&con_j,&con_b);
			*phi_data++ = con*con_b
			            * con_i*1.0/2.0*(5.0*pow(a_n,3.0)-3.0*a_n)
			            * con_j*1.0;
		} else if (d == 3 && p_b == 2) {
			const double con = 4.0/pow(2.0,0.25);

			set_scaling_basis_tet(0,0,0,b_n,c_n,&con_i,&con_j,&con_k,&con_b,&con_c);
			*phi_data++ = con*con_b*con_c
			            * con_i*1.0
			            * con_j*1.0
			            * con_k*1.0;

			set_scaling_basis_tet(0,0,1,b_n,c_n,&con_i,&con_j,&con_k,&con_b,&con_c);
			*phi_data++ = con*con_b*con_c
			            * con_i*1.0
			            * con_j*1.0
			            * con_k*1.0/2.0*(4.0*c_n+2.0);

			set_scaling_basis_tet(0,0,2,b_n,c_n,&con_i,&con_j,&con_k,&con_b,&con_c);
			*phi_data++ = con*con_b*con_c
			            * con_i*1.0
			            * con_j*1.0
			            * con_k*(15.0/4.0*pow(c_n-1.0,2.0)+10.0*(c_n-1.0)+6.0);

			set_scaling_basis_tet(0,1,0,b_n,c_n,&con_i,&con_j,&con_k,&con_b,&con_c);
			*phi_data++ = con*con_b*con_c
			            * con_i*1.0
			            * con_j*1.0/2.0*(3.0*b_n+1.0)
			            * con_k*1.0;

			set_scaling_basis_tet(0,1,1,b_n,c_n,&con_i,&con_j,&con_k,&con_b,&con_c);
			*phi_data++ = con*con_b*con_c
			            * con_i*1.0
			            * con_j*1.0/2.0*(3.0*b_n+1.0)
			            * con_k*1.0/2.0*(6.0*c_n+4.0);

			set_scaling_basis_tet(0,2,0,b_n,c_n,&con_i,&con_j,&con_k,&con_b,&con_c);
			*phi_data++ = con*con_b*con_c
			            * con_i*1.0
			            * con_j*(5.0/2.0*pow(b_n-1.0,2.0)+6.0*(b_n-1.0)+3.0)
			            * con_k*1.0;

			set_scaling_basis_tet(1,0,0,b_n,c_n,&con_i,&con_j,&con_k,&con_b,&con_c);
			*phi_data++ = con*con_b*con_c
			            * con_i*a_n
			            * con_j*1.0
			            * con_k*1.0;

			set_scaling_basis_tet(1,0,1,b_n,c_n,&con_i,&con_j,&con_k,&con_b,&con_c);
			*phi_data++ = con*con_b*con_c
			            * con_i*a_n
			            * con_j*1.0
			            * con_k*1.0/2.0*(6.0*c_n+4.0);

			set_scaling_basis_tet(1,1,0,b_n,c_n,&con_i,&con_j,&con_k,&con_b,&con_c);
			*phi_data++ = con*con_b*con_c
			            * con_i*a_n
			            * con_j*1.0/2.0*(5.0*b_n+3.0)
			            * con_k*1.0;

			set_scaling_basis_tet(2,0,0,b_n,c_n,&con_i,&con_j,&con_k,&con_b,&con_c);
			*phi_data++ = con*con_b*con_c
			            * con_i*1.0/2.0*(3.0*pow(a_n,2.0)-1.0)
			            * con_j*1.0
			            * con_k*1.0;
		} else {
			EXIT_UNSUPPORTED;
		}
	}
	destructor_const_Matrix_d(abc);

	return (const struct const_Matrix_d*) phi_rst;
}

const struct const_Multiarray_Matrix_d* constructor_grad_basis_si_orthonormal_def
	(const int p_b, const struct const_Matrix_d*const rst)
{
	assert(rst->layout == 'C');

	const ptrdiff_t d   = rst->ext_1,
	                n_n = rst->ext_0,
	                n_b = compute_n_basis(d,p_b,ST_SI);

	const struct const_Matrix_d*const abc = constructor_abc_from_rst_si(rst); // destructed
	const double*const a = get_col_const_Matrix_d(0,abc),
	            *const b = get_col_const_Matrix_d(1,abc),
	            *const c = ( d > 2 ? get_col_const_Matrix_d(2,abc) : NULL);

	struct Multiarray_Matrix_d* grad_phi_rst = constructor_empty_Multiarray_Matrix_d(false,1,&d); // returned

	double* grad_phi_data[d];
	for (int dim = 0; dim < d; ++dim) {
		grad_phi_rst->data[dim] = constructor_empty_Matrix_d('R',n_n,n_b); // keep
		grad_phi_data[dim] = grad_phi_rst->data[dim]->data;
	}

	double con_i = 0.0,
	       con_j = 0.0,
	       con_k = 0.0,
	       con_b = 0.0,
	       con_c = 0.0;
	for (ptrdiff_t n = 0; n < n_n; ++n) {
		const double a_n = a[n],
		             b_n = b[n],
		             c_n = ( d == 3 ? c[n] : -1.0 );

		int i = 0,
		    j = 0,
		    k = 0;
		if (d == 2 && p_b == 2) {
			const double con = 2.0/pow(3.0,0.25);
			i = 0; j = 0;
			set_scaling_basis_tri(i,j,b_n,&con_i,&con_j,&con_b);
			*grad_phi_data[0]++ = 0.0;
			*grad_phi_data[1]++ = 0.0;

			i = 0; j = 1;
			set_scaling_basis_tri(i,j,b_n,&con_i,&con_j,&con_b);
			*grad_phi_data[0]++ = 0.0;
			*grad_phi_data[1]++ = con*( 2.0/3.0*sqrt(3.0)*pow(1.0-b_n,i)
			                           *con_i*(1.0)
			                           *con_j*(3.0/2.0) );

			i = 0; j = 2;
			set_scaling_basis_tri(i,j,b_n,&con_i,&con_j,&con_b);
			*grad_phi_data[0]++ = 0.0;
			*grad_phi_data[1]++ = con*( 2.0/3.0*sqrt(3.0)*pow(1.0-b_n,i)
			                           *con_i*(1.0)
			                           *con_j*(5.0*(b_n-1.0)+6.0) );

			i = 1; j = 0;
			set_scaling_basis_tri(i,j,b_n,&con_i,&con_j,&con_b);
			*grad_phi_data[0]++ = con*2.0*pow(1.0-b_n,i-1.0)
			                     *con_i*(1.0)
			                     *con_j*(1.0);
			*grad_phi_data[1]++ = 0.0;

			i = 1; j = 1;
			set_scaling_basis_tri(i,j,b_n,&con_i,&con_j,&con_b);
			*grad_phi_data[0]++ = con*2.0*pow(1.0-b_n,i-1.0)
			                     *con_i*(1.0)
			                     *con_j*(1.0/2.0*(5.0*b_n+3.0));
			*grad_phi_data[1]++ = con*( -2.0/3.0*sqrt(3.0)*i*pow(1.0-b_n,i-1.0)
			                     *con_i*(a_n)
			                     *con_j*(1.0/2.0*(5.0*b_n+3.0))
			                     +2.0/3.0*sqrt(3.0)*a_n*pow(1.0-b_n,i-1.0)
			                     *con_i*(1.0)
			                     *con_j*(1.0/2.0*(5.0*b_n+3.0))
			                     +2.0/3.0*sqrt(3.0)*pow(1.0-b_n,i)
			                     *con_i*(a_n)
			                     *con_j*(5.0/2.0) );

			i = 2; j = 0;
			set_scaling_basis_tri(i,j,b_n,&con_i,&con_j,&con_b);
			*grad_phi_data[0]++ = con*2.0*pow(1.0-b_n,i-1.0)
			                     *con_i*(3.0*a_n)
			                     *con_j*(1.0);
			*grad_phi_data[1]++ = con*( -2.0/3.0*sqrt(3.0)*i*pow(1.0-b_n,i-1.0)
			                     *con_i*(1.0/2.0*(3.0*pow(a_n,2.0)-1.0))
			                     *con_j*(1.0)
			                     +2.0/3.0*sqrt(3.0)*a_n*pow(1.0-b_n,i-1.0)
			                     *con_i*(3.0*a_n)
			                     *con_j*(1.0) );
		} else if (d == 3 && p_b == 1) {
			const double con = 4.0/pow(2.0,0.25);

			i = 0; j = 0; k = 0;
			set_scaling_basis_tet(i,j,k,b_n,c_n,&con_i,&con_j,&con_k,&con_b,&con_c);
			*grad_phi_data[0]++ = 0.0;
			*grad_phi_data[1]++ = 0.0;
			*grad_phi_data[2]++ = 0.0;

			i = 0; j = 0; k = 1;
			set_scaling_basis_tet(i,j,k,b_n,c_n,&con_i,&con_j,&con_k,&con_b,&con_c);
			*grad_phi_data[0]++ = 0.0;
			*grad_phi_data[1]++ = 0.0;
			*grad_phi_data[2]++ = con*( 1.0/2.0*sqrt(6.0)*pow(1.0-b_n,i)*pow(1.0-c_n,i+j)
			                           *con_i*(1.0)
			                           *con_j*(1.0)
			                           *con_k*(2.0) );

			i = 0; j = 1; k = 0;
			set_scaling_basis_tet(i,j,k,b_n,c_n,&con_i,&con_j,&con_k,&con_b,&con_c);
			*grad_phi_data[0]++ = 0.0;
			*grad_phi_data[1]++ = con*( 4.0/3.0*sqrt(3.0)*pow(1.0-b_n,i)*pow(1.0-c_n,i+j-1.0)
			                           *con_i*(1.0)
			                           *con_j*(3.0/2.0)
			                           *con_k*(1.0) );
			*grad_phi_data[2]++ = 0.0;

			i = 1; j = 0; k = 0;
			set_scaling_basis_tet(i,j,k,b_n,c_n,&con_i,&con_j,&con_k,&con_b,&con_c);
			*grad_phi_data[0]++ = con*( 4.0*pow(1.0-b_n,i-1.0)*pow(1.0-c_n,i+j-1.0)
			                           *con_i*(1.0)
			                           *con_j*(1.0)
			                           *con_k*(1.0) );
			*grad_phi_data[1]++ = 0.0;
			*grad_phi_data[2]++ = 0.0;
		} else {
			EXIT_UNSUPPORTED;
		}
	}
	destructor_const_Matrix_d(abc);

	return (const struct const_Multiarray_Matrix_d*) grad_phi_rst;
}

// Pyramid Orthonormal ********************************************************************************************** //

const struct const_Matrix_d* constructor_basis_pyr_orthonormal_def
	(const int p_b, const struct const_Matrix_d*const rst)
{
	assert(rst->layout == 'C');

	const ptrdiff_t d   = rst->ext_1,
	                n_n = rst->ext_0,
	                n_b = compute_n_basis(d,p_b,ST_PYR);

	const struct const_Matrix_d*const abc = constructor_abc_from_rst_pyr(rst); // destructed
	const double*const a = get_col_const_Matrix_d(0,abc),
	            *const b = get_col_const_Matrix_d(1,abc),
	            *const c = get_col_const_Matrix_d(2,abc);

	struct Matrix_d* phi_rst = constructor_empty_Matrix_d('R',n_n,n_b); // returned
	double* phi_data = phi_rst->data;

	const double con = pow(2.0,1.25);
	double con_i = 0.0,
	       con_j = 0.0,
	       con_k = 0.0,
	       con_c = 0.0;
	for (ptrdiff_t n = 0; n < n_n; ++n) {
		const double a_n = a[n],
		             b_n = b[n],
		             c_n = c[n];

		if (d == 3 && p_b == 2) {
			set_scaling_basis_pyr(0,0,0,c_n,&con_i,&con_j,&con_k,&con_c);
			*phi_data++ = con*con_c
			            * con_i*(1.0)
			            * con_j*(1.0)
			            * con_k*(1.0);

			set_scaling_basis_pyr(0,0,1,c_n,&con_i,&con_j,&con_k,&con_c);
			*phi_data++ = con*con_c
			            * con_i*(1.0)
			            * con_j*(1.0)
			            * con_k*(1.0/2.0*(4.0*c_n+2.0));

			set_scaling_basis_pyr(0,0,2,c_n,&con_i,&con_j,&con_k,&con_c);
			*phi_data++ = con*con_c
			            * con_i*(1.0)
			            * con_j*(1.0)
			            * con_k*(15.0/4.0*pow(c_n-1.0,2.0)+10.0*(c_n-1.0)+6.0);

			set_scaling_basis_pyr(0,1,0,c_n,&con_i,&con_j,&con_k,&con_c);
			*phi_data++ = con*con_c
			            * con_i*(1.0)
			            * con_j*(b_n)
			            * con_k*(1.0);

			set_scaling_basis_pyr(0,1,1,c_n,&con_i,&con_j,&con_k,&con_c);
			*phi_data++ = con*con_c
			            * con_i*(1.0)
			            * con_j*(b_n)
			            * con_k*(1.0/2.0*(6.0*c_n+4.0));

			set_scaling_basis_pyr(0,2,0,c_n,&con_i,&con_j,&con_k,&con_c);
			*phi_data++ = con*con_c
			            * con_i*(1.0)
			            * con_j*(1.0/2.0*(3.0*pow(b_n,2.0)-1.0))
			            * con_k*(1.0);

			set_scaling_basis_pyr(1,0,0,c_n,&con_i,&con_j,&con_k,&con_c);
			*phi_data++ = con*con_c
			            * con_i*(a_n)
			            * con_j*(1.0)
			            * con_k*(1.0);

			set_scaling_basis_pyr(1,0,1,c_n,&con_i,&con_j,&con_k,&con_c);
			*phi_data++ = con*con_c
			            * con_i*(a_n)
			            * con_j*(1.0)
			            * con_k*(1.0/2.0*(6.0*c_n+4.0));

			set_scaling_basis_pyr(1,1,0,c_n,&con_i,&con_j,&con_k,&con_c);
			*phi_data++ = con*con_c
			            * con_i*(a_n)
			            * con_j*(b_n)
			            * con_k*(1.0);

			set_scaling_basis_pyr(1,1,1,c_n,&con_i,&con_j,&con_k,&con_c);
			*phi_data++ = con*con_c
			            * con_i*(a_n)
			            * con_j*(b_n)
			            * con_k*(1.0/2.0*(6.0*c_n+4.0));

			set_scaling_basis_pyr(1,2,0,c_n,&con_i,&con_j,&con_k,&con_c);
			*phi_data++ = con*con_c
			            * con_i*(a_n)
			            * con_j*(1.0/2.0*(3.0*pow(b_n,2.0)-1.0))
			            * con_k*(1.0);

			set_scaling_basis_pyr(2,0,0,c_n,&con_i,&con_j,&con_k,&con_c);
			*phi_data++ = con*con_c
			            * con_i*(1.0/2.0*(3.0*pow(a_n,2.0)-1.0))
			            * con_j*(1.0)
			            * con_k*(1.0);

			set_scaling_basis_pyr(2,1,0,c_n,&con_i,&con_j,&con_k,&con_c);
			*phi_data++ = con*con_c
			            * con_i*(1.0/2.0*(3.0*pow(a_n,2.0)-1.0))
			            * con_j*(b_n)
			            * con_k*(1.0);

			set_scaling_basis_pyr(2,2,0,c_n,&con_i,&con_j,&con_k,&con_c);
			*phi_data++ = con*con_c
			            * con_i*(1.0/2.0*(3.0*pow(a_n,2.0)-1.0))
			            * con_j*(1.0/2.0*(3.0*pow(b_n,2.0)-1.0))
			            * con_k*(1.0);
		} else {
			EXIT_UNSUPPORTED;
		}
	}
	destructor_const_Matrix_d(abc);

	return (const struct const_Matrix_d*) phi_rst;
}

const struct const_Multiarray_Matrix_d* constructor_grad_basis_pyr_orthonormal_def
	(const int p_b, const struct const_Matrix_d*const rst)
{
	assert(rst->layout == 'C');

	const ptrdiff_t d   = rst->ext_1,
	                n_n = rst->ext_0,
	                n_b = compute_n_basis(d,p_b,ST_PYR);

	const struct const_Matrix_d*const abc = constructor_abc_from_rst_pyr(rst); // destructed
	const double*const a = get_col_const_Matrix_d(0,abc),
	            *const b = get_col_const_Matrix_d(1,abc),
	            *const c = get_col_const_Matrix_d(2,abc);

	struct Multiarray_Matrix_d* grad_phi_rst = constructor_empty_Multiarray_Matrix_d(false,1,&d); // returned

	double* grad_phi_data[d];
	for (int dim = 0; dim < d; ++dim) {
		grad_phi_rst->data[dim] = constructor_empty_Matrix_d('R',n_n,n_b); // keep
		grad_phi_data[dim] = grad_phi_rst->data[dim]->data;
	}

	const double con = pow(2.0,1.25);

	double con_i = 0.0,
	       con_j = 0.0,
	       con_k = 0.0,
	       con_c = 0.0;
	for (ptrdiff_t n = 0; n < n_n; ++n) {
		const double a_n = a[n],
		             b_n = b[n],
		             c_n = c[n];

		if (d == 3 && p_b == 2) {
			set_scaling_basis_pyr(0,0,0,c_n,&con_i,&con_j,&con_k,&con_c);
			*grad_phi_data[0]++ = 0.0;
			*grad_phi_data[1]++ = 0.0;
			*grad_phi_data[2]++ = 0.0;

			set_scaling_basis_pyr(0,0,1,c_n,&con_i,&con_j,&con_k,&con_c);
			*grad_phi_data[0]++ = 0.0;
			*grad_phi_data[1]++ = 0.0;
			*grad_phi_data[2]++ = con*( sqrt(2.0)*pow(1.0-c_n,0.0)
			                          * con_i*(1.0)
			                          * con_j*(1.0)
			                          * con_k*(1.0/2.0*(4.0)));

			set_scaling_basis_pyr(0,0,2,c_n,&con_i,&con_j,&con_k,&con_c);
			*grad_phi_data[0]++ = 0.0;
			*grad_phi_data[1]++ = 0.0;
			*grad_phi_data[2]++ = con*( sqrt(2.0)*pow(1.0-c_n,0.0)
			                          * con_i*(1.0)
			                          * con_j*(1.0)
			                          * con_k*(15.0/4.0*2.0*pow(c_n-1.0,1.0)+10.0));

			set_scaling_basis_pyr(0,1,0,c_n,&con_i,&con_j,&con_k,&con_c);
			*grad_phi_data[0]++ = 0.0;
			*grad_phi_data[1]++ = con*( 2.0*pow(1.0-c_n,0.0)
			                          * con_i*(1.0)
			                          * con_j*(1.0)
			                          * con_k*(1.0));
			*grad_phi_data[2]++ = con*( sqrt(2.0)*b_n*pow(1.0-c_n,0.0)
			                          * con_i*(1.0)
			                          * con_j*(1.0)
			                          * con_k*(1.0)
			                          - 1.0*sqrt(2.0)*pow(1.0-c_n,0.0)
			                          * con_i*(1.0)
			                          * con_j*(b_n)
			                          * con_k*(1.0));

			set_scaling_basis_pyr(0,1,1,c_n,&con_i,&con_j,&con_k,&con_c);
			*grad_phi_data[0]++ = 0.0;
			*grad_phi_data[1]++ = con*( 2.0*pow(1.0-c_n,0.0)
			                          * con_i*(1.0)
			                          * con_j*(1.0)
			                          * con_k*(1.0/2.0*(6.0*c_n+4.0)));
			*grad_phi_data[2]++ = con*( sqrt(2.0)*b_n*pow(1.0-c_n,0.0)
			                          * con_i*(1.0)
			                          * con_j*(1.0)
			                          * con_k*(1.0/2.0*(6.0*c_n+4.0))
			                          + sqrt(2.0)*pow(1.0-c_n,1.0)
			                          * con_i*(1.0)
			                          * con_j*(b_n)
			                          * con_k*(1.0/2.0*(6.0))
			                          - 1.0*sqrt(2.0)*pow(1.0-c_n,0.0)
			                          * con_i*(1.0)
			                          * con_j*(b_n)
			                          * con_k*(1.0/2.0*(6.0*c_n+4.0)));

			set_scaling_basis_pyr(0,2,0,c_n,&con_i,&con_j,&con_k,&con_c);
			*grad_phi_data[0]++ = 0.0;
			*grad_phi_data[1]++ = con*( 2.0*pow(1.0-c_n,1.0)
			                          * con_i*(1.0)
			                          * con_j*(1.0/2.0*(3.0*2.0*b_n))
			                          * con_k*(1.0));
			*grad_phi_data[2]++ = con*( sqrt(2.0)*b_n*pow(1.0-c_n,1.0)
			                          * con_i*(1.0)
			                          * con_j*(1.0/2.0*(3.0*2.0*b_n))
			                          * con_k*(1.0)
			                          - 2.0*sqrt(2.0)*pow(1.0-c_n,1.0)
			                          * con_i*(1.0)
			                          * con_j*(1.0/2.0*(3.0*pow(b_n,2.0)-1.0))
			                          * con_k*(1.0));

			set_scaling_basis_pyr(1,0,0,c_n,&con_i,&con_j,&con_k,&con_c);
			*grad_phi_data[0]++ = con*( 2.0*pow(1.0-c_n,0.0)
			                          * con_i*(1.0)
			                          * con_j*(1.0)
			                          * con_k*(1.0));
			*grad_phi_data[1]++ = 0.0;
			*grad_phi_data[2]++ = con*( sqrt(2.0)*a_n*pow(1.0-c_n,0.0)
			                          * con_i*(1.0)
			                          * con_j*(1.0)
			                          * con_k*(1.0)
			                          - 1.0*sqrt(2.0)*pow(1.0-c_n,0.0)
			                          * con_i*(a_n)
			                          * con_j*(1.0)
			                          * con_k*(1.0));

			set_scaling_basis_pyr(1,0,1,c_n,&con_i,&con_j,&con_k,&con_c);
			*grad_phi_data[0]++ = con*( 2.0*pow(1.0-c_n,0.0)
			                          * con_i*(1.0)
			                          * con_j*(1.0)
			                          * con_k*(1.0/2.0*(6.0*c_n+4.0)));
			*grad_phi_data[1]++ = 0.0;
			*grad_phi_data[2]++ = con*( sqrt(2.0)*a_n*pow(1.0-c_n,0.0)
			                          * con_i*(1.0)
			                          * con_j*(1.0)
			                          * con_k*(1.0/2.0*(6.0*c_n+4.0))
			                          + sqrt(2.0)*pow(1.0-c_n,1.0)
			                          * con_i*(a_n)
			                          * con_j*(1.0)
			                          * con_k*(1.0/2.0*(6.0))
			                          - 1.0*sqrt(2.0)*pow(1.0-c_n,0.0)
			                          * con_i*(a_n)
			                          * con_j*(1.0)
			                          * con_k*(1.0/2.0*(6.0*c_n+4.0)));

			set_scaling_basis_pyr(1,1,0,c_n,&con_i,&con_j,&con_k,&con_c);
			*grad_phi_data[0]++ = con*( 2.0*pow(1.0-c_n,0.0)
			                          * con_i*(1.0)
			                          * con_j*(b_n)
			                          * con_k*(1.0));
			*grad_phi_data[1]++ = con*( 2.0*pow(1.0-c_n,0.0)
			                          * con_i*(a_n)
			                          * con_j*(1.0)
			                          * con_k*(1.0));
			*grad_phi_data[2]++ = con*( sqrt(2.0)*a_n*pow(1.0-c_n,0.0)
			                          * con_i*(1.0)
			                          * con_j*(b_n)
			                          * con_k*(1.0)
			                          + sqrt(2.0)*b_n*pow(1.0-c_n,0.0)
			                          * con_i*(a_n)
			                          * con_j*(1.0)
			                          * con_k*(1.0)
			                          - 1.0*sqrt(2.0)*pow(1.0-c_n,0.0)
			                          * con_i*(a_n)
			                          * con_j*(b_n)
			                          * con_k*(1.0));

			set_scaling_basis_pyr(1,1,1,c_n,&con_i,&con_j,&con_k,&con_c);
			*grad_phi_data[0]++ = con*( 2.0*pow(1.0-c_n,0.0)
			                          * con_i*(1.0)
			                          * con_j*(b_n)
			                          * con_k*(1.0/2.0*(6.0*c_n+4.0)));
			*grad_phi_data[1]++ = con*( 2.0*pow(1.0-c_n,0.0)
			                          * con_i*(a_n)
			                          * con_j*(1.0)
			                          * con_k*(1.0/2.0*(6.0*c_n+4.0)));
			*grad_phi_data[2]++ = con*( sqrt(2.0)*a_n*pow(1.0-c_n,0.0)
			                          * con_i*(1.0)
			                          * con_j*(b_n)
			                          * con_k*(1.0/2.0*(6.0*c_n+4.0))
			                          + sqrt(2.0)*b_n*pow(1.0-c_n,0.0)
			                          * con_i*(a_n)
			                          * con_j*(1.0)
			                          * con_k*(1.0/2.0*(6.0*c_n+4.0))
			                          + sqrt(2.0)*pow(1.0-c_n,1.0)
			                          * con_i*(a_n)
			                          * con_j*(b_n)
			                          * con_k*(1.0/2.0*(6.0))
			                          - 1.0*sqrt(2.0)*pow(1.0-c_n,0.0)
			                          * con_i*(a_n)
			                          * con_j*(b_n)
			                          * con_k*(1.0/2.0*(6.0*c_n+4.0)));

			set_scaling_basis_pyr(1,2,0,c_n,&con_i,&con_j,&con_k,&con_c);
			*grad_phi_data[0]++ = con*( 2.0*pow(1.0-c_n,1.0)
			                          * con_i*(1.0)
			                          * con_j*(1.0/2.0*(3.0*pow(b_n,2.0)-1.0))
			                          * con_k*(1.0));
			*grad_phi_data[1]++ = con*( 2.0*pow(1.0-c_n,1.0)
			                          * con_i*(a_n)
			                          * con_j*(1.0/2.0*(3.0*2.0*b_n))
			                          * con_k*(1.0));
			*grad_phi_data[2]++ = con*( sqrt(2.0)*a_n*pow(1.0-c_n,1.0)
			                          * con_i*(1.0)
			                          * con_j*(1.0/2.0*(3.0*pow(b_n,2.0)-1.0))
			                          * con_k*(1.0)
			                          + sqrt(2.0)*b_n*pow(1.0-c_n,1.0)
			                          * con_i*(a_n)
			                          * con_j*(1.0/2.0*(3.0*2.0*b_n))
			                          * con_k*(1.0)
			                          - 2.0*sqrt(2.0)*pow(1.0-c_n,1.0)
			                          * con_i*(a_n)
			                          * con_j*(1.0/2.0*(3.0*pow(b_n,2.0)-1.0))
			                          * con_k*(1.0));

			set_scaling_basis_pyr(2,0,0,c_n,&con_i,&con_j,&con_k,&con_c);
			*grad_phi_data[0]++ = con*( 2.0*pow(1.0-c_n,1.0)
			                          * con_i*(1.0/2.0*(3.0*2.0*a_n))
			                          * con_j*(1.0)
			                          * con_k*(1.0));
			*grad_phi_data[1]++ = 0.0;
			*grad_phi_data[2]++ = con*( sqrt(2.0)*a_n*pow(1.0-c_n,1.0)
			                          * con_i*(1.0/2.0*(3.0*2.0*a_n))
			                          * con_j*(1.0)
			                          * con_k*(1.0)
			                          - 2.0*sqrt(2.0)*pow(1.0-c_n,1.0)
			                          * con_i*(1.0/2.0*(3.0*pow(a_n,2.0)-1.0))
			                          * con_j*(1.0)
			                          * con_k*(1.0));

			set_scaling_basis_pyr(2,1,0,c_n,&con_i,&con_j,&con_k,&con_c);
			*grad_phi_data[0]++ = con*( 2.0*pow(1.0-c_n,1.0)
			                          * con_i*(1.0/2.0*(3.0*2.0*a_n))
			                          * con_j*(b_n)
			                          * con_k*(1.0));
			*grad_phi_data[1]++ = con*( 2.0*pow(1.0-c_n,1.0)
			                          * con_i*(1.0/2.0*(3.0*pow(a_n,2.0)-1.0))
			                          * con_j*(1.0)
			                          * con_k*(1.0));
			*grad_phi_data[2]++ = con*( sqrt(2.0)*a_n*pow(1.0-c_n,1.0)
			                          * con_i*(1.0/2.0*(3.0*2.0*a_n))
			                          * con_j*(b_n)
			                          * con_k*(1.0)
			                          + sqrt(2.0)*b_n*pow(1.0-c_n,1.0)
			                          * con_i*(1.0/2.0*(3.0*pow(a_n,2.0)-1.0))
			                          * con_j*(1.0)
			                          * con_k*(1.0)
			                          - 2.0*sqrt(2.0)*pow(1.0-c_n,1.0)
			                          * con_i*(1.0/2.0*(3.0*pow(a_n,2.0)-1.0))
			                          * con_j*(b_n)
			                          * con_k*(1.0));

			set_scaling_basis_pyr(2,2,0,c_n,&con_i,&con_j,&con_k,&con_c);
			*grad_phi_data[0]++ = con*( 2.0*pow(1.0-c_n,1.0)
			                          * con_i*(1.0/2.0*(3.0*2.0*a_n))
			                          * con_j*(1.0/2.0*(3.0*pow(b_n,2.0)-1.0))
			                          * con_k*(1.0));
			*grad_phi_data[1]++ = con*( 2.0*pow(1.0-c_n,1.0)
			                          * con_i*(1.0/2.0*(3.0*pow(a_n,2.0)-1.0))
			                          * con_j*(1.0/2.0*(3.0*2.0*b_n))
			                          * con_k*(1.0));
			*grad_phi_data[2]++ = con*( sqrt(2.0)*a_n*pow(1.0-c_n,1.0)
			                          * con_i*(1.0/2.0*(3.0*2.0*a_n))
			                          * con_j*(1.0/2.0*(3.0*pow(b_n,2.0)-1.0))
			                          * con_k*(1.0)
			                          + sqrt(2.0)*b_n*pow(1.0-c_n,1.0)
			                          * con_i*(1.0/2.0*(3.0*pow(a_n,2.0)-1.0))
			                          * con_j*(1.0/2.0*(3.0*2.0*b_n))
			                          * con_k*(1.0)
			                          - 2.0*sqrt(2.0)*pow(1.0-c_n,1.0)
			                          * con_i*(1.0/2.0*(3.0*pow(a_n,2.0)-1.0))
			                          * con_j*(1.0/2.0*(3.0*pow(b_n,2.0)-1.0))
			                          * con_k*(1.0));
		} else {
			EXIT_UNSUPPORTED;
		}
	}
	destructor_const_Matrix_d(abc);

	return (const struct const_Multiarray_Matrix_d*) grad_phi_rst;
}

// Tensor-Product Bezier ******************************************************************************************** //

const struct const_Matrix_d* constructor_basis_tp_bezier_def (const int p_b, const struct const_Matrix_d*const rst)
{
	assert(rst->layout == 'C');

	const ptrdiff_t d   = rst->ext_1,
	                n_n = rst->ext_0,
	                n_b = compute_n_basis(d,p_b,ST_TP);

	struct Matrix_d* phi_rst = constructor_empty_Matrix_d('R',n_n,n_b); // returned
	double* phi_data = phi_rst->data;

	const double*const r = get_col_const_Matrix_d(0,rst),
	            *const s = ( d > 1 ? get_col_const_Matrix_d(1,rst) : NULL),
	            *const t = ( d > 2 ? get_col_const_Matrix_d(2,rst) : NULL);

	for (ptrdiff_t i = 0; i < n_n; ++i) {
		const double r_i = r[i],
		             s_i = ( d > 1 ? s[i] : 0.0),
		             t_i = ( d > 2 ? t[i] : 0.0);

		const double b0r = (1.0-r_i)/2.0,
		             b1r = (1.0+r_i)/2.0,
		             b0s = ( d > 1 ? (1.0-s_i)/2.0 : 0.0 ),
		             b1s = ( d > 1 ? (1.0+s_i)/2.0 : 0.0 ),
		             b0t = ( d > 2 ? (1.0-t_i)/2.0 : 0.0 ),
		             b1t = ( d > 2 ? (1.0+t_i)/2.0 : 0.0 );

		if (d == 1 && p_b == 3) {
			*phi_data++ =     pow(b0r,3.0)*pow(b1r,0.0);
			*phi_data++ = 3.0*pow(b0r,2.0)*pow(b1r,1.0);
			*phi_data++ = 3.0*pow(b0r,1.0)*pow(b1r,2.0);
			*phi_data++ =     pow(b0r,0.0)*pow(b1r,3.0);
		} else if (d == 2 && p_b == 2) {
			*phi_data++ =     pow(b0r,2.0)*pow(b1r,0.0) *     pow(b0s,2.0)*pow(b1s,0.0);
			*phi_data++ = 2.0*pow(b0r,1.0)*pow(b1r,1.0) *     pow(b0s,2.0)*pow(b1s,0.0);
			*phi_data++ =     pow(b0r,0.0)*pow(b1r,2.0) *     pow(b0s,2.0)*pow(b1s,0.0);
			*phi_data++ =     pow(b0r,2.0)*pow(b1r,0.0) * 2.0*pow(b0s,1.0)*pow(b1s,1.0);
			*phi_data++ = 2.0*pow(b0r,1.0)*pow(b1r,1.0) * 2.0*pow(b0s,1.0)*pow(b1s,1.0);
			*phi_data++ =     pow(b0r,0.0)*pow(b1r,2.0) * 2.0*pow(b0s,1.0)*pow(b1s,1.0);
			*phi_data++ =     pow(b0r,2.0)*pow(b1r,0.0) *     pow(b0s,0.0)*pow(b1s,2.0);
			*phi_data++ = 2.0*pow(b0r,1.0)*pow(b1r,1.0) *     pow(b0s,0.0)*pow(b1s,2.0);
			*phi_data++ =     pow(b0r,0.0)*pow(b1r,2.0) *     pow(b0s,0.0)*pow(b1s,2.0);
		} else if (d == 3 && p_b == 2) {
			*phi_data++ =     pow(b0r,2.0)*pow(b1r,0.0) *     pow(b0s,2.0)*pow(b1s,0.0) *     pow(b0t,2.0)*pow(b1t,0.0);
			*phi_data++ = 2.0*pow(b0r,1.0)*pow(b1r,1.0) *     pow(b0s,2.0)*pow(b1s,0.0) *     pow(b0t,2.0)*pow(b1t,0.0);
			*phi_data++ =     pow(b0r,0.0)*pow(b1r,2.0) *     pow(b0s,2.0)*pow(b1s,0.0) *     pow(b0t,2.0)*pow(b1t,0.0);
			*phi_data++ =     pow(b0r,2.0)*pow(b1r,0.0) * 2.0*pow(b0s,1.0)*pow(b1s,1.0) *     pow(b0t,2.0)*pow(b1t,0.0);
			*phi_data++ = 2.0*pow(b0r,1.0)*pow(b1r,1.0) * 2.0*pow(b0s,1.0)*pow(b1s,1.0) *     pow(b0t,2.0)*pow(b1t,0.0);
			*phi_data++ =     pow(b0r,0.0)*pow(b1r,2.0) * 2.0*pow(b0s,1.0)*pow(b1s,1.0) *     pow(b0t,2.0)*pow(b1t,0.0);
			*phi_data++ =     pow(b0r,2.0)*pow(b1r,0.0) *     pow(b0s,0.0)*pow(b1s,2.0) *     pow(b0t,2.0)*pow(b1t,0.0);
			*phi_data++ = 2.0*pow(b0r,1.0)*pow(b1r,1.0) *     pow(b0s,0.0)*pow(b1s,2.0) *     pow(b0t,2.0)*pow(b1t,0.0);
			*phi_data++ =     pow(b0r,0.0)*pow(b1r,2.0) *     pow(b0s,0.0)*pow(b1s,2.0) *     pow(b0t,2.0)*pow(b1t,0.0);
			*phi_data++ =     pow(b0r,2.0)*pow(b1r,0.0) *     pow(b0s,2.0)*pow(b1s,0.0) * 2.0*pow(b0t,1.0)*pow(b1t,1.0);
			*phi_data++ = 2.0*pow(b0r,1.0)*pow(b1r,1.0) *     pow(b0s,2.0)*pow(b1s,0.0) * 2.0*pow(b0t,1.0)*pow(b1t,1.0);
			*phi_data++ =     pow(b0r,0.0)*pow(b1r,2.0) *     pow(b0s,2.0)*pow(b1s,0.0) * 2.0*pow(b0t,1.0)*pow(b1t,1.0);
			*phi_data++ =     pow(b0r,2.0)*pow(b1r,0.0) * 2.0*pow(b0s,1.0)*pow(b1s,1.0) * 2.0*pow(b0t,1.0)*pow(b1t,1.0);
			*phi_data++ = 2.0*pow(b0r,1.0)*pow(b1r,1.0) * 2.0*pow(b0s,1.0)*pow(b1s,1.0) * 2.0*pow(b0t,1.0)*pow(b1t,1.0);
			*phi_data++ =     pow(b0r,0.0)*pow(b1r,2.0) * 2.0*pow(b0s,1.0)*pow(b1s,1.0) * 2.0*pow(b0t,1.0)*pow(b1t,1.0);
			*phi_data++ =     pow(b0r,2.0)*pow(b1r,0.0) *     pow(b0s,0.0)*pow(b1s,2.0) * 2.0*pow(b0t,1.0)*pow(b1t,1.0);
			*phi_data++ = 2.0*pow(b0r,1.0)*pow(b1r,1.0) *     pow(b0s,0.0)*pow(b1s,2.0) * 2.0*pow(b0t,1.0)*pow(b1t,1.0);
			*phi_data++ =     pow(b0r,0.0)*pow(b1r,2.0) *     pow(b0s,0.0)*pow(b1s,2.0) * 2.0*pow(b0t,1.0)*pow(b1t,1.0);
			*phi_data++ =     pow(b0r,2.0)*pow(b1r,0.0) *     pow(b0s,2.0)*pow(b1s,0.0) *     pow(b0t,0.0)*pow(b1t,2.0);
			*phi_data++ = 2.0*pow(b0r,1.0)*pow(b1r,1.0) *     pow(b0s,2.0)*pow(b1s,0.0) *     pow(b0t,0.0)*pow(b1t,2.0);
			*phi_data++ =     pow(b0r,0.0)*pow(b1r,2.0) *     pow(b0s,2.0)*pow(b1s,0.0) *     pow(b0t,0.0)*pow(b1t,2.0);
			*phi_data++ =     pow(b0r,2.0)*pow(b1r,0.0) * 2.0*pow(b0s,1.0)*pow(b1s,1.0) *     pow(b0t,0.0)*pow(b1t,2.0);
			*phi_data++ = 2.0*pow(b0r,1.0)*pow(b1r,1.0) * 2.0*pow(b0s,1.0)*pow(b1s,1.0) *     pow(b0t,0.0)*pow(b1t,2.0);
			*phi_data++ =     pow(b0r,0.0)*pow(b1r,2.0) * 2.0*pow(b0s,1.0)*pow(b1s,1.0) *     pow(b0t,0.0)*pow(b1t,2.0);
			*phi_data++ =     pow(b0r,2.0)*pow(b1r,0.0) *     pow(b0s,0.0)*pow(b1s,2.0) *     pow(b0t,0.0)*pow(b1t,2.0);
			*phi_data++ = 2.0*pow(b0r,1.0)*pow(b1r,1.0) *     pow(b0s,0.0)*pow(b1s,2.0) *     pow(b0t,0.0)*pow(b1t,2.0);
			*phi_data++ =     pow(b0r,0.0)*pow(b1r,2.0) *     pow(b0s,0.0)*pow(b1s,2.0) *     pow(b0t,0.0)*pow(b1t,2.0);
		} else {
			EXIT_UNSUPPORTED;
		}
	}

	return (const struct const_Matrix_d*) phi_rst;
}

const struct const_Multiarray_Matrix_d* constructor_grad_basis_tp_bezier_def
	(const int p_b, const struct const_Matrix_d*const rst)
{
	assert(rst->layout == 'C');

	const ptrdiff_t d   = rst->ext_1,
	                n_n = rst->ext_0,
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

	for (ptrdiff_t i = 0; i < n_n; ++i) {
		const double r_i = r[i],
		             s_i = ( d > 1 ? s[i] : 0.0),
		             t_i = ( d > 2 ? t[i] : 0.0);

		const double b0r = (1.0-r_i)/2.0,
		             b1r = (1.0+r_i)/2.0,
		             b0s = ( d > 1 ? (1.0-s_i)/2.0 : 0.0 ),
		             b1s = ( d > 1 ? (1.0+s_i)/2.0 : 0.0 ),
		             b0t = ( d > 2 ? (1.0-t_i)/2.0 : 0.0 ),
		             b1t = ( d > 2 ? (1.0+t_i)/2.0 : 0.0 );
		if (d == 1 && p_b == 3) {
			*grad_phi_data[0]++ = 0.5*    (-3.0*pow(b0r,2.0)*pow(b1r,0.0) + 0.0);
			*grad_phi_data[0]++ = 0.5*3.0*(-2.0*pow(b0r,1.0)*pow(b1r,1.0) + 1.0*pow(b0r,2.0)*pow(b1r,0.0));
			*grad_phi_data[0]++ = 0.5*3.0*(-1.0*pow(b0r,0.0)*pow(b1r,2.0) + 2.0*pow(b0r,1.0)*pow(b1r,1.0));
			*grad_phi_data[0]++ = 0.5*    ( 0.0                           + 3.0*pow(b0r,0.0)*pow(b1r,2.0));
		} else if (d == 2 && p_b == 2) {
			*grad_phi_data[0]++ = 0.5*    (-2.0*pow(b0r,1.0)*pow(b1r,0.0) + 0.0)
			                    *               pow(b0s,2.0)*pow(b1s,0.0);
			*grad_phi_data[0]++ = 0.5*2.0*(-1.0*pow(b0r,0.0)*pow(b1r,1.0) + 1.0*pow(b0r,1.0)*pow(b1r,0.0))
			                    *               pow(b0s,2.0)*pow(b1s,0.0);
			*grad_phi_data[0]++ = 0.5*    ( 0.0                           + 2.0*pow(b0r,0.0)*pow(b1r,1.0))
			                    *               pow(b0s,2.0)*pow(b1s,0.0);
			*grad_phi_data[0]++ = 0.5*    (-2.0*pow(b0r,1.0)*pow(b1r,0.0) + 0.0)
			                    *           2.0*pow(b0s,1.0)*pow(b1s,1.0);
			*grad_phi_data[0]++ = 0.5*2.0*(-1.0*pow(b0r,0.0)*pow(b1r,1.0) + 1.0*pow(b0r,1.0)*pow(b1r,0.0))
			                    *           2.0*pow(b0s,1.0)*pow(b1s,1.0);
			*grad_phi_data[0]++ = 0.5*    ( 0.0                           + 2.0*pow(b0r,0.0)*pow(b1r,1.0))
			                    *           2.0*pow(b0s,1.0)*pow(b1s,1.0);
			*grad_phi_data[0]++ = 0.5*    (-2.0*pow(b0r,1.0)*pow(b1r,0.0) + 0.0)
			                    *               pow(b0s,0.0)*pow(b1s,2.0);
			*grad_phi_data[0]++ = 0.5*2.0*(-1.0*pow(b0r,0.0)*pow(b1r,1.0) + 1.0*pow(b0r,1.0)*pow(b1r,0.0))
			                    *               pow(b0s,0.0)*pow(b1s,2.0);
			*grad_phi_data[0]++ = 0.5*    ( 0.0                           + 2.0*pow(b0r,0.0)*pow(b1r,1.0))
			                    *               pow(b0s,0.0)*pow(b1s,2.0);

			*grad_phi_data[1]++ =               pow(b0r,2.0)*pow(b1r,0.0)
			                    * 0.5*    (-2.0*pow(b0s,1.0)*pow(b1s,0.0) + 0.0);
			*grad_phi_data[1]++ =           2.0*pow(b0r,1.0)*pow(b1r,1.0)
			                    * 0.5*    (-2.0*pow(b0s,1.0)*pow(b1s,0.0) + 0.0);
			*grad_phi_data[1]++ =               pow(b0r,0.0)*pow(b1r,2.0)
			                    * 0.5*    (-2.0*pow(b0s,1.0)*pow(b1s,0.0) + 0.0);
			*grad_phi_data[1]++ =               pow(b0r,2.0)*pow(b1r,0.0)
			                    * 0.5*2.0*(-1.0*pow(b0s,0.0)*pow(b1s,1.0) + 1.0*pow(b0s,1.0)*pow(b1s,0.0));
			*grad_phi_data[1]++ =           2.0*pow(b0r,1.0)*pow(b1r,1.0)
			                    * 0.5*2.0*(-1.0*pow(b0s,0.0)*pow(b1s,1.0) + 1.0*pow(b0s,1.0)*pow(b1s,0.0));
			*grad_phi_data[1]++ =               pow(b0r,0.0)*pow(b1r,2.0)
			                    * 0.5*2.0*(-1.0*pow(b0s,0.0)*pow(b1s,1.0) + 1.0*pow(b0s,1.0)*pow(b1s,0.0));
			*grad_phi_data[1]++ =               pow(b0r,2.0)*pow(b1r,0.0)
			                    * 0.5*    ( 0.0                           + 2.0*pow(b0s,0.0)*pow(b1s,1.0));
			*grad_phi_data[1]++ =           2.0*pow(b0r,1.0)*pow(b1r,1.0)
			                    * 0.5*    ( 0.0                           + 2.0*pow(b0s,0.0)*pow(b1s,1.0));
			*grad_phi_data[1]++ =               pow(b0r,0.0)*pow(b1r,2.0)
			                    * 0.5*    ( 0.0                           + 2.0*pow(b0s,0.0)*pow(b1s,1.0));
		} else if (d == 3 && p_b == 1) {
			*grad_phi_data[0]++ = 0.5*(-1.0*pow(b0r,0.0)*pow(b1r,0.0) + 0.0)
			                    *           pow(b0s,1.0)*pow(b1s,0.0)
			                    *           pow(b0t,1.0)*pow(b1t,0.0);
			*grad_phi_data[0]++ = 0.5*( 0.0                           + 1.0*pow(b0r,0.0)*pow(b1r,0.0))
			                    *           pow(b0s,1.0)*pow(b1s,0.0)
			                    *           pow(b0t,1.0)*pow(b1t,0.0);
			*grad_phi_data[0]++ = 0.5*(-1.0*pow(b0r,0.0)*pow(b1r,0.0) + 0.0)
			                    *           pow(b0s,0.0)*pow(b1s,1.0)
			                    *           pow(b0t,1.0)*pow(b1t,0.0);
			*grad_phi_data[0]++ = 0.5*( 0.0                           + 1.0*pow(b0r,0.0)*pow(b1r,0.0))
			                    *           pow(b0s,0.0)*pow(b1s,1.0)
			                    *           pow(b0t,1.0)*pow(b1t,0.0);
			*grad_phi_data[0]++ = 0.5*(-1.0*pow(b0r,0.0)*pow(b1r,0.0) + 0.0)
			                    *           pow(b0s,1.0)*pow(b1s,0.0)
			                    *           pow(b0t,0.0)*pow(b1t,1.0);
			*grad_phi_data[0]++ = 0.5*( 0.0                           + 1.0*pow(b0r,0.0)*pow(b1r,0.0))
			                    *           pow(b0s,1.0)*pow(b1s,0.0)
			                    *           pow(b0t,0.0)*pow(b1t,1.0);
			*grad_phi_data[0]++ = 0.5*(-1.0*pow(b0r,0.0)*pow(b1r,0.0) + 0.0)
			                    *           pow(b0s,0.0)*pow(b1s,1.0)
			                    *           pow(b0t,0.0)*pow(b1t,1.0);
			*grad_phi_data[0]++ = 0.5*( 0.0                           + 1.0*pow(b0r,0.0)*pow(b1r,0.0))
			                    *           pow(b0s,0.0)*pow(b1s,1.0)
			                    *           pow(b0t,0.0)*pow(b1t,1.0);

			*grad_phi_data[1]++ =           pow(b0r,1.0)*pow(b1r,0.0)
			                    * 0.5*(-1.0*pow(b0s,0.0)*pow(b1s,0.0) + 0.0)
			                    *           pow(b0t,1.0)*pow(b1t,0.0);
			*grad_phi_data[1]++ =           pow(b0r,0.0)*pow(b1r,1.0)
			                    * 0.5*(-1.0*pow(b0s,0.0)*pow(b1s,0.0) + 0.0)
			                    *           pow(b0t,1.0)*pow(b1t,0.0);
			*grad_phi_data[1]++ =           pow(b0r,1.0)*pow(b1r,0.0)
			                    * 0.5*( 0.0                           + 1.0*pow(b0s,0.0)*pow(b1s,0.0))
			                    *           pow(b0t,1.0)*pow(b1t,0.0);
			*grad_phi_data[1]++ =           pow(b0r,0.0)*pow(b1r,1.0)
			                    * 0.5*( 0.0                           + 1.0*pow(b0s,0.0)*pow(b1s,0.0))
			                    *           pow(b0t,1.0)*pow(b1t,0.0);
			*grad_phi_data[1]++ =           pow(b0r,1.0)*pow(b1r,0.0)
			                    * 0.5*(-1.0*pow(b0s,0.0)*pow(b1s,0.0) + 0.0)
			                    *           pow(b0t,0.0)*pow(b1t,1.0);
			*grad_phi_data[1]++ =           pow(b0r,0.0)*pow(b1r,1.0)
			                    * 0.5*(-1.0*pow(b0s,0.0)*pow(b1s,0.0) + 0.0)
			                    *           pow(b0t,0.0)*pow(b1t,1.0);
			*grad_phi_data[1]++ =           pow(b0r,1.0)*pow(b1r,0.0)
			                    * 0.5*( 0.0                           + 1.0*pow(b0s,0.0)*pow(b1s,0.0))
			                    *           pow(b0t,0.0)*pow(b1t,1.0);
			*grad_phi_data[1]++ =           pow(b0r,0.0)*pow(b1r,1.0)
			                    * 0.5*( 0.0                           + 1.0*pow(b0s,0.0)*pow(b1s,0.0))
			                    *           pow(b0t,0.0)*pow(b1t,1.0);

			*grad_phi_data[2]++ =           pow(b0r,1.0)*pow(b1r,0.0)
			                    *           pow(b0s,1.0)*pow(b1s,0.0)
			                    * 0.5*(-1.0*pow(b0t,0.0)*pow(b1t,0.0) + 0.0);
			*grad_phi_data[2]++ =           pow(b0r,0.0)*pow(b1r,1.0)
			                    *           pow(b0s,1.0)*pow(b1s,0.0)
			                    * 0.5*(-1.0*pow(b0t,0.0)*pow(b1t,0.0) + 0.0);
			*grad_phi_data[2]++ =           pow(b0r,1.0)*pow(b1r,0.0)
			                    *           pow(b0s,0.0)*pow(b1s,1.0)
			                    * 0.5*(-1.0*pow(b0t,0.0)*pow(b1t,0.0) + 0.0);
			*grad_phi_data[2]++ =           pow(b0r,0.0)*pow(b1r,1.0)
			                    *           pow(b0s,0.0)*pow(b1s,1.0)
			                    * 0.5*(-1.0*pow(b0t,0.0)*pow(b1t,0.0) + 0.0);
			*grad_phi_data[2]++ =           pow(b0r,1.0)*pow(b1r,0.0)
			                    *           pow(b0s,1.0)*pow(b1s,0.0)
			                    * 0.5*( 0.0                           + 1.0*pow(b0t,0.0)*pow(b1t,0.0));
			*grad_phi_data[2]++ =           pow(b0r,0.0)*pow(b1r,1.0)
			                    *           pow(b0s,1.0)*pow(b1s,0.0)
			                    * 0.5*( 0.0                           + 1.0*pow(b0t,0.0)*pow(b1t,0.0));
			*grad_phi_data[2]++ =           pow(b0r,1.0)*pow(b1r,0.0)
			                    *           pow(b0s,0.0)*pow(b1s,1.0)
			                    * 0.5*( 0.0                           + 1.0*pow(b0t,0.0)*pow(b1t,0.0));
			*grad_phi_data[2]++ =           pow(b0r,0.0)*pow(b1r,1.0)
			                    *           pow(b0s,0.0)*pow(b1s,1.0)
			                    * 0.5*( 0.0                           + 1.0*pow(b0t,0.0)*pow(b1t,0.0));
		} else {
			EXIT_UNSUPPORTED;
		}
	}

	return (const struct const_Multiarray_Matrix_d*) grad_phi_rst;
}

// Additional functions ********************************************************************************************* //

const struct const_Matrix_d* constructor_mass_orthonormal_def (const int d, const int p_b, const int super_type)
{
	const ptrdiff_t n_b = compute_n_basis(d,p_b,super_type);
	return constructor_identity_const_Matrix_d('R',n_b);
}

const struct const_Matrix_d* constructor_mass_orthonormal (const int d, const int p_b, const int super_type)
{
	int node_type  = 0;
	int node_order = 0;
	constructor_Nodes_fptr nodes_fun = NULL;
	constructor_basis_fptr basis_fun = NULL;
	if (super_type == ST_TP) {
		node_type  = NODES_GL;
		node_order = p_b;
		nodes_fun   = constructor_const_Nodes_tp;
		basis_fun = constructor_basis_tp_orthonormal;
	} else if (super_type == ST_SI) {
		node_type  = NODES_WV;
		node_order = 2*p_b;
		nodes_fun   = constructor_const_Nodes_si;
		basis_fun = constructor_basis_si_orthonormal;
	} else if (super_type == ST_PYR) {
		const int node_opt = 0;
		const int node_type_opt[]  = { NODES_GJW, NODES_GLW, NODES_GLLW, };
		const int node_order_opt[] = { p_b,       p_b+1,     p_b+2,      };

		node_type  = node_type_opt[node_opt];
		node_order = node_order_opt[node_opt];
		nodes_fun   = constructor_const_Nodes_pyr;
		basis_fun = constructor_basis_pyr_orthonormal;
	} else {
		EXIT_UNSUPPORTED;
	}

	const struct const_Nodes* nodes = nodes_fun(d,node_order,node_type); // destructed

	const struct const_Matrix_d*const phi   = basis_fun(p_b,nodes->rst),            // destructed
	                           *const phi_w = constructor_copy_const_Matrix_d(phi); // destructed

	transpose_const_Matrix_d(phi_w,false);
	scale_const_Matrix_by_Vector_d('R',1.0,phi_w,nodes->w,false);

	const struct const_Matrix_d* mass = constructor_mm_const_Matrix_d('N','N',1.0,phi_w,phi,'R'); // returned

	destructor_const_Matrix_d(phi_w);
	destructor_const_Matrix_d(phi);
	destructor_const_Nodes(nodes);

	return mass;
}

const struct const_Vector_d* constructor_part_unity_def (const ptrdiff_t n_n)
{
	struct Vector_d* dest = constructor_empty_Vector_d(n_n); // returned
	set_to_value_Vector_d(dest,1.0);
	return (const struct const_Vector_d*) dest;
}

const struct const_Vector_d* constructor_part_unity (const struct const_Matrix_d* phi_rst)
{
	assert(phi_rst->layout == 'R');
	return constructor_sum_const_Vector_d_const_Matrix_d('R',phi_rst); // returned
}


const struct const_Multiarray_d* constructor_grad_vals_computation_def
	(const int d, const int p_b, const char*const basis_name)
{
	int node_type  = 0;
	int node_order = 0;
	constructor_Nodes_fptr nodes_fun = NULL;
	if (strcmp(basis_name,"tp_ortho") == 0) {
		node_type  = NODES_GLL;
		node_order = p_b;
		nodes_fun   = constructor_const_Nodes_tp;
	} else if (strcmp(basis_name,"si_ortho") == 0) {
		node_type  = NODES_AO;
		node_order = p_b;
		nodes_fun   = constructor_const_Nodes_si;
	} else if (strcmp(basis_name,"pyr_ortho") == 0) {
		node_type  = NODES_GLL;
		node_order = p_b;
		nodes_fun   = constructor_const_Nodes_pyr;
	} else if (strcmp(basis_name,"tp_bezier") == 0) {
		node_type  = NODES_GLL;
		node_order = p_b;
		nodes_fun   = constructor_const_Nodes_tp;
	} else {
		EXIT_UNSUPPORTED;
	}

	const struct const_Nodes* nodes = nodes_fun(d,node_order,node_type); // destructed
	const struct const_Multiarray_d* grad_f_vc = constructor_grad_f_rst(nodes->rst); // returned

	destructor_const_Nodes(nodes);
	return grad_f_vc;
}

const struct const_Multiarray_d* constructor_grad_vals_computation
	(const int d, const int p_b, const char*const basis_name)
{
	int node_type  = 0;
	int node_order = 0;
	constructor_Nodes_fptr nodes_fun           = NULL;
	constructor_basis_fptr basis_fun           = NULL;
	constructor_grad_basis_fptr grad_basis_fun = NULL;
	if (strcmp(basis_name,"tp_ortho") == 0) {
		node_type  = NODES_GLL;
		node_order = p_b;
		nodes_fun   = constructor_const_Nodes_tp;
		basis_fun      = constructor_basis_tp_orthonormal;
		grad_basis_fun = constructor_grad_basis_tp_orthonormal;
	} else if (strcmp(basis_name,"si_ortho") == 0) {
		node_type  = NODES_AO;
		node_order = p_b;
		nodes_fun   = constructor_const_Nodes_si;
		basis_fun      = constructor_basis_si_orthonormal;
		grad_basis_fun = constructor_grad_basis_si_orthonormal;
	} else if (strcmp(basis_name,"pyr_ortho") == 0) {
		node_type  = NODES_GLL;
		node_order = p_b;
		nodes_fun   = constructor_const_Nodes_pyr;
		basis_fun      = constructor_basis_pyr_orthonormal;
		grad_basis_fun = constructor_grad_basis_pyr_orthonormal;
	} else if (strcmp(basis_name,"tp_bezier") == 0) {
		node_type  = NODES_GLL;
		node_order = p_b;
		nodes_fun   = constructor_const_Nodes_tp;
		basis_fun      = constructor_basis_tp_bezier;
		grad_basis_fun = constructor_grad_basis_tp_bezier;
	} else {
		EXIT_UNSUPPORTED;
	}

	const struct const_Nodes* nodes = nodes_fun(d,node_order,node_type); // destructed
	const struct const_Multiarray_d* f_vc_MA = constructor_f_rst(nodes->rst); // destructed

	const struct const_Vector_d* f_vc =
		constructor_set_const_Vector_d_Multiarray_d(f_vc_MA,(ptrdiff_t[]){0}); // destructed

	const struct const_Matrix_d*const phi = basis_fun(p_b,nodes->rst); // destructed
	const struct const_Vector_d* f_vc_coef = constructor_sgesv_const_Vector_d(phi,f_vc); // destructed
	destructor_const_Matrix_d(phi);
	destructor_const_Vector_d(f_vc);
	destructor_const_Multiarray_d(f_vc_MA);

	const struct const_Multiarray_Matrix_d*const grad_phi = grad_basis_fun(p_b,nodes->rst); // destructed
	const struct const_Multiarray_d* grad_f =
		constructor_MaM1_V_const_Multiarray_d('C','N',1.0,0.0,grad_phi,f_vc_coef); // returned
	destructor_const_Multiarray_Matrix_d(grad_phi);
	destructor_const_Vector_d(f_vc_coef);
	destructor_const_Nodes(nodes);

	return (const struct const_Multiarray_d*) grad_f;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void set_scaling_basis_tri
	(const int i, const int j, const double b, double* con_i, double* con_j, double* con_b)
{
	*con_i = sqrt((2.0*i+1.0)/2.0);
	*con_j = sqrt((i+j+1.0)/pow(2.0,2.0*i+1.0));
	*con_b = pow(1.0-b,i);
}

static void set_scaling_basis_tet
	(const int i, const int j, const int k, const double b, const double c,
         double* con_i, double* con_j, double* con_k, double* con_b, double* con_c)
{
	set_scaling_basis_tri(i,j,b,con_i,con_j,con_b);
	*con_k = sqrt((2.0*(i+j+k)+3.0)/pow(2.0,2.0*(i+j)+3.0));
	*con_c = pow(1.0-c,i+j);
}

static void set_scaling_basis_pyr
	(const int i, const int j, const int k, const double c,
	 double* con_i, double* con_j, double* con_k, double* con_c)
{
	const double mu_ij = GSL_MAX(i,j);

	*con_i = sqrt((2.0*i+1.0)/2.0);
	*con_j = sqrt((2.0*j+1.0)/2.0);
	*con_k = sqrt((2.0*(mu_ij+k)+3.0)/pow(2.0,2.0*mu_ij+3.0));
	*con_c = pow(1.0-c,mu_ij);
}

static double poly_rst (const double r, const double s, const double t, const int derivative_index)
{
	assert((0 <= derivative_index) && (derivative_index <= 3));

	const int d     = 3;
	const int n_der = 2;

	const double rst[]  = { r, s, t, };
	const int len[]     = { 3, 3, 3, };
	const double c[][DMAX] = { { 1.0, 2.0, 5.0 },
	                           { 1.0, 3.0, 6.0 },
	                           { 1.0, 4.0, 7.0 }, };

	double dc[n_der];

	double result = 1;
	for (int dim = 0; dim < d; ++dim) {
		gsl_poly_eval_derivs(c[dim],len[dim],rst[dim],dc,n_der);
		const int ind_dc = ( dim+1 != derivative_index ? 0 : 1);
		result *= dc[ind_dc];
	}

	return result;
}

static const struct const_Multiarray_d* constructor_f_rst (const struct const_Matrix_d* rst)
{
	assert(rst->layout == 'C');
	const ptrdiff_t n_n = rst->ext_0,
	                d   = rst->ext_1;

	const double*const r_ptr = get_col_const_Matrix_d(0,rst),
	            *const s_ptr = ( d > 1 ? get_col_const_Matrix_d(1,rst) : NULL),
	            *const t_ptr = ( d > 2 ? get_col_const_Matrix_d(2,rst) : NULL);

	struct Matrix_d* f_vc = constructor_empty_Matrix_d('C',n_n,1); // destructed

	double* data = f_vc->data;
	for (int n = 0; n < n_n; ++n) {
		const double r = r_ptr[n],
		             s = ( d > 1 ? s_ptr[n] : 0.0),
		             t = ( d > 2 ? t_ptr[n] : 0.0);

		*data++ = poly_rst(r,s,t,0);
	}

	const struct const_Multiarray_d* dest =
		constructor_move_const_Multiarray_d_Matrix_d((const struct const_Matrix_d*)f_vc);
	destructor_Matrix_d(f_vc);

	return dest;
}

static const struct const_Multiarray_d* constructor_grad_f_rst (const struct const_Matrix_d* rst)
{
	assert(rst->layout == 'C');
	const ptrdiff_t n_n = rst->ext_0,
	                d   = rst->ext_1;

	const double*const r_ptr = get_col_const_Matrix_d(0,rst),
	            *const s_ptr = ( d > 1 ? get_col_const_Matrix_d(1,rst) : NULL),
	            *const t_ptr = ( d > 2 ? get_col_const_Matrix_d(2,rst) : NULL);

	struct Matrix_d* grad_f_vc = constructor_empty_Matrix_d('C',n_n,d); // destructed

	double* data = grad_f_vc->data;
	for (int dim = 0; dim < d; ++dim) {
	for (int n = 0; n < n_n; ++n) {
		const double r = r_ptr[n],
		             s = ( d > 1 ? s_ptr[n] : 0.0),
		             t = ( d > 2 ? t_ptr[n] : 0.0);

		*data++ = poly_rst(r,s,t,dim+1);
	}}

	const struct const_Multiarray_d* dest =
		constructor_move_const_Multiarray_d_Matrix_d((const struct const_Matrix_d*)grad_f_vc); // returned
	destructor_Matrix_d(grad_f_vc);

	return dest;
}
