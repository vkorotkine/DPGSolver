// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/** \file
 */

#include "test_support_bases.h"

#include <stdlib.h>
#include <math.h>
#include "gsl/gsl_math.h"

#include "test_support_matrix.h"

#include "macros.h"
#include "definitions_cubature.h"
#include "definitions_elements.h"

#include "matrix.h"
#include "vector.h"

#include "bases.h"
#include "cubature.h"

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

// Interface functions ********************************************************************************************** //
// Constructor functions ******************************************************************************************** //
// Tensor-Product Orthonormal *************************************************************************************** //

const struct const_Matrix_d* constructor_basis_tp_orthonormal_def (const int p_b, const struct const_Matrix_d*const rst)
{
	if (rst->layout != 'C')
		EXIT_UNSUPPORTED;

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

// Simplex Orthonormal ********************************************************************************************** //

const struct const_Matrix_d* constructor_basis_si_orthonormal_def (const int p_b, const struct const_Matrix_d*const rst)
{
	if (rst->layout != 'C')
		EXIT_UNSUPPORTED;

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

// Pyramid Orthonormal ********************************************************************************************** //

const struct const_Matrix_d* constructor_basis_pyr_orthonormal_def
	(const int p_b, const struct const_Matrix_d*const rst)
{
	if (rst->layout != 'C')
		EXIT_UNSUPPORTED;

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

// Tensor-Product Bezier ******************************************************************************************** //

const struct const_Matrix_d* constructor_basis_tp_bezier_def (const int p_b, const struct const_Matrix_d*const rst)
{
	if (rst->layout != 'C')
		EXIT_UNSUPPORTED;

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

// Orthogonality functions ****************************************************************************************** //

const struct const_Matrix_d* constructor_mass_orthonormal_def (const int d, const int p_b, const int super_type)
{
	const ptrdiff_t n_b = compute_n_basis(d,p_b,super_type);
	return constructor_identity_const_Matrix_d('R',n_b);
}

const struct const_Matrix_d* constructor_mass_orthonormal (const int d, const int p_b, const int super_type)
{
	int cub_type  = 0;
	int cub_order = 0;
	cubature_fptr cub_fun   = NULL;
	basis_fptr    basis_fun = NULL;
	if (super_type == ST_TP) {
		cub_type  = CUB_GL;
		cub_order = p_b;
		cub_fun   = constructor_const_Cubature_tp;
		basis_fun = constructor_basis_tp_orthonormal;
	} else if (super_type == ST_SI) {
		cub_type  = CUB_WV;
		cub_order = 2*p_b;
		cub_fun   = constructor_const_Cubature_si;
		basis_fun = constructor_basis_si_orthonormal;
	} else if (super_type == ST_PYR) {
		const int cub_opt = 0;
		const int cub_type_opt[]  = { CUB_GJW, CUB_GLW, CUB_GLLW, };
		const int cub_order_opt[] = { p_b,     p_b+1,   p_b+2,    };

		cub_type  = cub_type_opt[cub_opt];
		cub_order = cub_order_opt[cub_opt];
		cub_fun   = constructor_const_Cubature_pyr;
		basis_fun = constructor_basis_pyr_orthonormal;
	} else {
		EXIT_UNSUPPORTED;
	}

	const struct const_Cubature* cub = cub_fun(d,cub_order,cub_type); // destructed

	const struct const_Matrix_d*const phi   = basis_fun(p_b,cub->rst),              // destructed
	                           *const phi_w = constructor_copy_const_Matrix_d(phi); // destructed

	transpose_const_Matrix_d(phi_w,false);
	scale_const_Matrix_by_Vector_d('R',1.0,phi_w,cub->w);

	const struct const_Matrix_d* mass = constructor_mm_const_Matrix_d('N','N',1.0,0.0,phi_w,phi,'R'); // returned

	destructor_const_Matrix_d(phi_w);
	destructor_const_Matrix_d(phi);
	destructor_const_Cubature(cub);

	return mass;
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
