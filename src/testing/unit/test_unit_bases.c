// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/** \file
 */

#include "test_unit_bases.h"

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

#include "test_base.h"
#include "test_support.h"
#include "test_support_matrix.h"

#include "macros.h"
#include "definitions_alloc.h"
#include "definitions_cubature.h"
#include "definitions_tol.h"
#include "definitions_elements.h"

#include "matrix.h"

#include "bases.h"
#include "cubature.h"

// Static function declarations ************************************************************************************* //

/// \brief Provides unit tests for the orthonormal tensor-product basis functions.
static void test_unit_basis_tensor_product_orthonormal
	(struct Test_Info*const test_info ///< \ref Test_Info.
	);

/// \brief Provides unit tests for the orthonormal simplex basis functions.
static void test_unit_basis_simplex_orthonormal
	(struct Test_Info*const test_info ///< \ref Test_Info.
	);

// Interface functions ********************************************************************************************** //

void test_unit_bases (struct Test_Info*const test_info)
{
	test_unit_basis_tensor_product_orthonormal(test_info);
//	test_unit_basis_tensor_product_bezier(test_info);
	test_unit_basis_simplex_orthonormal(test_info);
//	test_unit_basis_simplex_bezier(test_info);
//	test_unit_basis_pyramid(test_info);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
// Tensor-Product Orthonormal *************************************************************************************** //

/// Container for tensor-product orthonormal basis data to be tested for comparison with expected values.
struct Basis_Data_TP_Ortho {
	const struct const_Matrix_d* phi13, ///< The 1d basis functions of order 3.
	                           * phi22, ///< The 2d basis functions of order 2.
	                           * phi31; ///< The 3d basis functions of order 1.
};

/** \brief Constructor for \ref Basis_Data_TP_Ortho.
 *  \return Standard. */
static struct Basis_Data_TP_Ortho* constructor_Basis_Data_TP_Ortho
	(const char eval_type ///< Method to use to obtain the data. Options: 'a'nalytical, 'c'ompute.
	);

/// \brief Destructor for \ref Basis_Data_TP_Ortho.
static void destructor_Basis_Data_TP_Ortho
	(struct Basis_Data_TP_Ortho* b_data ///< Standard.
	);

static void test_unit_basis_tensor_product_orthonormal (struct Test_Info*const test_info)
{
	char* test_name = "Bases - basis_tp_orthonormal";
	bool pass = true;

	struct Basis_Data_TP_Ortho* b_data_a = constructor_Basis_Data_TP_Ortho('a'), // destructed
	                          * b_data_c = constructor_Basis_Data_TP_Ortho('c'); // destructed

	double tol[]       = { EPS, 2*EPS, EPS, };
	bool differences[] =
		{ diff_const_Matrix_d(b_data_a->phi13,b_data_c->phi13,tol[0]),
		  diff_const_Matrix_d(b_data_a->phi22,b_data_c->phi22,tol[1]),
		  diff_const_Matrix_d(b_data_a->phi31,b_data_c->phi31,tol[2]),
		};

	const bool diff = check_diff(sizeof(differences)/sizeof(*differences),differences);

	if (diff) {
		pass = false;

		if (differences[0]) print_diff_const_Matrix_d(b_data_a->phi13,b_data_c->phi13,tol[0]);
		if (differences[1]) print_diff_const_Matrix_d(b_data_a->phi22,b_data_c->phi22,tol[1]);
		if (differences[2]) print_diff_const_Matrix_d(b_data_a->phi31,b_data_c->phi31,tol[2]);
	}

	destructor_Basis_Data_TP_Ortho(b_data_a);
	destructor_Basis_Data_TP_Ortho(b_data_c);

	test_increment_and_print(test_info,pass,test_name);
}

// Tensor-Product Bezier ******************************************************************************************** //
// Simplex Orthonormal ********************************************************************************************** //

/// Container for tensor-product orthonormal basis data to be tested for comparison with expected values.
struct Basis_Data_SI_Ortho {
	const struct const_Matrix_d* phi22, ///< The 2d basis functions of order 2.
	                           * phi23, ///< The 2d basis functions of order 3.
	                           * phi32; ///< The 3d basis functions of order 2.
};

/** \brief Constructor for \ref Basis_Data_SI_Ortho.
 *  \return Standard. */
static struct Basis_Data_SI_Ortho* constructor_Basis_Data_SI_Ortho
	(const char eval_type ///< Method to use to obtain the data. Options: 'a'nalytical, 'c'ompute.
	);

/// \brief Destructor for \ref Basis_Data_SI_Ortho.
static void destructor_Basis_Data_SI_Ortho
	(struct Basis_Data_SI_Ortho* b_data ///< Standard.
	);

static void test_unit_basis_simplex_orthonormal (struct Test_Info*const test_info)
{
	char* test_name = "Bases - basis_si_orthonormal";
	bool pass = true;

	struct Basis_Data_SI_Ortho* b_data_a = constructor_Basis_Data_SI_Ortho('a'), // destructed
	                          * b_data_c = constructor_Basis_Data_SI_Ortho('c'); // destructed

	double tol[]       = { 10*EPS, 20*EPS, 10*EPS, };
	bool differences[] =
		{ diff_const_Matrix_d(b_data_a->phi22,b_data_c->phi22,tol[0]),
		  diff_const_Matrix_d(b_data_a->phi23,b_data_c->phi23,tol[1]),
		  diff_const_Matrix_d(b_data_a->phi32,b_data_c->phi32,tol[2]),
		};

	const bool diff = check_diff(sizeof(differences)/sizeof(*differences),differences);

	if (diff) {
		pass = false;

		if (differences[0]) print_diff_const_Matrix_d(b_data_a->phi22,b_data_c->phi22,tol[0]);
		if (differences[1]) print_diff_const_Matrix_d(b_data_a->phi23,b_data_c->phi23,tol[1]);
		if (differences[2]) print_diff_const_Matrix_d(b_data_a->phi32,b_data_c->phi32,tol[2]);
	}

	destructor_Basis_Data_SI_Ortho(b_data_a);
	destructor_Basis_Data_SI_Ortho(b_data_c);

	test_increment_and_print(test_info,pass,test_name);
}

// Pyramid Orthonormal ********************************************************************************************** //


// Level 1 ********************************************************************************************************** //
// Tensor-Product Orthonormal *************************************************************************************** //

/** \brief Constructor for a polynomial operator for the tensor-product orthonomal basis from the basis function
 *         definitions.
 *  \return Standard. */
static const struct const_Matrix_d* constructor_basis_tp_orthonormal_def
	(const int p_b,                        ///< The order of the basis.
	 const struct const_Matrix_d*const rst ///< The nodes at which the basis functions are evaluated.
	);

static struct Basis_Data_TP_Ortho* constructor_Basis_Data_TP_Ortho (const char eval_type)
{
	struct Basis_Data_TP_Ortho* b_data = calloc(1,sizeof *b_data); // returned

	const struct const_Cubature* d1_p4_GLL = constructor_const_Cubature_tp(1,4,CUB_GLL), // destructed
	                           * d2_p4_GLL = constructor_const_Cubature_tp(2,4,CUB_GLL), // destructed
	                           * d3_p4_GLL = constructor_const_Cubature_tp(3,4,CUB_GLL); // destructed
	if (eval_type == 'a') {
		b_data->phi13 = constructor_basis_tp_orthonormal_def(3,d1_p4_GLL->rst); // keep
		b_data->phi22 = constructor_basis_tp_orthonormal_def(2,d2_p4_GLL->rst); // keep
		b_data->phi31 = constructor_basis_tp_orthonormal_def(1,d3_p4_GLL->rst); // keep
	} else if (eval_type == 'c') {
		b_data->phi13 = constructor_basis_tp_orthonormal(3,d1_p4_GLL->rst); // keep
		b_data->phi22 = constructor_basis_tp_orthonormal(2,d2_p4_GLL->rst); // keep
		b_data->phi31 = constructor_basis_tp_orthonormal(1,d3_p4_GLL->rst); // keep
	} else {
		EXIT_UNSUPPORTED;
	}
	destructor_const_Cubature(d1_p4_GLL);
	destructor_const_Cubature(d2_p4_GLL);
	destructor_const_Cubature(d3_p4_GLL);

	return b_data;
}

static void destructor_Basis_Data_TP_Ortho (struct Basis_Data_TP_Ortho* b_data)
{
	destructor_const_Matrix_d(b_data->phi13);
	destructor_const_Matrix_d(b_data->phi22);
	destructor_const_Matrix_d(b_data->phi31);
	free(b_data);
}

// Simplex Orthonormal ********************************************************************************************** //

/** \brief Constructor for a polynomial operator for the simplex orthonomal basis from the basis function definitions.
 *  \return Standard. */
static const struct const_Matrix_d* constructor_basis_si_orthonormal_def
	(const int p_b,                        ///< The order of the basis.
	 const struct const_Matrix_d*const rst ///< The nodes at which the basis functions are evaluated.
	);

static struct Basis_Data_SI_Ortho* constructor_Basis_Data_SI_Ortho (const char eval_type)
{
	struct Basis_Data_SI_Ortho* b_data = calloc(1,sizeof *b_data); // returned

	const struct const_Cubature* d2_p4_AO = constructor_const_Cubature_si(2,4,CUB_AO), // destructed
	                           * d3_p4_AO = constructor_const_Cubature_si(3,4,CUB_AO); // destructed
	if (eval_type == 'a') {
		b_data->phi22 = constructor_basis_si_orthonormal_def(2,d2_p4_AO->rst); // keep
		b_data->phi23 = constructor_basis_si_orthonormal_def(3,d2_p4_AO->rst); // keep
		b_data->phi32 = constructor_basis_si_orthonormal_def(2,d3_p4_AO->rst); // keep
	} else if (eval_type == 'c') {
		b_data->phi22 = constructor_basis_si_orthonormal(2,d2_p4_AO->rst); // keep
		b_data->phi23 = constructor_basis_si_orthonormal(3,d2_p4_AO->rst); // keep
		b_data->phi32 = constructor_basis_si_orthonormal(2,d3_p4_AO->rst); // keep
	} else {
		EXIT_UNSUPPORTED;
	}
	destructor_const_Cubature(d2_p4_AO);
	destructor_const_Cubature(d3_p4_AO);

	return b_data;
}

static void destructor_Basis_Data_SI_Ortho (struct Basis_Data_SI_Ortho* b_data)
{
	destructor_const_Matrix_d(b_data->phi22);
	destructor_const_Matrix_d(b_data->phi23);
	destructor_const_Matrix_d(b_data->phi32);
	free(b_data);
}

// Level 2 ********************************************************************************************************** //
// Tensor-Product Orthonormal *************************************************************************************** //

static const struct const_Matrix_d* constructor_basis_tp_orthonormal_def
	(const int p_b, const struct const_Matrix_d*const rst)
{
	if (rst->layout != 'C')
		EXIT_UNSUPPORTED;

	const ptrdiff_t d   = rst->ext_1,
	                n_n = rst->ext_0,
	                n_b = pow(p_b+1,d);

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

/// \brief Set the scaling parameters for the tri orthonormal basis functions.
void set_scaling_basis_tri
	(const int i,    ///< Index of the `a` polynomial component.
	 const int j,    ///< Index of the `b` polynomial component.
	 const double b, ///< Value of `b` (from abc coordinates).
	 double *con_i,  ///< Constant multiplier for the `a` polynomial component.
	 double *con_j,  ///< Constant multiplier for the `b` polynomial component.
	 double *con_b   ///< Constant multiplier for the additional `b` component.
	);

/// \brief Set the scaling parameters for the tet orthonormal basis functions.
void set_scaling_basis_tet
	(const int i,    ///< Index of the `a` polynomial component.
	 const int j,    ///< Index of the `b` polynomial component.
	 const int k,    ///< Index of the `c` polynomial component.
	 const double b, ///< Value of `b` (from abc coordinates).
	 const double c, ///< Value of `c` (from abc coordinates).
         double *con_i,  ///< Constant multiplier for the `a` polynomial component.
	 double *con_j,  ///< Constant multiplier for the `b` polynomial component.
	 double *con_k,  ///< Constant multiplier for the `c` polynomial component.
	 double *con_b,  ///< Constant multiplier for the additional `b` component.
	 double *con_c   ///< Constant multiplier for the additional `c` component.
	);

static const struct const_Matrix_d* constructor_basis_si_orthonormal_def
	(const int p_b, const struct const_Matrix_d*const rst)
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

// Level 3 ********************************************************************************************************** //

void set_scaling_basis_tri (const int i, const int j, const double b, double *con_i, double *con_j, double *con_b)
{
	*con_i = sqrt((2.0*i+1.0)/2.0);
	*con_j = sqrt((i+j+1.0)/pow(2.0,2.0*i+1.0));
	*con_b = pow(1.0-b,i);
}

void set_scaling_basis_tet
	(const int i, const int j, const int k, const double b, const double c,
         double *con_i, double *con_j, double *con_k, double *con_b, double *con_c)
{
	set_scaling_basis_tri(i,j,b,con_i,con_j,con_b);
	*con_k = sqrt((2.0*(i+j+k)+3.0)/pow(2.0,2.0*(i+j)+3.0));
	*con_c = pow(1.0-c,i+j);
}
