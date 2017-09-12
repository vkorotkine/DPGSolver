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

#include "macros.h"
#include "definitions_alloc.h"
#include "definitions_cubature.h"

#include "matrix.h"

#include "cubature.h"

// Static function declarations ************************************************************************************* //

/// \brief Provides unit tests for the orthonormal tensor-product basis functions.
static void test_unit_basis_tensor_product_orthonormal
	(struct Test_Info*const test_info ///< \ref Test_Info.
	);

// Interface functions ********************************************************************************************** //

void test_unit_bases (struct Test_Info*const test_info)
{
	test_unit_basis_tensor_product_orthonormal(test_info);
//	test_unit_basis_tensor_product_bezier(test_info);
//	test_unit_basis_simplex(test_info);
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

// Tensor-Product Bezier ******************************************************************************************** //
// Simplex Orthonormal ********************************************************************************************** //
// Pyramid Orthonormal ********************************************************************************************** //

static void test_unit_basis_tensor_product_orthonormal (struct Test_Info*const test_info)
{
	char* test_name = "basis_tp_orthonormal";
	char test_name_full[STRLEN_MAX];
	strcpy(test_name_full,"Bases - ");
	strcat(test_name_full,test_name);

	bool pass = true;

	struct Basis_Data_TP_Ortho* b_data_a = constructor_Basis_Data_TP_Ortho('a'), // destructed
	                          * b_data_c = constructor_Basis_Data_TP_Ortho('c'); // destructed

	pass = false;

	destructor_Basis_Data_TP_Ortho(b_data_a);
	destructor_Basis_Data_TP_Ortho(b_data_c);

	test_increment_and_print(test_info,pass,test_name);
}

// Level 1 ********************************************************************************************************** //
// Tensor-Product Orthonormal *************************************************************************************** //

/** \brief Constructor for an polynomial operator for the tensor-product orthonomal basis.
 *  \return Standard. */
static const struct const_Matrix_d* constructor_basis_tp_orthonormal_13
	(const struct const_Matrix_d*const rst ///< The nodes at which the basis functions are evaluated.
	);

static struct Basis_Data_TP_Ortho* constructor_Basis_Data_TP_Ortho (const char eval_type)
{
	struct Basis_Data_TP_Ortho* b_data = calloc(1,sizeof *b_data); // returned

	const struct const_Cubature* d1_p4_GLL = constructor_const_Cubature_tp(1,4,CUB_GLL), // destructed
	                           * d2_p4_GLL = constructor_const_Cubature_tp(2,4,CUB_GLL), // destructed
	                           * d3_p4_GLL = constructor_const_Cubature_tp(3,4,CUB_GLL); // destructed
	if (eval_type == 'a') {
		b_data->phi13 = constructor_basis_tp_orthonormal_13(d1_p4_GLL->rst); // keep
//		b_data->phi22 = constructor_basis_tp_orthonormal_22(d2_p4_GLL->rst); // keep
//		b_data->phi31 = constructor_basis_tp_orthonormal_31(d3_p4_GLL->rst); // keep
	} else if (eval_type == 'c') {
//		b_data->phi13 = constructor_basis_tp_orthonormal(3,d1_p3_GLL->rst); // keep
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
//	destructor_const_Matrix_d(b_data->phi22);
//	destructor_const_Matrix_d(b_data->phi31);
	free(b_data);
}

// Level 2 ********************************************************************************************************** //
// Tensor-Product Orthonormal *************************************************************************************** //

static const struct const_Matrix_d* constructor_basis_tp_orthonormal_13 (const struct const_Matrix_d*const rst)
{
	const int d = 1;
	const ptrdiff_t n_b = pow(3+1,d);

	if (rst->layout != 'C' || rst->ext_1 != d)
		EXIT_UNSUPPORTED;

	const int n_n = rst->ext_0;
	const double*const r = get_col_const_Matrix_d(0,rst);

	struct Matrix_d* phi_rst = constructor_empty_Matrix_d('R',n_n,n_b); // returned

	double* phi_data = phi_rst->data;
	for (ptrdiff_t i = 0; i < n_n; ++i) {
		const double r_i = r[i];

		*phi_data++ = sqrt(1.0/2.0)*1.0;
		*phi_data++ = sqrt(3.0/2.0)*r_i;
		*phi_data++ = sqrt(5.0/2.0)*1.0/2.0*(3.0*pow(r_i,2.0)-1);
		*phi_data++ = sqrt(7.0/2.0)*1.0/2.0*(5.0*pow(r_i,3.0)-3.0*r_i);
	}

	return (const struct const_Matrix_d*) phi_rst;
}

static const struct const_Matrix_d* constructor_basis_tp_orthonormal_22 (const struct const_Matrix_d*const rst)
{
	const int d = 2;
	const ptrdiff_t n_b = pow(2+1,d);

	if (rst->layout != 'C' || rst->ext_1 != d)
		EXIT_UNSUPPORTED;

	const int n_n = rst->ext_0;
	const double*const r = get_col_const_Matrix_d(0,rst);
	            *const s = get_col_const_Matrix_d(1,rst);

	struct Matrix_d* phi_rst = constructor_empty_Matrix_d('R',n_n,n_b); // returned

	double* phi_data = phi_rst->data;
	for (ptrdiff_t i = 0; i < n_n; ++i) {
		const double r_i = r[i],
		             s_i = s[i];

		*phi_data++ = sqrt(1.0/2.0)*1.0*sqrt(1.0/2.0)*1.0;
		*phi_data++ = sqrt(3.0/2.0)*r_i*sqrt(1.0/2.0)*1.0;
		*phi_data++ = sqrt(5.0/2.0)*1.0/2.0*(3.0*pow(r_i,2.0)-1.0)*sqrt(1.0/2.0)*1.0;
		*phi_data++ = sqrt(1.0/2.0)*1.0*sqrt(3.0/2.0)*s_i;
		*phi_data++ = sqrt(3.0/2.0)*r_i*sqrt(3.0/2.0)*s_i;
		*phi_data++ = sqrt(5.0/2.0)*1.0/2.0*(3.0*pow(r_i,2.0)-1.0)*sqrt(3.0/2.0)*s_i;
		*phi_data++ = sqrt(1.0/2.0)*1.0*sqrt(5.0/2.0)*1.0/2.0*(3.0*pow(s_i,2.0)-1.0);
		*phi_data++ = sqrt(3.0/2.0)*r_i*sqrt(5.0/2.0)*1.0/2.0*(3.0*pow(s_i,2.0)-1.0);
		*phi_data++ = sqrt(5.0/2.0)*1.0/2.0*(3.0*pow(r_i,2.0)-1.0)*sqrt(5.0/2.0)*1.0/2.0*(3.0*pow(s_i,2.0)-1.0);
	}

	return (const struct const_Matrix_d*) phi_rst;
}
