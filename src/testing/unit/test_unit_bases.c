// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/** \file
 */

#include "test_unit_bases.h"

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>

#include "test_base.h"
#include "test_support.h"
#include "test_support_bases.h"
#include "test_support_multiarray.h"
#include "test_support_matrix.h"
#include "test_support_vector.h"

#include "macros.h"
#include "definitions_alloc.h"
#include "definitions_cubature.h"
#include "definitions_tol.h"
#include "definitions_elements.h"

#include "multiarray.h"
#include "matrix.h"
#include "vector.h"

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

/// \brief Provides unit tests for the orthonormal pyramid basis functions.
static void test_unit_basis_pyramid_orthonormal
	(struct Test_Info*const test_info ///< \ref Test_Info.
	);

/// \brief Provides unit tests for the bezier tensor-product basis functions.
static void test_unit_basis_tensor_product_bezier
	(struct Test_Info*const test_info ///< \ref Test_Info.
	);

// Interface functions ********************************************************************************************** //

void test_unit_bases (struct Test_Info*const test_info)
{
	test_unit_basis_tensor_product_orthonormal(test_info);
	test_unit_basis_simplex_orthonormal(test_info);
	test_unit_basis_pyramid_orthonormal(test_info);
	test_unit_basis_tensor_product_bezier(test_info);
//	test_unit_basis_simplex_bezier(test_info);
//	test_unit_basis_pyramid_bezier(test_info);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
// Tensor-Product Orthonormal *************************************************************************************** //

/// Container for tensor-product orthonormal basis data to be tested for comparison with expected values.
struct Basis_Data_TP_Ortho {
	const struct const_Matrix_d* phi13, ///< The 1d basis functions of order 3.
	                           * phi22, ///< The 2d basis functions of order 2.
	                           * phi31; ///< The 3d basis functions of order 1.

	const struct const_Multiarray_Matrix_d* grad_phi13, ///< The 1d basis gradient functions of order 3.
	                                      * grad_phi22, ///< The 2d basis gradient functions of order 2.
	                                      * grad_phi31; ///< The 3d basis gradient functions of order 1.

	const struct const_Matrix_d* m_14, ///< The 1d basis mass matrix of order 4.
	                           * m_24, ///< The 2d basis mass matrix of order 4.
	                           * m_34; ///< The 3d basis mass matrix of order 4.

	/// Gradient coefficients for the approximation of a function for the 1d basis of order 4.
	const struct const_Multiarray_d* grad_coef_14;
	/// Gradient coefficients for the approximation of a function for the 2d basis of order 4.
	const struct const_Multiarray_d* grad_coef_24;
	/// Gradient coefficients for the approximation of a function for the 3d basis of order 4.
	const struct const_Multiarray_d* grad_coef_34;
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
	bool    pass        = false;
	double* tol         = NULL;
	bool*   differences = NULL;

	struct Basis_Data_TP_Ortho* b_data_a = constructor_Basis_Data_TP_Ortho('a'), // destructed
	                          * b_data_c = constructor_Basis_Data_TP_Ortho('c'); // destructed

	tol = (double[]) { EPS, 2*EPS, EPS, };
	differences = (bool[])
		{ diff_const_Matrix_d(b_data_a->phi13,b_data_c->phi13,tol[0]),
		  diff_const_Matrix_d(b_data_a->phi22,b_data_c->phi22,tol[1]),
		  diff_const_Matrix_d(b_data_a->phi31,b_data_c->phi31,tol[2]),
		};
	if (check_diff(3,differences,&pass)) {
		if (differences[0]) print_diff_const_Matrix_d(b_data_a->phi13,b_data_c->phi13,tol[0]);
		if (differences[1]) print_diff_const_Matrix_d(b_data_a->phi22,b_data_c->phi22,tol[1]);
		if (differences[2]) print_diff_const_Matrix_d(b_data_a->phi31,b_data_c->phi31,tol[2]);
	}
	test_increment_and_print(test_info,pass,"Bases - tp_orthonormal_basis");

	tol = (double[]) { EPS, 2*EPS, EPS, };
	differences = (bool[])
		{ diff_const_Multiarray_Matrix_d(b_data_a->grad_phi13,b_data_c->grad_phi13,tol[0]),
		  diff_const_Multiarray_Matrix_d(b_data_a->grad_phi22,b_data_c->grad_phi22,tol[1]),
		  diff_const_Multiarray_Matrix_d(b_data_a->grad_phi31,b_data_c->grad_phi31,tol[2]),
		};
	if (check_diff(3,differences,&pass)) {
		if (differences[0])
			print_diff_const_Multiarray_Matrix_d(b_data_a->grad_phi13,b_data_c->grad_phi13,tol[0]);
		if (differences[1])
			print_diff_const_Multiarray_Matrix_d(b_data_a->grad_phi22,b_data_c->grad_phi22,tol[1]);
		if (differences[2])
			print_diff_const_Multiarray_Matrix_d(b_data_a->grad_phi31,b_data_c->grad_phi31,tol[2]);
	}
	test_increment_and_print(test_info,pass,"Bases - tp_orthonormal_grad_basis");

	tol = (double[]) { EPS, EPS, 2*EPS, };
	differences = (bool[])
		{ diff_const_Matrix_d(b_data_a->m_14,b_data_c->m_14,tol[0]),
		  diff_const_Matrix_d(b_data_a->m_24,b_data_c->m_24,tol[1]),
		  diff_const_Matrix_d(b_data_a->m_34,b_data_c->m_34,tol[2]),
		};
	if (check_diff(3,differences,&pass)) {
		if (differences[0]) print_diff_const_Matrix_d(b_data_a->m_14,b_data_c->m_14,tol[0]);
		if (differences[1]) print_diff_const_Matrix_d(b_data_a->m_24,b_data_c->m_24,tol[1]);
		if (differences[2]) print_diff_const_Matrix_d(b_data_a->m_34,b_data_c->m_34,tol[2]);
	}
	test_increment_and_print(test_info,pass,"Bases - tp_orthonormal_mass_matrix");

	tol = (double[]) { EPS, 5*EPS, 30*EPS, };
	differences = (bool[])
		{ diff_const_Multiarray_d(b_data_a->grad_coef_14,b_data_c->grad_coef_14,tol[0]),
		  diff_const_Multiarray_d(b_data_a->grad_coef_24,b_data_c->grad_coef_24,tol[1]),
		  diff_const_Multiarray_d(b_data_a->grad_coef_34,b_data_c->grad_coef_34,tol[2]),
		};
	if (check_diff(3,differences,&pass)) {
		if (differences[0]) print_diff_const_Multiarray_d(b_data_a->grad_coef_14,b_data_c->grad_coef_14,tol[0]);
		if (differences[1]) print_diff_const_Multiarray_d(b_data_a->grad_coef_24,b_data_c->grad_coef_24,tol[1]);
		if (differences[2]) print_diff_const_Multiarray_d(b_data_a->grad_coef_34,b_data_c->grad_coef_34,tol[2]);
	}
	test_increment_and_print(test_info,pass,"Bases - tp_orthonormal_grad_basis_computation");

	destructor_Basis_Data_TP_Ortho(b_data_a);
	destructor_Basis_Data_TP_Ortho(b_data_c);
}

// Simplex Orthonormal ********************************************************************************************** //

/// Container for simplex orthonormal basis data to be tested for comparison with expected values.
struct Basis_Data_SI_Ortho {
	const struct const_Matrix_d* phi22, ///< The 2d basis functions of order 2.
	                           * phi23, ///< The 2d basis functions of order 3.
	                           * phi32; ///< The 3d basis functions of order 2.

	const struct const_Multiarray_Matrix_d* grad_phi22, ///< The 2d basis gradient functions of order 2.
	                                      * grad_phi31; ///< The 3d basis gradient functions of order 1.

	const struct const_Matrix_d* m_24, ///< The 2d basis mass matrix of order 4.
	                           * m_34; ///< The 3d basis mass matrix of order 4.
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
	bool    pass        = false;
	double* tol         = NULL;
	bool*   differences = NULL;

	struct Basis_Data_SI_Ortho* b_data_a = constructor_Basis_Data_SI_Ortho('a'), // destructed
	                          * b_data_c = constructor_Basis_Data_SI_Ortho('c'); // destructed

	tol = (double[]) { 10*EPS, 20*EPS, 10*EPS, };
	differences = (bool[])
		{ diff_const_Matrix_d(b_data_a->phi22,b_data_c->phi22,tol[0]),
		  diff_const_Matrix_d(b_data_a->phi23,b_data_c->phi23,tol[1]),
		  diff_const_Matrix_d(b_data_a->phi32,b_data_c->phi32,tol[2]),
		};
	if (check_diff(3,differences,&pass)) {
		if (differences[0]) print_diff_const_Matrix_d(b_data_a->phi22,b_data_c->phi22,tol[0]);
		if (differences[1]) print_diff_const_Matrix_d(b_data_a->phi23,b_data_c->phi23,tol[1]);
		if (differences[2]) print_diff_const_Matrix_d(b_data_a->phi32,b_data_c->phi32,tol[2]);
	}
	test_increment_and_print(test_info,pass,"Bases - si_orthonormal_basis");

	tol = (double[]) { 2*EPS, EPS,};
	differences = (bool[])
		{ diff_const_Multiarray_Matrix_d(b_data_a->grad_phi22,b_data_c->grad_phi22,tol[0]),
		  diff_const_Multiarray_Matrix_d(b_data_a->grad_phi31,b_data_c->grad_phi31,tol[1]),
		};
	if (check_diff(2,differences,&pass)) {
		if (differences[0])
			print_diff_const_Multiarray_Matrix_d(b_data_a->grad_phi22,b_data_c->grad_phi22,tol[0]);
		if (differences[1])
			print_diff_const_Multiarray_Matrix_d(b_data_a->grad_phi31,b_data_c->grad_phi31,tol[1]);
	}
	test_increment_and_print(test_info,pass,"Bases - si_orthonormal_grad_basis");

	tol = (double[]) { 20*EPS, 30*EPS, };
	differences = (bool[])
		{ diff_const_Matrix_d(b_data_a->m_24, b_data_c->m_24, tol[0]),
		  diff_const_Matrix_d(b_data_a->m_34, b_data_c->m_34, tol[1]),
		};
	if (check_diff(2,differences,&pass)) {
		if (differences[0]) print_diff_const_Matrix_d(b_data_a->m_24,b_data_c->m_24,tol[0]);
		if (differences[1]) print_diff_const_Matrix_d(b_data_a->m_34,b_data_c->m_34,tol[1]);
	}
	test_increment_and_print(test_info,pass,"Bases - si_orthonormal_mass_matrix");

	test_print_warning(test_info,"Not performing grad_basis_si projection test.");

	destructor_Basis_Data_SI_Ortho(b_data_a);
	destructor_Basis_Data_SI_Ortho(b_data_c);
}

// Pyramid Orthonormal ********************************************************************************************** //

/// Container for pyramid orthonormal basis data to be tested for comparison with expected values.
struct Basis_Data_PYR_Ortho {
	const struct const_Matrix_d* phi32; ///< The 3d basis functions of order 2.

	const struct const_Multiarray_Matrix_d* grad_phi32; ///< The 3d basis gradient functions of order 2.

	const struct const_Matrix_d* m_34; ///< The 3d basis mass matrix of order 4.
};

/** \brief Constructor for \ref Basis_Data_PYR_Ortho.
 *  \return Standard. */
static struct Basis_Data_PYR_Ortho* constructor_Basis_Data_PYR_Ortho
	(const char eval_type ///< Method to use to obtain the data. Options: 'a'nalytical, 'c'ompute.
	);

/// \brief Destructor for \ref Basis_Data_PYR_Ortho.
static void destructor_Basis_Data_PYR_Ortho
	(struct Basis_Data_PYR_Ortho* b_data ///< Standard.
	);

static void test_unit_basis_pyramid_orthonormal (struct Test_Info*const test_info)
{
	bool    pass        = false;
	double* tol         = NULL;
	bool*   differences = NULL;

	struct Basis_Data_PYR_Ortho* b_data_a = constructor_Basis_Data_PYR_Ortho('a'), // destructed
	                           * b_data_c = constructor_Basis_Data_PYR_Ortho('c'); // destructed

	tol = (double[]) { 3*EPS, };
	differences = (bool[])
		{ diff_const_Matrix_d(b_data_a->phi32,b_data_c->phi32,tol[0]),
		};
	if (check_diff(1,differences,&pass)) {
		if (differences[0]) print_diff_const_Matrix_d(b_data_a->phi32,b_data_c->phi32,tol[0]);
	}
	test_increment_and_print(test_info,pass,"Bases - pyr_orthonormal_basis");

	tol = (double[]) { 2*EPS, };
	differences = (bool[])
		{ diff_const_Multiarray_Matrix_d(b_data_a->grad_phi32,b_data_c->grad_phi32,tol[1]),
		};
	if (check_diff(1,differences,&pass)) {
		if (differences[0])
			print_diff_const_Multiarray_Matrix_d(b_data_a->grad_phi32,b_data_c->grad_phi32,tol[0]);
	}
	test_increment_and_print(test_info,pass,"Bases - pyr_orthonormal_grad_basis");

	tol = (double[]) { 2e2*EPS };
	differences = (bool[])
		{ diff_const_Matrix_d(b_data_a->m_34, b_data_c->m_34, tol[0]),
		};
	if (check_diff(1,differences,&pass)) {
		if (differences[0]) print_diff_const_Matrix_d(b_data_a->m_34,b_data_c->m_34,tol[0]);
	}
	test_increment_and_print(test_info,pass,"Bases - pyr_orthonormal_mass_matrix");

	test_print_warning(test_info,"Not performing grad_basis_pyr projection test.");

	destructor_Basis_Data_PYR_Ortho(b_data_a);
	destructor_Basis_Data_PYR_Ortho(b_data_c);
}

// Tensor-Product Bezier ******************************************************************************************** //

/// Container for tensor-product bezier basis data to be tested for comparison with expected values.
struct Basis_Data_TP_Bezier {
	const struct const_Matrix_d* phi13, ///< The 1d basis functions of order 3.
	                           * phi22, ///< The 2d basis functions of order 2.
	                           * phi32; ///< The 3d basis functions of order 1.

	const struct const_Multiarray_Matrix_d* grad_phi13, ///< The 1d basis gradient functions of order 3.
	                                      * grad_phi22, ///< The 2d basis gradient functions of order 2.
	                                      * grad_phi31; ///< The 3d basis gradient functions of order 1.

	const struct const_Vector_d* p_14, ///< The 1d basis partition of unity (summed) vector of order 4.
	                           * p_24, ///< The 2d basis partition of unity (summed) vector of order 4.
	                           * p_34; ///< The 3d basis partition of unity (summed) vector of order 4.
};

/** \brief Constructor for \ref Basis_Data_TP_Bezier.
 *  \return Standard. */
static struct Basis_Data_TP_Bezier* constructor_Basis_Data_TP_Bezier
	(const char eval_type ///< Method to use to obtain the data. Options: 'a'nalytical, 'c'ompute.
	);

/// \brief Destructor for \ref Basis_Data_TP_Bezier.
static void destructor_Basis_Data_TP_Bezier
	(struct Basis_Data_TP_Bezier* b_data ///< Standard.
	);

static void test_unit_basis_tensor_product_bezier (struct Test_Info*const test_info)
{
	bool    pass        = false;
	double* tol         = NULL;
	bool*   differences = NULL;

	struct Basis_Data_TP_Bezier* b_data_a = constructor_Basis_Data_TP_Bezier('a'), // destructed
	                           * b_data_c = constructor_Basis_Data_TP_Bezier('c'); // destructed

	tol = (double[]) { 2*EPS, 2*EPS, 3*EPS, };
	differences = (bool[])
		{ diff_const_Matrix_d(b_data_a->phi13,b_data_c->phi13,tol[0]),
		  diff_const_Matrix_d(b_data_a->phi22,b_data_c->phi22,tol[1]),
		  diff_const_Matrix_d(b_data_a->phi32,b_data_c->phi32,tol[2]),
		};
	if (check_diff(3,differences,&pass)) {
		if (differences[0]) print_diff_const_Matrix_d(b_data_a->phi13,b_data_c->phi13,tol[0]);
		if (differences[1]) print_diff_const_Matrix_d(b_data_a->phi22,b_data_c->phi22,tol[1]);
		if (differences[2]) print_diff_const_Matrix_d(b_data_a->phi32,b_data_c->phi32,tol[2]);
	}
	test_increment_and_print(test_info,pass,"Bases - tp_bezier_basis");

	tol = (double[]) { EPS, 2*EPS, EPS, };
	differences = (bool[])
		{ diff_const_Multiarray_Matrix_d(b_data_a->grad_phi13,b_data_c->grad_phi13,tol[0]),
		  diff_const_Multiarray_Matrix_d(b_data_a->grad_phi22,b_data_c->grad_phi22,tol[1]),
		  diff_const_Multiarray_Matrix_d(b_data_a->grad_phi31,b_data_c->grad_phi31,tol[2]),
		};
	if (check_diff(3,differences,&pass)) {
		if (differences[0])
			print_diff_const_Multiarray_Matrix_d(b_data_a->grad_phi13,b_data_c->grad_phi13,tol[0]);
		if (differences[1])
			print_diff_const_Multiarray_Matrix_d(b_data_a->grad_phi22,b_data_c->grad_phi22,tol[1]);
		if (differences[2])
			print_diff_const_Multiarray_Matrix_d(b_data_a->grad_phi31,b_data_c->grad_phi31,tol[2]);
	}
	test_increment_and_print(test_info,pass,"Bases - tp_bezier_grad_basis");

	tol = (double[]) { EPS, EPS, EPS, };
	differences = (bool[])
		{ diff_const_Vector_d(b_data_a->p_14,b_data_c->p_14,tol[0]),
		  diff_const_Vector_d(b_data_a->p_24,b_data_c->p_24,tol[1]),
		  diff_const_Vector_d(b_data_a->p_34,b_data_c->p_34,tol[2]),
		};
	if (check_diff(3,differences,&pass)) {
		if (differences[0]) print_diff_const_Vector_d(b_data_a->p_14,b_data_c->p_14,tol[0]);
		if (differences[1]) print_diff_const_Vector_d(b_data_a->p_24,b_data_c->p_24,tol[1]);
		if (differences[2]) print_diff_const_Vector_d(b_data_a->p_34,b_data_c->p_34,tol[2]);
	}
	test_increment_and_print(test_info,pass,"Bases - tp_bezier_part_unity");

	destructor_Basis_Data_TP_Bezier(b_data_a);
	destructor_Basis_Data_TP_Bezier(b_data_c);

	test_print_warning(test_info,"Not performing grad_basis_tp_bezier projection test.");
test_print_warning(test_info,"Ensure that the basis on the std element is being used.");
}

// Level 1 ********************************************************************************************************** //
// Tensor-Product Orthonormal *************************************************************************************** //

static struct Basis_Data_TP_Ortho* constructor_Basis_Data_TP_Ortho (const char eval_type)
{
	struct Basis_Data_TP_Ortho* b_data = calloc(1,sizeof *b_data); // returned

	const struct const_Cubature* d1_p4_GLL = constructor_const_Cubature_tp(1,4,CUB_GLL), // destructed
	                           * d2_p4_GLL = constructor_const_Cubature_tp(2,4,CUB_GLL), // destructed
	                           * d3_p4_GLL = constructor_const_Cubature_tp(3,4,CUB_GLL); // destructed
	const int   super_type = ST_TP;
	const char* basis_name = "tp_ortho";
	if (eval_type == 'a') {
		b_data->phi13        = constructor_basis_tp_orthonormal_def(3,d1_p4_GLL->rst);      // keep
		b_data->phi22        = constructor_basis_tp_orthonormal_def(2,d2_p4_GLL->rst);      // keep
		b_data->phi31        = constructor_basis_tp_orthonormal_def(1,d3_p4_GLL->rst);      // keep
		b_data->grad_phi13   = constructor_grad_basis_tp_orthonormal_def(3,d1_p4_GLL->rst); // keep
		b_data->grad_phi22   = constructor_grad_basis_tp_orthonormal_def(2,d2_p4_GLL->rst); // keep
		b_data->grad_phi31   = constructor_grad_basis_tp_orthonormal_def(1,d3_p4_GLL->rst); // keep
		b_data->m_14         = constructor_mass_orthonormal_def(1,4,super_type);            // keep
		b_data->m_24         = constructor_mass_orthonormal_def(2,4,super_type);            // keep
		b_data->m_34         = constructor_mass_orthonormal_def(3,4,super_type);            // keep
		b_data->grad_coef_14 = constructor_grad_vals_computation_def(1,4,basis_name);       // keep
		b_data->grad_coef_24 = constructor_grad_vals_computation_def(2,4,basis_name);       // keep
		b_data->grad_coef_34 = constructor_grad_vals_computation_def(3,4,basis_name);       // keep
	} else if (eval_type == 'c') {
		b_data->phi13        = constructor_basis_tp_orthonormal(3,d1_p4_GLL->rst);      // keep
		b_data->phi22        = constructor_basis_tp_orthonormal(2,d2_p4_GLL->rst);      // keep
		b_data->phi31        = constructor_basis_tp_orthonormal(1,d3_p4_GLL->rst);      // keep
		b_data->grad_phi13   = constructor_grad_basis_tp_orthonormal(3,d1_p4_GLL->rst); // keep
		b_data->grad_phi22   = constructor_grad_basis_tp_orthonormal(2,d2_p4_GLL->rst); // keep
		b_data->grad_phi31   = constructor_grad_basis_tp_orthonormal(1,d3_p4_GLL->rst); // keep
		b_data->m_14         = constructor_mass_orthonormal(1,4,super_type);            // keep
		b_data->m_24         = constructor_mass_orthonormal(2,4,super_type);            // keep
		b_data->m_34         = constructor_mass_orthonormal(3,4,super_type);            // keep
		b_data->grad_coef_14 = constructor_grad_vals_computation(1,4,basis_name);       // keep
		b_data->grad_coef_24 = constructor_grad_vals_computation(2,4,basis_name);       // keep
		b_data->grad_coef_34 = constructor_grad_vals_computation(3,4,basis_name);       // keep
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
	destructor_const_Multiarray_Matrix_d(b_data->grad_phi13);
	destructor_const_Multiarray_Matrix_d(b_data->grad_phi22);
	destructor_const_Multiarray_Matrix_d(b_data->grad_phi31);
	destructor_const_Matrix_d(b_data->m_14);
	destructor_const_Matrix_d(b_data->m_24);
	destructor_const_Matrix_d(b_data->m_34);
	destructor_const_Multiarray_d(b_data->grad_coef_14);
	destructor_const_Multiarray_d(b_data->grad_coef_24);
	destructor_const_Multiarray_d(b_data->grad_coef_34);
	free(b_data);
}

// Simplex Orthonormal ********************************************************************************************** //

static struct Basis_Data_SI_Ortho* constructor_Basis_Data_SI_Ortho (const char eval_type)
{
	struct Basis_Data_SI_Ortho* b_data = calloc(1,sizeof *b_data); // returned

	const struct const_Cubature* d2_p4_AO = constructor_const_Cubature_si(2,4,CUB_AO), // destructed
	                           * d3_p4_AO = constructor_const_Cubature_si(3,4,CUB_AO); // destructed
	const int super_type = ST_SI;
	const int p_mass = 4;
	if (eval_type == 'a') {
		b_data->phi22      = constructor_basis_si_orthonormal_def(2,d2_p4_AO->rst);      // keep
		b_data->phi23      = constructor_basis_si_orthonormal_def(3,d2_p4_AO->rst);      // keep
		b_data->phi32      = constructor_basis_si_orthonormal_def(2,d3_p4_AO->rst);      // keep
		b_data->grad_phi22 = constructor_grad_basis_si_orthonormal_def(2,d2_p4_AO->rst); // keep
		b_data->grad_phi31 = constructor_grad_basis_si_orthonormal_def(1,d3_p4_AO->rst); // keep
		b_data->m_24       = constructor_mass_orthonormal_def(2,p_mass,super_type);      // keep
		b_data->m_34       = constructor_mass_orthonormal_def(3,p_mass,super_type);      // keep
	} else if (eval_type == 'c') {
		b_data->phi22      = constructor_basis_si_orthonormal(2,d2_p4_AO->rst);      // keep
		b_data->phi23      = constructor_basis_si_orthonormal(3,d2_p4_AO->rst);      // keep
		b_data->phi32      = constructor_basis_si_orthonormal(2,d3_p4_AO->rst);      // keep
		b_data->grad_phi22 = constructor_grad_basis_si_orthonormal(2,d2_p4_AO->rst); // keep
		b_data->grad_phi31 = constructor_grad_basis_si_orthonormal(1,d3_p4_AO->rst); // keep
		b_data->m_24       = constructor_mass_orthonormal(2,p_mass,super_type);      // keep
		b_data->m_34       = constructor_mass_orthonormal(3,p_mass,super_type);      // keep
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
	destructor_const_Multiarray_Matrix_d(b_data->grad_phi22);
	destructor_const_Multiarray_Matrix_d(b_data->grad_phi31);
	destructor_const_Matrix_d(b_data->m_24);
	destructor_const_Matrix_d(b_data->m_34);
	free(b_data);
}

// Pyramid Orthonormal ********************************************************************************************** //

static struct Basis_Data_PYR_Ortho* constructor_Basis_Data_PYR_Ortho (const char eval_type)
{
	struct Basis_Data_PYR_Ortho* b_data = calloc(1,sizeof *b_data); // returned

	const struct const_Cubature* d3_p4_GLL = constructor_const_Cubature_pyr(3,4,CUB_GLL); // destructed
	const int super_type = ST_PYR;
	const int p_mass = 4;
	if (eval_type == 'a') {
		b_data->phi32      = constructor_basis_pyr_orthonormal_def(2,d3_p4_GLL->rst);      // keep
		b_data->grad_phi32 = constructor_grad_basis_pyr_orthonormal_def(2,d3_p4_GLL->rst); // keep
		b_data->m_34       = constructor_mass_orthonormal_def(3,p_mass,super_type);        // keep
	} else if (eval_type == 'c') {
		b_data->phi32      = constructor_basis_pyr_orthonormal(2,d3_p4_GLL->rst);      // keep
		b_data->grad_phi32 = constructor_grad_basis_pyr_orthonormal(2,d3_p4_GLL->rst); // keep
		b_data->m_34       = constructor_mass_orthonormal(3,p_mass,super_type);        // keep
	} else {
		EXIT_UNSUPPORTED;
	}
	destructor_const_Cubature(d3_p4_GLL);

	return b_data;
}

static void destructor_Basis_Data_PYR_Ortho (struct Basis_Data_PYR_Ortho* b_data)
{
	destructor_const_Matrix_d(b_data->phi32);
	destructor_const_Multiarray_Matrix_d(b_data->grad_phi32);
	destructor_const_Matrix_d(b_data->m_34);
	free(b_data);
}

// Tensor-Product Orthonormal *************************************************************************************** //

static struct Basis_Data_TP_Bezier* constructor_Basis_Data_TP_Bezier (const char eval_type)
{
	struct Basis_Data_TP_Bezier* b_data = calloc(1,sizeof *b_data); // returned

	const struct const_Cubature* d1_p4_GLL = constructor_const_Cubature_tp(1,4,CUB_GLL), // destructed
	                           * d2_p4_GLL = constructor_const_Cubature_tp(2,4,CUB_GLL), // destructed
	                           * d3_p4_GLL = constructor_const_Cubature_tp(3,4,CUB_GLL); // destructed
	if (eval_type == 'a') {
		b_data->phi13      = constructor_basis_tp_bezier_def(3,d1_p4_GLL->rst);      // keep
		b_data->phi22      = constructor_basis_tp_bezier_def(2,d2_p4_GLL->rst);      // keep
		b_data->phi32      = constructor_basis_tp_bezier_def(2,d3_p4_GLL->rst);      // keep
		b_data->grad_phi13 = constructor_grad_basis_tp_bezier_def(3,d1_p4_GLL->rst); // keep
		b_data->grad_phi22 = constructor_grad_basis_tp_bezier_def(2,d2_p4_GLL->rst); // keep
		b_data->grad_phi31 = constructor_grad_basis_tp_bezier_def(1,d3_p4_GLL->rst); // keep
		b_data->p_14       = constructor_part_unity_def(d1_p4_GLL->rst->ext_0);      // keep
		b_data->p_24       = constructor_part_unity_def(d2_p4_GLL->rst->ext_0);      // keep
		b_data->p_34       = constructor_part_unity_def(d3_p4_GLL->rst->ext_0);      // keep
	} else if (eval_type == 'c') {
		b_data->phi13      = constructor_basis_tp_bezier(3,d1_p4_GLL->rst);      // keep
		b_data->phi22      = constructor_basis_tp_bezier(2,d2_p4_GLL->rst);      // keep
		b_data->phi32      = constructor_basis_tp_bezier(2,d3_p4_GLL->rst);      // keep
		b_data->grad_phi13 = constructor_grad_basis_tp_bezier(3,d1_p4_GLL->rst); // keep
		b_data->grad_phi22 = constructor_grad_basis_tp_bezier(2,d2_p4_GLL->rst); // keep
		b_data->grad_phi31 = constructor_grad_basis_tp_bezier(1,d3_p4_GLL->rst); // keep
		b_data->p_14       = constructor_part_unity(b_data->phi13);              // keep
		b_data->p_24       = constructor_part_unity(b_data->phi22);              // keep
		b_data->p_34       = constructor_part_unity(b_data->phi32);              // keep
	} else {
		EXIT_UNSUPPORTED;
	}
	destructor_const_Cubature(d1_p4_GLL);
	destructor_const_Cubature(d2_p4_GLL);
	destructor_const_Cubature(d3_p4_GLL);

	return b_data;
}

static void destructor_Basis_Data_TP_Bezier (struct Basis_Data_TP_Bezier* b_data)
{
	destructor_const_Matrix_d(b_data->phi13);
	destructor_const_Matrix_d(b_data->phi22);
	destructor_const_Matrix_d(b_data->phi32);
	destructor_const_Multiarray_Matrix_d(b_data->grad_phi13);
	destructor_const_Multiarray_Matrix_d(b_data->grad_phi22);
	destructor_const_Multiarray_Matrix_d(b_data->grad_phi31);
	destructor_const_Vector_d(b_data->p_14);
	destructor_const_Vector_d(b_data->p_24);
	destructor_const_Vector_d(b_data->p_34);
	free(b_data);
}
