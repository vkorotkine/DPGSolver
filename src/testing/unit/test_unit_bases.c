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
#include "definitions_elements.h"
#include "definitions_nodes.h"
#include "definitions_tol.h"

#include "multiarray.h"
#include "matrix.h"
#include "vector.h"

#include "bases.h"
#include "nodes.h"

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

/// \brief Provides unit tests for the bezier simplex basis functions.
static void test_unit_basis_simplex_bezier
	(struct Test_Info*const test_info ///< \ref Test_Info.
	);

// Interface functions ********************************************************************************************** //

/** \test Performs unit testing for the bases (\ref test_unit_bases.c).
 *  \return 0 on success. */
int main
	(int nargc,  ///< Standard.
	 char** argv ///< Standard.
	)
{
	assert_condition_message(nargc == 2,"Invalid number of input arguments");
	const char* test_name = argv[1];

	struct Test_Info test_info = { .n_warn = 0, };
	if (strcmp(test_name,"tp_orthonormal") == 0)
		test_unit_basis_tensor_product_orthonormal(&test_info);
	else if (strcmp(test_name,"si_orthonormal") == 0)
		test_unit_basis_simplex_orthonormal(&test_info);
	else if (strcmp(test_name,"pyr_orthonormal") == 0)
		test_unit_basis_pyramid_orthonormal(&test_info);
	else if (strcmp(test_name,"tp_bezier") == 0)
		test_unit_basis_tensor_product_bezier(&test_info);
	else if (strcmp(test_name,"si_bezier") == 0)
		test_unit_basis_simplex_bezier(&test_info);
	else
		EXIT_ERROR("Invalid test name: %s\n",test_name);
//	test_unit_basis_pyramid_bezier(&test_info);

	output_warning_count(&test_info);
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

	/// Gradient coefficients for the approximation of order 3 of a p2 tensor-product polynomial for the 1d basis.
	const struct const_Multiarray_d* grad_coef_13;
	/// Gradient coefficients for the approximation of order 3 of a p2 tensor-product polynomial for the 2d basis.
	const struct const_Multiarray_d* grad_coef_23;
	/// Gradient coefficients for the approximation of order 3 of a p2 tensor-product polynomial for the 3d basis.
	const struct const_Multiarray_d* grad_coef_33;
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
	UNUSED(test_info);
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
	expect_condition(pass,"basis");

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
	expect_condition(pass,"grad basis");

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
	expect_condition(pass,"mass matrix");

	tol = (double[]) { EPS, 2*EPS, 20*EPS, };
	differences = (bool[])
		{ diff_const_Multiarray_d(b_data_a->grad_coef_13,b_data_c->grad_coef_13,tol[0]),
		  diff_const_Multiarray_d(b_data_a->grad_coef_23,b_data_c->grad_coef_23,tol[1]),
		  diff_const_Multiarray_d(b_data_a->grad_coef_33,b_data_c->grad_coef_33,tol[2]),
		};
	if (check_diff(3,differences,&pass)) {
		if (differences[0]) print_diff_const_Multiarray_d(b_data_a->grad_coef_13,b_data_c->grad_coef_13,tol[0]);
		if (differences[1]) print_diff_const_Multiarray_d(b_data_a->grad_coef_23,b_data_c->grad_coef_23,tol[1]);
		if (differences[2]) print_diff_const_Multiarray_d(b_data_a->grad_coef_33,b_data_c->grad_coef_33,tol[2]);
	}
	expect_condition(pass,"gradient evaluation");

	destructor_Basis_Data_TP_Ortho(b_data_a);
	destructor_Basis_Data_TP_Ortho(b_data_c);

	assert_condition(pass);
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

	/// Gradient coefficients for the approximation of order 5 of a p2 tensor-product polynomial for the 2d basis.
	const struct const_Multiarray_d* grad_coef_25;
	/// Gradient coefficients for the approximation of order 7 of a p2 tensor-product polynomial for the 3d basis.
	const struct const_Multiarray_d* grad_coef_37;
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
	UNUSED(test_info);
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
	expect_condition(pass,"basis");

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
	expect_condition(pass,"grad basis");

	tol = (double[]) { 20*EPS, 30*EPS, };
	differences = (bool[])
		{ diff_const_Matrix_d(b_data_a->m_24, b_data_c->m_24, tol[0]),
		  diff_const_Matrix_d(b_data_a->m_34, b_data_c->m_34, tol[1]),
		};
	if (check_diff(2,differences,&pass)) {
		if (differences[0]) print_diff_const_Matrix_d(b_data_a->m_24,b_data_c->m_24,tol[0]);
		if (differences[1]) print_diff_const_Matrix_d(b_data_a->m_34,b_data_c->m_34,tol[1]);
	}
	expect_condition(pass,"mass matrix");

	tol = (double[]) { 9*EPS, 8e4*EPS, };
	differences = (bool[])
		{ diff_const_Multiarray_d(b_data_a->grad_coef_25,b_data_c->grad_coef_25,tol[0]),
		  diff_const_Multiarray_d(b_data_a->grad_coef_37,b_data_c->grad_coef_37,tol[1]),
		};
	if (check_diff(2,differences,&pass)) {
		if (differences[0]) print_diff_const_Multiarray_d(b_data_a->grad_coef_25,b_data_c->grad_coef_25,tol[0]);
		if (differences[1]) print_diff_const_Multiarray_d(b_data_a->grad_coef_37,b_data_c->grad_coef_37,tol[1]);
	}
	expect_condition(pass,"gradient evaluation");

	destructor_Basis_Data_SI_Ortho(b_data_a);
	destructor_Basis_Data_SI_Ortho(b_data_c);

	assert_condition(pass);
}

// Pyramid Orthonormal ********************************************************************************************** //

/// Container for pyramid orthonormal basis data to be tested for comparison with expected values.
struct Basis_Data_PYR_Ortho {
	const struct const_Matrix_d* phi32; ///< The 3d basis functions of order 2.

	const struct const_Multiarray_Matrix_d* grad_phi32; ///< The 3d basis gradient functions of order 2.

	const struct const_Matrix_d* m_34; ///< The 3d basis mass matrix of order 4.

	/// Gradient coefficients for the approximation of order 6 of a p2 tensor-product polynomial for the 3d basis.
	const struct const_Multiarray_d* grad_coef_36;
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
	UNUSED(test_info);
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
	expect_condition(pass,"basis");

	tol = (double[]) { 2*EPS, };
	differences = (bool[])
		{ diff_const_Multiarray_Matrix_d(b_data_a->grad_phi32,b_data_c->grad_phi32,tol[0]),
		};
	if (check_diff(1,differences,&pass)) {
		if (differences[0])
			print_diff_const_Multiarray_Matrix_d(b_data_a->grad_phi32,b_data_c->grad_phi32,tol[0]);
	}
	expect_condition(pass,"grad basis");

	tol = (double[]) { 2e2*EPS };
	differences = (bool[])
		{ diff_const_Matrix_d(b_data_a->m_34, b_data_c->m_34, tol[0]),
		};
	if (check_diff(1,differences,&pass)) {
		if (differences[0]) print_diff_const_Matrix_d(b_data_a->m_34,b_data_c->m_34,tol[0]);
	}
	expect_condition(pass,"mass matrix");

	tol = (double[]) { 5e6*EPS, };
	differences = (bool[])
		{ diff_const_Multiarray_d(b_data_a->grad_coef_36,b_data_c->grad_coef_36,tol[0]),
		};
	if (check_diff(1,differences,&pass)) {
		if (differences[0]) print_diff_const_Multiarray_d(b_data_a->grad_coef_36,b_data_c->grad_coef_36,tol[0]);
	}
	expect_condition(pass,"gradient evaluation");

	destructor_Basis_Data_PYR_Ortho(b_data_a);
	destructor_Basis_Data_PYR_Ortho(b_data_c);

	assert_condition(pass);
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

	/// Gradient coefficients for the approximation of order 3 of a p2 tensor-product polynomial for the 1d basis.
	const struct const_Multiarray_d* grad_coef_13;
	/// Gradient coefficients for the approximation of order 3 of a p2 tensor-product polynomial for the 2d basis.
	const struct const_Multiarray_d* grad_coef_23;
	/// Gradient coefficients for the approximation of order 3 of a p2 tensor-product polynomial for the 3d basis.
	const struct const_Multiarray_d* grad_coef_33;
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
	UNUSED(test_info);
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
	expect_condition(pass,"basis");

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
	expect_condition(pass,"grad basis");

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
	expect_condition(pass,"partition of unity");

	tol = (double[]) { EPS, 20*EPS, 5*EPS, };
	differences = (bool[])
		{ diff_const_Multiarray_d(b_data_a->grad_coef_13,b_data_c->grad_coef_13,tol[0]),
		  diff_const_Multiarray_d(b_data_a->grad_coef_23,b_data_c->grad_coef_23,tol[1]),
		  diff_const_Multiarray_d(b_data_a->grad_coef_33,b_data_c->grad_coef_33,tol[2]),
		};
	if (check_diff(3,differences,&pass)) {
		if (differences[0]) print_diff_const_Multiarray_d(b_data_a->grad_coef_13,b_data_c->grad_coef_13,tol[0]);
		if (differences[1]) print_diff_const_Multiarray_d(b_data_a->grad_coef_23,b_data_c->grad_coef_23,tol[1]);
		if (differences[2]) print_diff_const_Multiarray_d(b_data_a->grad_coef_33,b_data_c->grad_coef_33,tol[2]);
	}
	expect_condition(pass,"gradient evaluation");

	destructor_Basis_Data_TP_Bezier(b_data_a);
	destructor_Basis_Data_TP_Bezier(b_data_c);

	assert_condition(pass);
}

// Simplex Bezier *************************************************************************************************** //

/// Container for simplex bezier basis data to be tested for comparison with expected values.
struct Basis_Data_SI_Bezier {
	const struct const_Matrix_d* phi22, ///< The 2d basis functions of order 2.
	                           * phi23, ///< The 2d basis functions of order 3.
	                           * phi32; ///< The 3d basis functions of order 2.

	const struct const_Multiarray_Matrix_d* grad_phi22, ///< The 2d basis gradient functions of order 2.
	                                      * grad_phi31; ///< The 3d basis gradient functions of order 1.

	const struct const_Vector_d* p_24, ///< The 2d basis partition of unity (summed) vector of order 4.
	                           * p_34; ///< The 3d basis partition of unity (summed) vector of order 4.

	/// Gradient coefficients for the approximation of order 5 of a p2 tensor-product polynomial for the 2d basis.
	const struct const_Multiarray_d* grad_coef_25;
	/// Gradient coefficients for the approximation of order 7 of a p2 tensor-product polynomial for the 3d basis.
	const struct const_Multiarray_d* grad_coef_37;
};

/** \brief Constructor for \ref Basis_Data_SI_Bezier.
 *  \return Standard. */
static struct Basis_Data_SI_Bezier* constructor_Basis_Data_SI_Bezier
	(const char eval_type ///< Method to use to obtain the data. Options: 'a'nalytical, 'c'ompute.
	);

/// \brief Destructor for \ref Basis_Data_SI_Bezier.
static void destructor_Basis_Data_SI_Bezier
	(struct Basis_Data_SI_Bezier* b_data ///< Standard.
	);

static void test_unit_basis_simplex_bezier (struct Test_Info*const test_info)
{
	UNUSED(test_info);
	bool    pass        = false;
	double* tol         = NULL;
	bool*   differences = NULL;

	struct Basis_Data_SI_Bezier* b_data_a = constructor_Basis_Data_SI_Bezier('a'), // destructed
	                           * b_data_c = constructor_Basis_Data_SI_Bezier('c'); // destructed

	tol = (double[]) { EPS, 2*EPS, EPS, };
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
	expect_condition(pass,"basis");

	tol = (double[]) { EPS, EPS, };
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
	expect_condition(pass,"grad basis");

	tol = (double[]) { EPS, EPS, };
	differences = (bool[])
		{ diff_const_Vector_d(b_data_a->p_24,b_data_c->p_24,tol[0]),
		  diff_const_Vector_d(b_data_a->p_34,b_data_c->p_34,tol[1]),
		};
	if (check_diff(2,differences,&pass)) {
		if (differences[0]) print_diff_const_Vector_d(b_data_a->p_24,b_data_c->p_24,tol[0]);
		if (differences[1]) print_diff_const_Vector_d(b_data_a->p_34,b_data_c->p_34,tol[1]);
	}
	expect_condition(pass,"partition of unity");

	tol = (double[]) { 2e1*EPS, 11e4*EPS, };
	differences = (bool[])
		{ diff_const_Multiarray_d(b_data_a->grad_coef_25,b_data_c->grad_coef_25,tol[0]),
		  diff_const_Multiarray_d(b_data_a->grad_coef_37,b_data_c->grad_coef_37,tol[1]),
		};
	if (check_diff(2,differences,&pass)) {
		if (differences[0]) print_diff_const_Multiarray_d(b_data_a->grad_coef_25,b_data_c->grad_coef_25,tol[0]);
		if (differences[1]) print_diff_const_Multiarray_d(b_data_a->grad_coef_37,b_data_c->grad_coef_37,tol[1]);
	}
	expect_condition(pass,"gradient evaluation");

	destructor_Basis_Data_SI_Bezier(b_data_a);
	destructor_Basis_Data_SI_Bezier(b_data_c);

	assert_condition(pass);
}

// Level 1 ********************************************************************************************************** //
// Tensor-Product Orthonormal *************************************************************************************** //

static struct Basis_Data_TP_Ortho* constructor_Basis_Data_TP_Ortho (const char eval_type)
{
	struct Basis_Data_TP_Ortho* b_data = calloc(1,sizeof *b_data); // returned

	const struct const_Nodes* d1_p4_GLL = constructor_const_Nodes_tp(1,4,NODES_GLL), // destructed
	                        * d2_p4_GLL = constructor_const_Nodes_tp(2,4,NODES_GLL), // destructed
	                        * d3_p4_GLL = constructor_const_Nodes_tp(3,4,NODES_GLL); // destructed
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
		b_data->grad_coef_13 = constructor_grad_vals_computation_def(1,3,basis_name);       // keep
		b_data->grad_coef_23 = constructor_grad_vals_computation_def(2,3,basis_name);       // keep
		b_data->grad_coef_33 = constructor_grad_vals_computation_def(3,3,basis_name);       // keep
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
		b_data->grad_coef_13 = constructor_grad_vals_computation(1,3,basis_name);       // keep
		b_data->grad_coef_23 = constructor_grad_vals_computation(2,3,basis_name);       // keep
		b_data->grad_coef_33 = constructor_grad_vals_computation(3,3,basis_name);       // keep
	} else {
		EXIT_UNSUPPORTED;
	}
	destructor_const_Nodes(d1_p4_GLL);
	destructor_const_Nodes(d2_p4_GLL);
	destructor_const_Nodes(d3_p4_GLL);

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
	destructor_const_Multiarray_d(b_data->grad_coef_13);
	destructor_const_Multiarray_d(b_data->grad_coef_23);
	destructor_const_Multiarray_d(b_data->grad_coef_33);
	free(b_data);
}

// Simplex Orthonormal ********************************************************************************************** //

static struct Basis_Data_SI_Ortho* constructor_Basis_Data_SI_Ortho (const char eval_type)
{
	struct Basis_Data_SI_Ortho* b_data = calloc(1,sizeof *b_data); // returned

	const struct const_Nodes* d2_p4_AO = constructor_const_Nodes_si(2,4,NODES_AO), // destructed
	                        * d3_p4_AO = constructor_const_Nodes_si(3,4,NODES_AO); // destructed
	const int super_type = ST_SI;
	const int p_mass = 4;
	const char* basis_name = "si_ortho";
	if (eval_type == 'a') {
		b_data->phi22        = constructor_basis_si_orthonormal_def(2,d2_p4_AO->rst);      // keep
		b_data->phi23        = constructor_basis_si_orthonormal_def(3,d2_p4_AO->rst);      // keep
		b_data->phi32        = constructor_basis_si_orthonormal_def(2,d3_p4_AO->rst);      // keep
		b_data->grad_phi22   = constructor_grad_basis_si_orthonormal_def(2,d2_p4_AO->rst); // keep
		b_data->grad_phi31   = constructor_grad_basis_si_orthonormal_def(1,d3_p4_AO->rst); // keep
		b_data->m_24         = constructor_mass_orthonormal_def(2,p_mass,super_type);      // keep
		b_data->m_34         = constructor_mass_orthonormal_def(3,p_mass,super_type);      // keep
		b_data->grad_coef_25 = constructor_grad_vals_computation_def(2,5,basis_name);      // keep
		b_data->grad_coef_37 = constructor_grad_vals_computation_def(3,7,basis_name);      // keep
	} else if (eval_type == 'c') {
		b_data->phi22        = constructor_basis_si_orthonormal(2,d2_p4_AO->rst);      // keep
		b_data->phi23        = constructor_basis_si_orthonormal(3,d2_p4_AO->rst);      // keep
		b_data->phi32        = constructor_basis_si_orthonormal(2,d3_p4_AO->rst);      // keep
		b_data->grad_phi22   = constructor_grad_basis_si_orthonormal(2,d2_p4_AO->rst); // keep
		b_data->grad_phi31   = constructor_grad_basis_si_orthonormal(1,d3_p4_AO->rst); // keep
		b_data->m_24         = constructor_mass_orthonormal(2,p_mass,super_type);      // keep
		b_data->m_34         = constructor_mass_orthonormal(3,p_mass,super_type);      // keep
		b_data->grad_coef_25 = constructor_grad_vals_computation(2,5,basis_name);      // keep
		b_data->grad_coef_37 = constructor_grad_vals_computation(3,7,basis_name);      // keep
	} else {
		EXIT_UNSUPPORTED;
	}
	destructor_const_Nodes(d2_p4_AO);
	destructor_const_Nodes(d3_p4_AO);

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
	destructor_const_Multiarray_d(b_data->grad_coef_25);
	destructor_const_Multiarray_d(b_data->grad_coef_37);
	free(b_data);
}

// Pyramid Orthonormal ********************************************************************************************** //

static struct Basis_Data_PYR_Ortho* constructor_Basis_Data_PYR_Ortho (const char eval_type)
{
	struct Basis_Data_PYR_Ortho* b_data = calloc(1,sizeof *b_data); // returned

	const struct const_Nodes* d3_p4_GLL = constructor_const_Nodes_pyr(3,4,NODES_GLL); // destructed
	const int super_type = ST_PYR;
	const int p_mass = 4;
	const char* basis_name = "pyr_ortho";
	if (eval_type == 'a') {
		b_data->phi32        = constructor_basis_pyr_orthonormal_def(2,d3_p4_GLL->rst);      // keep
		b_data->grad_phi32   = constructor_grad_basis_pyr_orthonormal_def(2,d3_p4_GLL->rst); // keep
		b_data->m_34         = constructor_mass_orthonormal_def(3,p_mass,super_type);        // keep
		b_data->grad_coef_36 = constructor_grad_vals_computation_def(3,6,basis_name);        // keep
	} else if (eval_type == 'c') {
		b_data->phi32        = constructor_basis_pyr_orthonormal(2,d3_p4_GLL->rst);      // keep
		b_data->grad_phi32   = constructor_grad_basis_pyr_orthonormal(2,d3_p4_GLL->rst); // keep
		b_data->m_34         = constructor_mass_orthonormal(3,p_mass,super_type);        // keep
		b_data->grad_coef_36 = constructor_grad_vals_computation(3,6,basis_name);        // keep
	} else {
		EXIT_UNSUPPORTED;
	}
	destructor_const_Nodes(d3_p4_GLL);

	return b_data;
}

static void destructor_Basis_Data_PYR_Ortho (struct Basis_Data_PYR_Ortho* b_data)
{
	destructor_const_Matrix_d(b_data->phi32);
	destructor_const_Multiarray_Matrix_d(b_data->grad_phi32);
	destructor_const_Matrix_d(b_data->m_34);
	destructor_const_Multiarray_d(b_data->grad_coef_36);
	free(b_data);
}

// Tensor-Product Bezier********************************************************************************************* //

static struct Basis_Data_TP_Bezier* constructor_Basis_Data_TP_Bezier (const char eval_type)
{
	struct Basis_Data_TP_Bezier* b_data = calloc(1,sizeof *b_data); // returned

	const struct const_Nodes* d1_p4_GLL = constructor_const_Nodes_tp(1,4,NODES_GLL), // destructed
	                        * d2_p4_GLL = constructor_const_Nodes_tp(2,4,NODES_GLL), // destructed
	                        * d3_p4_GLL = constructor_const_Nodes_tp(3,4,NODES_GLL); // destructed
	const char* basis_name = "tp_bezier";
	if (eval_type == 'a') {
		b_data->phi13        = constructor_basis_tp_bezier_def(3,d1_p4_GLL->rst);      // keep
		b_data->phi22        = constructor_basis_tp_bezier_def(2,d2_p4_GLL->rst);      // keep
		b_data->phi32        = constructor_basis_tp_bezier_def(2,d3_p4_GLL->rst);      // keep
		b_data->grad_phi13   = constructor_grad_basis_tp_bezier_def(3,d1_p4_GLL->rst); // keep
		b_data->grad_phi22   = constructor_grad_basis_tp_bezier_def(2,d2_p4_GLL->rst); // keep
		b_data->grad_phi31   = constructor_grad_basis_tp_bezier_def(1,d3_p4_GLL->rst); // keep
		b_data->p_14         = constructor_part_unity_def(d1_p4_GLL->rst->ext_0);      // keep
		b_data->p_24         = constructor_part_unity_def(d2_p4_GLL->rst->ext_0);      // keep
		b_data->p_34         = constructor_part_unity_def(d3_p4_GLL->rst->ext_0);      // keep
		b_data->grad_coef_13 = constructor_grad_vals_computation_def(1,3,basis_name);  // keep
		b_data->grad_coef_23 = constructor_grad_vals_computation_def(2,3,basis_name);  // keep
		b_data->grad_coef_33 = constructor_grad_vals_computation_def(3,3,basis_name);  // keep
	} else if (eval_type == 'c') {
		b_data->phi13        = constructor_basis_tp_bezier(3,d1_p4_GLL->rst);      // keep
		b_data->phi22        = constructor_basis_tp_bezier(2,d2_p4_GLL->rst);      // keep
		b_data->phi32        = constructor_basis_tp_bezier(2,d3_p4_GLL->rst);      // keep
		b_data->grad_phi13   = constructor_grad_basis_tp_bezier(3,d1_p4_GLL->rst); // keep
		b_data->grad_phi22   = constructor_grad_basis_tp_bezier(2,d2_p4_GLL->rst); // keep
		b_data->grad_phi31   = constructor_grad_basis_tp_bezier(1,d3_p4_GLL->rst); // keep
		b_data->p_14         = constructor_part_unity(b_data->phi13);              // keep
		b_data->p_24         = constructor_part_unity(b_data->phi22);              // keep
		b_data->p_34         = constructor_part_unity(b_data->phi32);              // keep
		b_data->grad_coef_13 = constructor_grad_vals_computation(1,3,basis_name);  // keep
		b_data->grad_coef_23 = constructor_grad_vals_computation(2,3,basis_name);  // keep
		b_data->grad_coef_33 = constructor_grad_vals_computation(3,3,basis_name);  // keep
	} else {
		EXIT_UNSUPPORTED;
	}
	destructor_const_Nodes(d1_p4_GLL);
	destructor_const_Nodes(d2_p4_GLL);
	destructor_const_Nodes(d3_p4_GLL);

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
	destructor_const_Multiarray_d(b_data->grad_coef_13);
	destructor_const_Multiarray_d(b_data->grad_coef_23);
	destructor_const_Multiarray_d(b_data->grad_coef_33);
	free(b_data);
}

// Simplex Bezier *************************************************************************************************** //

static struct Basis_Data_SI_Bezier* constructor_Basis_Data_SI_Bezier (const char eval_type)
{
	struct Basis_Data_SI_Bezier* b_data = calloc(1,sizeof *b_data); // returned

	const struct const_Nodes* d2_p4_AO = constructor_const_Nodes_si(2,4,NODES_AO), // destructed
	                        * d3_p4_AO = constructor_const_Nodes_si(3,4,NODES_AO); // destructed
	const char* basis_name = "si_bezier";
	if (eval_type == 'a') {
		b_data->phi22        = constructor_basis_si_bezier_def(2,d2_p4_AO->rst);      // keep
		b_data->phi23        = constructor_basis_si_bezier_def(3,d2_p4_AO->rst);      // keep
		b_data->phi32        = constructor_basis_si_bezier_def(2,d3_p4_AO->rst);      // keep
		b_data->grad_phi22   = constructor_grad_basis_si_bezier_def(2,d2_p4_AO->rst); // keep
		b_data->grad_phi31   = constructor_grad_basis_si_bezier_def(1,d3_p4_AO->rst); // keep
		b_data->p_24         = constructor_part_unity_def(d2_p4_AO->rst->ext_0);      // keep
		b_data->p_34         = constructor_part_unity_def(d3_p4_AO->rst->ext_0);      // keep
		b_data->grad_coef_25 = constructor_grad_vals_computation_def(2,5,basis_name); // keep
		b_data->grad_coef_37 = constructor_grad_vals_computation_def(3,7,basis_name); // keep
	} else if (eval_type == 'c') {
		b_data->phi22        = constructor_basis_si_bezier(2,d2_p4_AO->rst);      // keep
		b_data->phi23        = constructor_basis_si_bezier(3,d2_p4_AO->rst);      // keep
		b_data->phi32        = constructor_basis_si_bezier(2,d3_p4_AO->rst);      // keep
		b_data->grad_phi22   = constructor_grad_basis_si_bezier(2,d2_p4_AO->rst); // keep
		b_data->grad_phi31   = constructor_grad_basis_si_bezier(1,d3_p4_AO->rst); // keep
		b_data->p_24         = constructor_part_unity(b_data->phi22);             // keep
		b_data->p_34         = constructor_part_unity(b_data->phi32);             // keep
		b_data->grad_coef_25 = constructor_grad_vals_computation(2,5,basis_name); // keep
		b_data->grad_coef_37 = constructor_grad_vals_computation(3,7,basis_name); // keep
	} else {
		EXIT_UNSUPPORTED;
	}
	destructor_const_Nodes(d2_p4_AO);
	destructor_const_Nodes(d3_p4_AO);

	return b_data;
}

static void destructor_Basis_Data_SI_Bezier (struct Basis_Data_SI_Bezier* b_data)
{
	destructor_const_Matrix_d(b_data->phi22);
	destructor_const_Matrix_d(b_data->phi23);
	destructor_const_Matrix_d(b_data->phi32);
	destructor_const_Multiarray_Matrix_d(b_data->grad_phi22);
	destructor_const_Multiarray_Matrix_d(b_data->grad_phi31);
	destructor_const_Vector_d(b_data->p_24);
	destructor_const_Vector_d(b_data->p_34);
	destructor_const_Multiarray_d(b_data->grad_coef_25);
	destructor_const_Multiarray_d(b_data->grad_coef_37);
	free(b_data);
}
