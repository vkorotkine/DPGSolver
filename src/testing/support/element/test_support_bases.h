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

#ifndef DPG__test_support_bases_h__INCLUDED
#define DPG__test_support_bases_h__INCLUDED
/** \file
 *  \brief Provides support functions for testing relating to the bases.
 */

#include <stddef.h>

// Constructor functions ******************************************************************************************** //

/** \brief Constructor for a polynomial operator for the tensor-product orthonomal basis from the basis function
 *         definitions.
 *  \return Standard. */
const struct const_Matrix_d* constructor_basis_tp_orthonormal_def
	(const int p_b,                        ///< Defined in \ref constructor_basis_fptr.
	 const struct const_Matrix_d*const rst ///< Defined in \ref constructor_basis_fptr.
	);

/** \brief Constructor for polynomial operator(s) for the gradient(s) of the tensor-product orthonomal basis from the
 *         basis function definitions.
 *  \return Standard. */
const struct const_Multiarray_Matrix_d* constructor_grad_basis_tp_orthonormal_def
	(const int p_b,                        ///< Defined in \ref constructor_basis_fptr.
	 const struct const_Matrix_d*const rst ///< Defined in \ref constructor_basis_fptr.
	);

/** \brief Constructor for a polynomial operator for the simplex orthonomal basis from the basis function definitions.
 *  \return Standard. */
const struct const_Matrix_d* constructor_basis_si_orthonormal_def
	(const int p_b,                        ///< Defined in \ref constructor_basis_fptr.
	 const struct const_Matrix_d*const rst ///< Defined in \ref constructor_basis_fptr.
	);

/** \brief Constructor for polynomial operators for the gradients of the simplex orthonomal basis from the basis
 *         function definitions.
 *  \return Standard. */
const struct const_Multiarray_Matrix_d* constructor_grad_basis_si_orthonormal_def
	(const int p_b,                        ///< Defined in \ref constructor_basis_fptr.
	 const struct const_Matrix_d*const rst ///< Defined in \ref constructor_basis_fptr.
	);

/** \brief Constructor for a polynomial operator for the pyramid orthonomal basis from the basis function definitions.
 *  \return Standard. */
const struct const_Matrix_d* constructor_basis_pyr_orthonormal_def
	(const int p_b,                        ///< Defined in \ref constructor_basis_fptr.
	 const struct const_Matrix_d*const rst ///< Defined in \ref constructor_basis_fptr.
	);

/** \brief Constructor for polynomial operators for the gradients of the pyramid orthonomal basis from the basis
 *         function definitions.
 *  \return Standard. */
const struct const_Multiarray_Matrix_d* constructor_grad_basis_pyr_orthonormal_def
	(const int p_b,                        ///< Defined in \ref constructor_basis_fptr.
	 const struct const_Matrix_d*const rst ///< Defined in \ref constructor_basis_fptr.
	);

/** \brief Constructor for a polynomial operator for the tensor-product bezier basis from the basis function
 *         definitions.
 *  \return Standard. */
const struct const_Matrix_d* constructor_basis_tp_bezier_def
	(const int p_b,                        ///< Defined in \ref constructor_basis_fptr.
	 const struct const_Matrix_d*const rst ///< Defined in \ref constructor_basis_fptr.
	);

/** \brief Constructor for polynomial operator(s) for the gradient(s) of the tensor-product bezier basis from the basis
 *         function definitions.
 *  \return Standard. */
const struct const_Multiarray_Matrix_d* constructor_grad_basis_tp_bezier_def
	(const int p_b,                        ///< Defined in \ref constructor_basis_fptr.
	 const struct const_Matrix_d*const rst ///< Defined in \ref constructor_basis_fptr.
	);

// Additional functions ********************************************************************************************* //

/** \brief Constructor for the mass matrix of an orthonormal basis from the definition (i.e. the identity matrix).
 *  \return Standard. */
const struct const_Matrix_d* constructor_mass_orthonormal_def
	(const int d,         ///< The dimension.
	 const int p_b,       ///< The basis order.
	 const int super_type ///< Defined in \ref definitions_elements.h.
	);

/** \brief Constructor for the mass matrix of an orthonormal basis using cubature.
 *  \return Standard. */
const struct const_Matrix_d* constructor_mass_orthonormal
	(const int d,         ///< The dimension.
	 const int p_b,       ///< The basis order.
	 const int super_type ///< Defined in \ref definitions_elements.h.
	);

/** \brief Constructor for basis sums from the definition (i.e. the vector with all entries = 1.0).
 *  \return Standard. */
const struct const_Vector_d* constructor_part_unity_def
	(const ptrdiff_t n_n ///< The number of nodes at which the basis functions were evaluated.
	);

/** \brief Constructor for basis sums from the sum of the basis functions at each node.
 *  \return Standard. */
const struct const_Vector_d* constructor_part_unity
	(const struct const_Matrix_d* phi_rst ///< The basis functions evaluated at the rst coordinates.
	);

/** \brief Constructor for values of the gradient of a specified polynomial from the definition.
 *  \return Standard. */
const struct const_Multiarray_d* constructor_grad_vals_computation_def
	(const int d,                ///< The dimension.
	 const int p_b,              ///< The basis order.
	 const char*const basis_name ///< The name of the basis to use.
	);

/** \brief Constructor for values of the gradient of a specified polynomial using the derivative operator.
 *  \return Standard. */
const struct const_Multiarray_d* constructor_grad_vals_computation
	(const int d,                ///< The dimension.
	 const int p_b,              ///< The basis order.
	 const char*const basis_name ///< The name of the basis to use.
	);

#endif // DPG__test_support_bases_h__INCLUDED
