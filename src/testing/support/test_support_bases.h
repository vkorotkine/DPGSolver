// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__test_support_bases_h__INCLUDED
#define DPG__test_support_bases_h__INCLUDED
/** \file
 *  \brief Provides support functions for testing relating to the bases.
 */

// Constructor functions ******************************************************************************************** //

/** \brief Constructor for a polynomial operator for the tensor-product orthonomal basis from the basis function
 *         definitions.
 *  \return Standard. */
const struct const_Matrix_d* constructor_basis_tp_orthonormal_def
	(const int p_b,                        ///< The order of the basis.
	 const struct const_Matrix_d*const rst ///< The nodes at which the basis functions are evaluated.
	);

/** \brief Constructor for a polynomial operator for the simplex orthonomal basis from the basis function definitions.
 *  \return Standard. */
const struct const_Matrix_d* constructor_basis_si_orthonormal_def
	(const int p_b,                        ///< The order of the basis.
	 const struct const_Matrix_d*const rst ///< The nodes at which the basis functions are evaluated.
	);

/** \brief Constructor for a polynomial operator for the pyramid orthonomal basis from the basis function definitions.
 *  \return Standard. */
const struct const_Matrix_d* constructor_basis_pyr_orthonormal_def
	(const int p_b,                        ///< The order of the basis.
	 const struct const_Matrix_d*const rst ///< The nodes at which the basis functions are evaluated.
	);

/** \brief Constructor for a polynomial operator for the tensor-product bezier basis from the basis function
 *         definitions.
 *  \return Standard. */
const struct const_Matrix_d* constructor_basis_tp_bezier_def
	(const int p_b,                        ///< The order of the basis.
	 const struct const_Matrix_d*const rst ///< The nodes at which the basis functions are evaluated.
	);

// Orthogonality functions ****************************************************************************************** //

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

#endif // DPG__test_support_bases_h__INCLUDED
