// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__bases_h__INCLUDED
#define DPG__bases_h__INCLUDED
/** \file
 *  \brief Provides the interface to functions evaluating the basis functions at the input reference coordinates.
 */

#include <stddef.h>

// Constructor functions ******************************************************************************************** //

/** \brief Constructor for a polynomial operator for the tensor-product orthonomal basis.
 *  \return Standard.
 *
 *  The dimension of the basis functions is determined according to the dimension of the input nodes.
 */
const struct const_Matrix_d* constructor_basis_tp_orthonormal
	(const int p_b,                        ///< The order of the basis.
	 const struct const_Matrix_d*const rst ///< The nodes at which the basis functions are evaluated.
	);

/** \brief Constructor for a polynomial operator for the simplex orthonomal basis.
 *  \return Standard.
 *
 *  The dimension of the basis functions is determined according to the dimension of the input nodes.
 */
const struct const_Matrix_d* constructor_basis_si_orthonormal
	(const int p_b,                        ///< The order of the basis.
	 const struct const_Matrix_d*const rst ///< The nodes at which the basis functions are evaluated.
	);

// Helper functions ************************************************************************************************* //

/** \brief Compute the number of basis functions for the given element type.
 *  \return See brief. */
ptrdiff_t compute_n_basis
	(const int d,         ///< The dimension.
	 const int p_b,       ///< The basis order.
	 const int super_type ///< Defined in \ref definitions_elements.h.
	);

/** \brief Constructor for abc coordinates from rst coordinates using a Duffy-type transform.
 *  \return Standard.
 *
 *  <!-- References: -->
 *  \todo Add reference.
 */
const struct const_Matrix_d* constructor_abc_from_rst_si
	(const struct const_Matrix_d*const rst ///< The input coordinates on the reference simplex.
	);

#endif // DPG__bases_h__INCLUDED
