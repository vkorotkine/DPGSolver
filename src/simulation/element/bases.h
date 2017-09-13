// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__bases_h__INCLUDED
#define DPG__bases_h__INCLUDED
/** \file
 *  \brief Provides the interface to functions evaluating the basis functions at the input reference coordinates.
 */

#include <stddef.h>

// Constructor functions ******************************************************************************************** //

/** \brief Function pointer to basis functions.
 *  \param p_b The order of the basis.
 *  \param rst The nodes at which the basis functions are evaluated.
 */
typedef const struct const_Matrix_d* (*basis_fptr)
	(const int p_b,
	 const struct const_Matrix_d*const rst
	);

/** \brief Constructor for a polynomial operator for the tensor-product orthonomal basis (normalized Legendre
 *         polynomials).
 *  \return Standard.
 *
 *  The dimension of the basis functions is determined according to the dimension of the input nodes.
 */
const struct const_Matrix_d* constructor_basis_tp_orthonormal
	(const int p_b,                        ///< Defined in \ref basis_fptr.
	 const struct const_Matrix_d*const rst ///< Defined in \ref basis_fptr.
	);

/** \brief Constructor for a polynomial operator for the gradient of the tensor-product orthonomal basis (derivative of
 *         normalized Legendre polynomials).
 *  \return Standard. */
const struct const_Multiarray_Matrix_d* constructor_grad_basis_tp_orthonormal
	(const int p_b,                        ///< Defined in \ref basis_fptr.
	 const struct const_Matrix_d*const rst ///< Defined in \ref basis_fptr.
	);

/** \brief Constructor for a polynomial operator for the simplex orthonomal basis (appendix A.1, \cite Hesthaven2007).
 *  \return Standard.
 *
 *  The dimension of the basis functions is determined according to the dimension of the input nodes.
 */
const struct const_Matrix_d* constructor_basis_si_orthonormal
	(const int p_b,                        ///< Defined in \ref basis_fptr.
	 const struct const_Matrix_d*const rst ///< Defined in \ref basis_fptr.
	);

/** \brief Constructor for a polynomial operator for the pyramid orthonomal basis (eq. (2.1), \cite Chan2016).
 *  \return Standard.
 *
 *  The dimension of the basis functions is determined according to the dimension of the input nodes.
 */
const struct const_Matrix_d* constructor_basis_pyr_orthonormal
	(const int p_b,                        ///< Defined in \ref basis_fptr.
	 const struct const_Matrix_d*const rst ///< Defined in \ref basis_fptr.
	);

/** \brief Constructor for a polynomial operator for the tensor-product bezier basis (Bernstein polynomials, (section
 *         2.1 \cite Prautzsch2002)).
 *  \return Standard.
 *
 *  The dimension of the basis functions is determined according to the dimension of the input nodes.
 */
const struct const_Matrix_d* constructor_basis_tp_bezier
	(const int p_b,                        ///< Defined in \ref basis_fptr.
	 const struct const_Matrix_d*const rst ///< Defined in \ref basis_fptr.
	);

// Helper functions ************************************************************************************************* //

/** \brief Compute the number of basis functions for the given element type.
 *  \return See brief. */
ptrdiff_t compute_n_basis
	(const int d,         ///< The dimension.
	 const int p_b,       ///< The basis order.
	 const int super_type ///< Defined in \ref definitions_elements.h.
	);

/** \brief Constructor for simplex abc coordinates from rst coordinates using a Duffy-type transform.
 *  \return Standard.
 *
 *  The transform is that of Hesthaven et al. (appendix A.1, \cite Hesthaven2007) using the symmetric regular simplices
 *  instead of the bi-unit right simplices.
 */
const struct const_Matrix_d* constructor_abc_from_rst_si
	(const struct const_Matrix_d*const rst ///< The input coordinates on the reference simplex.
	);

/** \brief Constructor for pyramid abc coordinates from rst coordinates using a Duffy-type transform.
 *  \return Standard.
 *
 *  The transform is that of Chan et al. (section 2, \cite Chan2016) using the symmetric regular pyramid instead of the
 *  bi-unit right pyramid.
 */
const struct const_Matrix_d* constructor_abc_from_rst_pyr
	(const struct const_Matrix_d*const rst ///< The input coordinates on the reference simplex.
	);

#endif // DPG__bases_h__INCLUDED
