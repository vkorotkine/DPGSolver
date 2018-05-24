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

#ifndef DPG__bases_h__INCLUDED
#define DPG__bases_h__INCLUDED
/** \file
 *  \brief Provides the interface to functions evaluating the basis functions at the input reference coordinates.
 */

#include <stddef.h>

/// Container for the set of basis functions used to represent the geometry or the solution
struct Basis {
	int basis_type;	          ///< The basis type as defined in \ref definition_bases
	int order;                ///< The order of the basis.
	int n_functions;          ///< The number of basis functions in the 1D case. For polynomials, n_bf = order + 1.

	int n_knot;               ///< The number of knots in the knot vector. Note: n_bf = n_k - (p_b+1)
	struct Vector_d* knots;   ///< Knot vector used to form B-splines and NURBS
	struct Vector_d* weights; ///< Weights associated with each basis to form NURBS. Weights of 1 result in B-spline.
};

/** \brief Function pointer to basis constructor function.
 *  \param p_b The order of the basis.
 *  \param rst The nodes at which the basis functions are evaluated.
 *
 *  The dimension of the basis functions is determined according to the dimension of the input nodes.
 */
typedef const struct const_Matrix_d* (*constructor_basis_fptr)
	(const struct Basis*const basis,
	 const struct const_Matrix_d*const rst
	);

/** \brief Function pointer to basis gradient constructor function.
 *  \param p_b Defined in \ref constructor_basis_fptr.
 *  \param rst Defined in \ref constructor_basis_fptr.
 */
typedef const struct const_Multiarray_Matrix_d* (*constructor_grad_basis_fptr)
	(const struct Basis*const basis,
	 const struct const_Matrix_d*const rst
	);

// Constructor functions ******************************************************************************************** //

/** \brief Constructor for a polynomial operator for a general basis of the given super type and basis type.
 *  \return Standard. */
const struct const_Matrix_d* constructor_basis
	(const struct Basis*const basis,        ///< Defined in \ref constructor_basis_fptr.
	 const struct const_Matrix_d*const rst, ///< Defined in \ref constructor_basis_fptr.
	 const int s_type,                      ///< \ref Element::s_type.
	 const char*const basis_type            ///< The basis type.
	);

/** \brief Constructor for polynomial operator(s) for the gradient(s) of the general basis of the given super type and
 *         basis type.
 *  \return Standard. */
const struct const_Multiarray_Matrix_d* constructor_grad_basis
	(const struct Basis*const basis,        ///< Defined in \ref constructor_basis_fptr.
	 const struct const_Matrix_d*const rst, ///< Defined in \ref constructor_basis_fptr.
	 const int s_type,                      ///< \ref Element::s_type.
	 const char*const basis_type            ///< The basis type.
	);

// Orthonormal basis ************************************************************************************************ //

/** \brief Constructor for a polynomial operator for the tensor-product orthonomal basis (normalized Legendre
 *         polynomials).
 *  \return Standard.
 *
 *  The dimension of the basis functions is determined according to the dimension of the input nodes.
 */
const struct const_Matrix_d* constructor_basis_tp_orthonormal
	(const struct Basis*const basis,       ///< Defined in \ref constructor_basis_fptr.
	 const struct const_Matrix_d*const rst ///< Defined in \ref constructor_basis_fptr.
	);

/** \brief Constructor for polynomial operator(s) for the gradient(s) of the tensor-product orthonomal basis
 *         (derivative of normalized Legendre polynomials).
 *  \return Standard. */
const struct const_Multiarray_Matrix_d* constructor_grad_basis_tp_orthonormal
	(const struct Basis*const basis,                        ///< Defined in \ref constructor_basis_fptr.
	 const struct const_Matrix_d*const rst ///< Defined in \ref constructor_basis_fptr.
	);

/** \brief Constructor for a polynomial operator for the simplex orthonomal basis (appendix A.1, \cite Hesthaven2007).
 *  \return Standard.
 *
 *  The dimension of the basis functions is determined according to the dimension of the input nodes.
 */
const struct const_Matrix_d* constructor_basis_si_orthonormal
	(const struct Basis*const basis,                        ///< Defined in \ref constructor_basis_fptr.
	 const struct const_Matrix_d*const rst ///< Defined in \ref constructor_basis_fptr.
	);

/** \brief Constructor for polynomial operators for the gradients of the simplex orthonomal basis.
 *  \return Standard. */
const struct const_Multiarray_Matrix_d* constructor_grad_basis_si_orthonormal
	(const struct Basis*const basis,                        ///< Defined in \ref constructor_basis_fptr.
	 const struct const_Matrix_d*const rst ///< Defined in \ref constructor_basis_fptr.
	);

/** \brief Constructor for a polynomial operator for the pyramid orthonomal basis (eq. (2.1), \cite Chan2016).
 *  \return Standard.
 *
 *  The dimension of the basis functions is determined according to the dimension of the input nodes.
 */
const struct const_Matrix_d* constructor_basis_pyr_orthonormal
	(const struct Basis*const basis,                        ///< Defined in \ref constructor_basis_fptr.
	 const struct const_Matrix_d*const rst ///< Defined in \ref constructor_basis_fptr.
	);

/** \brief Constructor for polynomial operators for the gradients of the pyramid orthonomal basis.
 *  \return Standard. */
const struct const_Multiarray_Matrix_d* constructor_grad_basis_pyr_orthonormal
	(const struct Basis*const basis,                        ///< Defined in \ref constructor_basis_fptr.
	 const struct const_Matrix_d*const rst ///< Defined in \ref constructor_basis_fptr.
	);

// Bezier basis ***************************************************************************************************** //

/** \brief Constructor for a polynomial operator for the tensor-product bezier basis (Bernstein polynomials, (section
 *         2.1 \cite Prautzsch2002)).
 *  \return Standard.
 *
 *  The dimension of the basis functions is determined according to the dimension of the input nodes.
 */
const struct const_Matrix_d* constructor_basis_tp_bezier
	(const struct Basis*const basis,                        ///< Defined in \ref constructor_basis_fptr.
	 const struct const_Matrix_d*const rst ///< Defined in \ref constructor_basis_fptr.
	);

/** \brief Constructor for polynomial operator(s) for the gradient(s) of the tensor-product bezier basis (derivatives of
 *         Bernstein polynomials, (section 2.4 \cite Prautzsch2002)).
 *  \return Standard. */
const struct const_Multiarray_Matrix_d* constructor_grad_basis_tp_bezier
	(const struct Basis*const basis,                        ///< Defined in \ref constructor_basis_fptr.
	 const struct const_Matrix_d*const rst ///< Defined in \ref constructor_basis_fptr.
	);

/** \brief Version of \ref constructor_basis_fptr for the simplex bezier basis (section 2, \cite Chan2016_bez).
 *  \return Standard.
 *
 *  Note that the extension of the extension to tetrahedral elements is straightforward after noting the general
 *  extension of the Duffy-type transform for simplices (eq. (3.3), \cite Ainsworth2011).
 */
const struct const_Matrix_d* constructor_basis_si_bezier
	(const struct Basis*const basis,                        ///< See brief.
	 const struct const_Matrix_d*const rst ///< See brief.
	);

/** \brief Version of \ref constructor_grad_basis_fptr for the simplex bezier basis.
 *  \return Standard. */
const struct const_Multiarray_Matrix_d* constructor_grad_basis_si_bezier
	(const struct Basis*const basis,                        ///< See brief.
	 const struct const_Matrix_d*const rst ///< See brief.
	);

/** \brief Version of \ref constructor_basis_fptr for the pyramid bezier basis (section 3, \cite Chan2016_bez).
 *  \return Standard. */
const struct const_Matrix_d* constructor_basis_pyr_bezier
	(const struct Basis*const basis,                        ///< See brief.
	 const struct const_Matrix_d*const rst ///< See brief.
	);

/** \brief Version of \ref constructor_grad_basis_fptr for the pyramid bezier basis.
 *  \return Standard. */
const struct const_Multiarray_Matrix_d* constructor_grad_basis_pyr_bezier
	(const struct Basis*const basis,                        ///< See brief.
	 const struct const_Matrix_d*const rst ///< See brief.
	);

// NURBS basis ***************************************************************************************************** //

/** \brief Constructor for a polynomial operator for the tensor-product NURBS basis (\cite Piegl1995).
 *  \return Standard.
 *
 *  The dimension of the basis functions is determined according to the dimension of the knot vector.
 */
const struct const_Matrix_d* constructor_basis_tp_nurbs
	(const struct Basis*const basis,       ///< Defined in \ref constructor_basis_fptr.
	 const struct const_Matrix_d*const rst ///< Defined in \ref constructor_basis_fptr.
	);

/** \brief Constructor for polynomial operator(s) for the gradient(s) of the tensor-product bezier basis (derivatives of
 *         Bernstein polynomials, (section 2.4 \cite Prautzsch2002)).
 *  \return Standard. */
const struct const_Multiarray_Matrix_d* constructor_grad_basis_tp_nurbs
	(const struct Basis*const basis,       ///< Defined in \ref constructor_basis_fptr.
	 const struct const_Matrix_d*const rst ///< Defined in \ref constructor_basis_fptr.
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

/** \brief Constructor for simplex barycentric coordinates from rst coordinates.
 *  \return Standard.
 *
 *  The barycentric coordinates are given by Szabo et al. (eq. (6.6), eq. (13.20), \cite Szabo1991) with a correction to
 *  the constant terms such that the reference simplices are centered at the origin.
 */
const struct const_Matrix_d* constructor_bcoords_from_rst_si
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

/** \brief Get a pointer to the appropriate basis constructor function based on the input element super type.
 *  \return See brief. */
constructor_basis_fptr get_constructor_basis_by_super_type
	(const int s_type,               ///< \ref Element::s_type.
	 const char*const ref_basis_name ///< The name of the reference basis to be used.
	);

/** \brief Version of \ref get_constructor_basis_by_super_type for the bezier basis.
 *  \return See brief. */
constructor_basis_fptr get_constructor_basis_bezier_by_super_type
	(const int s_type ///< See brief.
	);

/** \brief Same as \ref get_constructor_basis_by_super_type but with `int` input.
 *  \return See brief. */
constructor_basis_fptr get_constructor_basis_by_super_type_i
	(const int s_type,   ///< \ref Element::s_type.
	 const int ind_basis ///< The index of the basis to be used.
	);

/** \brief Get a pointer to the appropriate basis gradient constructor function based on the input element super type.
 *  \return See brief. */
constructor_grad_basis_fptr get_constructor_grad_basis_by_super_type
	(const int s_type,               ///< \ref Element::s_type.
	 const char*const ref_basis_name ///< The name of the reference basis to be used.
	);

/** \brief Get the integer index of the basis based on the string input.
 *  \return See brief. */
int get_basis_i_from_s
 	(const char*const basis_name_s ///< The `char*` name of the basis.
	);

#endif // DPG__bases_h__INCLUDED
