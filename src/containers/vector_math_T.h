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
 *  \brief Provides Vector_\* math functions.
 */

struct Vector_T;
struct Vector_R;
struct const_Vector_T;
struct const_Vector_R;

/// \brief Invert each of the entries of the input \ref Vector_T\*.
void invert_Vector_T
	(struct Vector_T* a ///< Input vector.
	);

/// \brief Add the input value to all entries of the \ref Vector_T\*.
void add_val_to_Vector_T
	(struct Vector_T*const a, ///< To be added to.
	 const Type val           ///< Value to add.
	);

/// \brief Add to a \ref Vector_T\*.
void add_to_Vector_T_T
	(struct Vector_T* a, ///< To be added to.
	 const Type* b     ///< Data to add.
	);

/// \brief Add a \ref Vector_T\* to a \ref Vector_T\*.
void add_to_Vector_T
	(struct Vector_T* a,            ///< To be added to.
	 const struct const_Vector_T* b ///< To be added.
	);

/** \brief Return the dot product of the two input vectors, optionally scaled by the scaling constant.
 *  \return See brief. */
Type dot_product_Vector_T
	(const Type alpha,                    ///< Scaling constant.
	 const struct const_Vector_T*const a, ///< The 1st input.
	 const struct const_Vector_T*const b  ///< The 2nd input.
	);

/** \brief Return the dot product of the two input vectors, optionally scaled by the scaling constant ('R'eal, 'T'ype).
 *  \return See brief. */
Type dot_product_Vector_RT
	(const Type alpha,                    ///< Scaling constant.
	 const struct const_Vector_R*const a, ///< The 1st input.
	 const struct const_Vector_T*const b  ///< The 2nd input.
	);

/** \brief Store the element-wise multiplication of the input \ref const_Vector_T\* entries in the output \ref
 *         Vector_T\*. */
void dot_mult_Vector_T
	(const Type alpha,                    ///< Scaling constant.
	 const struct const_Vector_T*const a, ///< 1st input vector.
	 const struct const_Vector_T*const b, ///< 2nd input vector.
	 struct Vector_T*const c              ///< output vector.
	);
