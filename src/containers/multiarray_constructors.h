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

#ifndef DPG__multiarray_constructors_h__INCLUDED
#define DPG__multiarray_constructors_h__INCLUDED
/** \file
 *  \brief Provides real Multiarray_\* constructors and destructors.
 *
 *  \section s1_Multi_con Naming Convention
 *
 *  Function names are chosen according to the following template: `constructor_{0}_(1)_{2}_(3)_(4)` where elements
 *  in curly braces `{}` are required and those in round brackets `()` are optional.
 *  	- {0} : Type of constructor
 *  		- default: reserve memory only for the container itself (not for the data).
 *  		- empty:   reserve memory storage for the data but do not set it.
 *  		- zero:    same as empty, but with memory initially set to zero.
 *  		- copy:    reserve memory storage for the data and set it to a copy of the input data.
 *  		- move:    move data to the container being constructed.
 *  		- set:     construct a container and set its entries to the input container but use `owns_data = false`.
 *  	- (1) : Optional `const` specifier
 *  	- {2} : Type of container to be returned
 *  		- Multiarray_\*: d, Vector_i
 *  	- (3) : Level of dereferencing if not equal to 1.
 *  	- (4) : Type of input from which the container is constructed
 *
 *  \section s2_Multi_con Const Versions of Multiarray Containers
 *
 *  In order to avoid unintentionally overwriting data which should be constant, `const` versions of the containers
 *  are also provided where relevant. A complication which arises as a result of declaring objects `const` is that
 *  it is not possible to define them! To overcome this difficulty, casts are used to set the data.
 *
 *  \warning This procedure exhibits undefined behaviour (relating to changing a `const`-qualified type) unless the
 *  memory for the struct is **dynamically allocated**. See [this SO answer][SO_dyn_const_struct] for a detailed
 *  explanation of the problem and possible approaches used to overcome it.
 *
 *  \warning It is **required** that `const` and non-`const` versions of the containers have identical memory layout
 *           such that casts can be used to convert between them.
 *
 *  Noting that containers should be dynamically allocated, we have the implication that containers with no
 *  dereferencing should never be used. Taken together with the upper limitation on the level of dereferencing, it
 *  is thus required that **containers have exactly one level of dereferencing**. The only place in which an
 *  additional level of dereferencing is permitted is as a list of pointers to containers in a multiarray.
 *
 *  \subsection s21_Multi_con Constructors for `const` Containers
 *
 *  The `const_constructor_move_*` functions are used to *define* the `const` equivalent of the container. The
 *  functions move **all** container members through casts.
 *
 *  \note This means that a single destructor call should be made for both the `src` and `dest` variables.
 *
 *  \section s3_Multi_con Compound Literals
 *
 *  In the interest of greater generic programming, compound literals should be passed as multiarray constructor
 *  arguments when `order > 2`. Example implementations can be found in questions and answers of
 *  [this][SO_compound1] and [this][SO_compound2] SO discussion.
 *
 *  <!-- References: -->
 *  [SO_dyn_const_struct]: https://stackoverflow.com/q/2219001/5983549
 *  [SO_compound1]: https://stackoverflow.com/q/17374792/5983549
 *  [SO_compound2]: https://stackoverflow.com/a/11268460/5983549
 */

#include <stddef.h>

struct Vector_T;
struct Matrix_T;
struct Multiarray_Vector_T;
struct Multiarray_Matrix_T;
struct const_Vector_i;
struct const_Vector_T;
struct const_Matrix_T;
struct const_Multiarray_T;
struct const_Multiarray_Vector_T;
struct const_Multiarray_Matrix_T;

#include "def_templates_type_d.h"
#include "def_templates_matrix_d.h"
#include "def_templates_multiarray_d.h"
#include "def_templates_vector_d.h"
#include "multiarray_constructors_T.h"
#include "undef_templates_type.h"
#include "undef_templates_matrix.h"
#include "undef_templates_multiarray.h"
#include "undef_templates_vector.h"

#include "def_templates_type_i.h"
#include "def_templates_matrix_i.h"
#include "def_templates_multiarray_i.h"
#include "def_templates_vector_i.h"
#include "multiarray_constructors_T.h"
#include "undef_templates_type.h"
#include "undef_templates_matrix.h"
#include "undef_templates_multiarray.h"
#include "undef_templates_vector.h"

/** \brief Allocated and set the `extents` for a `Multiarray_*`.
 *  \return See brief. */
ptrdiff_t* allocate_and_set_extents
	(const int order,                ///< Defined in \ref Multiarray_d.
	 const ptrdiff_t*const extents_i ///< The input extents.
	);

#endif // DPG__multiarray_constructors_h__INCLUDED
