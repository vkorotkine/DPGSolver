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

#ifndef DPG__operator_h__INCLUDED
#define DPG__operator_h__INCLUDED
/** \file
 *  \brief Provides the real \ref Operator\* container and related functions.
 */

#include <stddef.h>
#include "definitions_core.h"

struct Operator;

#include "def_templates_type_d.h"
#include "operator_T.h"
#include "undef_templates_type.h"

#include "def_templates_type_dc.h"
#include "operator_T.h"
#include "undef_templates_type.h"

struct const_Multiarray_d;
struct Multiarray_d;
struct const_Multiarray_Matrix_c;

/// Container holding the operator matrices in various formats.
struct Operator {
	const struct const_Matrix_d*const op_std;             ///< The standard dense matrix operator.
	const struct const_Multiarray_Matrix_d*const  ops_tp; ///< The multiarray of tensor-product sub-operators.
	const struct const_Matrix_CSR_d*const op_csr;         ///< The sparse matrix operator in CSR format.
};

/// `mutable` version of \ref Operator
struct mutable_Operator {
	struct Matrix_d* op_std;            ///< The standard dense matrix operator.
	struct Multiarray_Matrix_d* ops_tp; ///< The multiarray of tensor-product sub-operators.
	struct Matrix_CSR_d* op_csr;        ///< The sparse matrix operator in CSR format.
};

// Templated functions ********************************************************************************************** //

// Interface functions ********************************************************************************************** //
// Destructors ****************************************************************************************************** //

/// \brief Destructor for a \ref Operator.
void destructor_Operator
	(const struct Operator* op ///< Standard.
	);

/// \brief Destructor for a \ref mutable_Operator.
void destructor_mutable_Operator
	(struct mutable_Operator* op ///< Standard.
	);

// Printing functions *********************************************************************************************** //

/// \brief Print a \ref Operator\* to the terminal displaying entries below the default tolerance as 0.0.
void print_Operator
	(const struct Operator*const a ///< Standard.
	);

/// \brief Print a \ref Operator\* to the terminal displaying entries below the tolerance as 0.0.
void print_Operator_tol
	(const struct Operator*const a, ///< Standard.
	 const double tol               ///< The tolerance.
	);

#endif // DPG__operator_h__INCLUDED
