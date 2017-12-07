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

#include "operator.h"
#include <assert.h>

#include "macros.h"
#include "definitions_core.h"
#include "definitions_tol.h"

#include "multiarray.h"
#include "matrix.h"

// Templated functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "operator_T.c"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //
// Destructors ****************************************************************************************************** //

void destructor_Operator (const struct Operator* op)
{
	destructor_mutable_Operator((struct mutable_Operator*)op);
}

void destructor_mutable_Operator (struct mutable_Operator* op)
{
	if (op->op_std)
		destructor_Matrix_d(op->op_std);
	if (op->ops_tp) {
		assert(op->ops_tp->owns_data == false);
		destructor_Multiarray_Matrix_d(op->ops_tp);
	}
	if (op->op_csr)
		EXIT_ADD_SUPPORT;
//		destructor_Matrix_CSR_d(op->op_csr);
	free(op);
}

// Printing functions *********************************************************************************************** //

void print_Operator (const struct Operator*const a)
{
	print_Operator_tol(a,EPS);
}

void print_Operator_tol (const struct Operator*const a, const double tol)
{
	printf("\n");
	printf("%-35s","\tdense operator:");
	if (a->op_std) {
		printf("\n\n");
		print_const_Matrix_d_tol(a->op_std,tol);
	} else {
		printf("*** NULL ***\n");
	}

	printf("%-35s","\ttensor-product sub-operators:");
	if (a->ops_tp) {
		printf("\n{\n\n");
		print_const_Multiarray_Matrix_d_tol(a->ops_tp,tol);
		printf("}\n");
	} else {
		printf("*** NULL ***\n");
	}

	printf("%-35s","\tsparse (CSR) operator:");
	if (a->op_csr) {
		printf("\n\n");
		EXIT_ADD_SUPPORT;
	} else {
		printf("*** NULL ***\n");
	}
	printf("\n");
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
