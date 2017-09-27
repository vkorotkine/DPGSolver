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
/// \file

#include "element_operators_tp.h"

#include <assert.h>

#include "macros.h"
#include "definitions_core.h"
#include "definitions_element_operators.h"

#include "multiarray.h"
#include "matrix.h"
#include "vector.h"

#include "simulation.h"
#include "element.h"
#include "element_operators.h"

// Static function declarations ************************************************************************************* //

/// Container for information needed to set the sub-operators.
struct Sub_Operator_Info {
	const int d;      ///< The dimension of the input.
	const int s_type; ///< The \ref Element::s_type for which the operator is being constructed.

	const struct const_Matrix_d* sub_ops_Md[DMAX]; ///< The array of matrix sub-operators.
};

/** \brief Constructor for a \ref const_Matrix_d\* of sub-operators used for the construction of the associated
 *         tensor-product operator.
 *  \return See brief. */
static const struct const_Matrix_d* constructor_operator_tp
	(const struct Sub_Operator_Info* sub_op_info ///< \ref Sub_Operator_Info.
	);

/// \brief Set the pointers \in ref Sub_Operator_Info::sub_ops_Md to the appropriate sub-operator matrices.
static void set_sub_operator_info
	(struct Sub_Operator_Info* sub_op_info, ///< \ref Sub_Operator_Info.
	 const struct Operator_Info* op_info,   ///< \ref Sub_Operator_Info.
	 const struct Operators_TP* ops_tp,     ///< \ref Operators_TP.
	 const int ind_values                   ///< The index of the row of \ref Operator_Info::values_op to use.
	);

// Interface functions ********************************************************************************************** //

void set_operators_tp
	(struct Operators_TP* ops_tp,
	 const struct const_Multiarray_Matrix_d* op_00, const struct const_Multiarray_Matrix_d* op_01,
	 const struct const_Multiarray_Matrix_d* op_10, const struct const_Multiarray_Matrix_d* op_11
	)
{
	assert(op_00 != NULL);
	assert(op_01 != NULL);

	ops_tp->op[0][0] = op_00;
	ops_tp->op[0][1] = op_01;
	if (op_10 != NULL) {
		assert(op_11 != NULL);

		ops_tp->op[1][0] = op_10;
		ops_tp->op[1][1] = op_11;
	} else {
		assert(op_11 == NULL);

		ops_tp->op[1][0] = op_00;
		ops_tp->op[1][1] = op_01;
	}
}

const struct const_Multiarray_Matrix_d* constructor_operators_tp
	(const char*const name_type, const char*const name_in, const char*const name_out, const char*const name_range,
	 const struct const_Element* element, const struct Simulation* sim, const struct Operators_TP* ops_tp)
{
	const int p_ref[2] = { -1, -1 };
	struct Operator_Info* op_info =
		constructor_Operator_Info(name_type,name_in,name_out,name_range,p_ref,element); // destructed

	const struct Op_IO* op_io = op_info->op_io;

	const char ce_i = op_io[OP_IND_I].ce;

	assert(ce_i == 'v'); // May potentially be made flexible in future.

printf("ex_op_tp\n");
print_const_Vector_i(op_info->extents_op);
	const struct const_Multiarray_Matrix_d* op =
		constructor_empty_const_Multiarray_Matrix_d_V(false,op_info->extents_op); // returned
print_const_Matrix_i(op_info->values_op);

	const ptrdiff_t row_max = op_info->values_op->ext_0;
	for (ptrdiff_t row = 0; row < row_max; ++row) {
		const char ce_o = op_io[OP_IND_O].ce;

		struct Sub_Operator_Info sub_op_info =
			{ .d          = element->d,
			  .s_type     = element->s_type,
			  .sub_ops_Md = { NULL }, };

		set_sub_operator_info(&sub_op_info,op_info,ops_tp,row);
		switch (ce_o) {
		case 'v':

			break;
		default:
			EXIT_ERROR("Unsupported: %c.\n",ce_o);
			break;
		}

		const struct const_Matrix_d* op_ioN = constructor_operator_tp(&sub_op_info); // moved

		const int* op_values = get_row_const_Matrix_i(row,op_info->values_op);

		const int order_op = op_info->extents_op->ext_0;
		const struct const_Vector_i* indices_op =
			constructor_indices_Vector_i(order_op,op_values,NULL); // destructed
		const ptrdiff_t ind_op = compute_index_sub_container_pi(op->order,0,op->extents,indices_op->data);

		const_constructor_move_const_Matrix_d(&op->data[ind_op],op_ioN); // keep
		destructor_const_Vector_i(indices_op);
	}

UNUSED(sim);

	destructor_Operator_Info(op_info);

	return op;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Set the input `op_MMd[3]` to point to the appropriate operators from `ops_tp`.
 *  Only one of `op_val_d` and `op_val_ce` may be valid. */
static void set_op_MMd
	(const struct const_Multiarray_Matrix_d* op_MMd[3], ///< The array of \ref Multiarray_Matrix_d operators.
	 const struct Operators_TP* ops_tp,                 ///< \ref Operators_TP.
	 const int op_val_d,                                ///< The differentiation index of the current operator.
	 const int op_val_ce,                               ///< The computational element index of the current operator.
	 const int e_type                                   ///< \ref Element::type.
	);

static void set_sub_operator_info
	(struct Sub_Operator_Info* sub_op_info, const struct Operator_Info* op_info, const struct Operators_TP* ops_tp,
	 const int ind_values)
{
	const struct Op_IO* op_io = op_info->op_io;
	assert(op_io[OP_IND_I].ce == 'v');

	const int* op_values = get_row_const_Matrix_i(ind_values,op_info->values_op);
	if (check_op_info_loss(op_values)
		EXIT_ERROR("Ensure that all is working as expected.\n");

	// See comments in \ref element_operators_tp.h
	assert((op_values[OP_IND_D] == OP_INVALID_IND) ||
	       (op_values[OP_IND_CE] == OP_INVALID_IND && op_values[OP_IND_CE+1] == OP_INVALID_IND));

	const int op_val_d  = op_values[OP_IND_D];
	const int op_val_ce = op_values[OP_IND_CE+1];
	const int e_type    = op_info->element->type;

	const struct const_Multiarray_Matrix_d* op_MMd[] = { NULL, NULL, NULL };
	set_op_MMd(op_MMd,ops_tp,op_val_d,op_val_ce,e_type);



EXIT_UNSUPPORTED;

// Make a vector_i of indices:
// Need an index for the computational element index (related to ce_o)
// Need an index for which h_io indices to use.
// Need an index for which p_io indices to use.

	const int order_op = op_info->extents_op->ext_0;
	const struct const_Vector_i* indices_op = constructor_indices_Vector_i(order_op,op_values,NULL); // destructed
	const ptrdiff_t ind_op_0 = compute_index_sub_container_pi(op_MMd[0]->order,0,op_MMd[0]->extents,indices_op->data);
UNUSED(ind_op_0);
	destructor_const_Vector_i(indices_op);
UNUSED(sub_op_info);
}

static const struct const_Matrix_d* constructor_operator_tp (const struct Sub_Operator_Info* sub_op_info)
{
	// Think about how this might be done by applying the single sum-factorized operators to some form of transposed
	// input operator instead of using sparse-dense with potential permutation.
UNUSED(sub_op_info);
return NULL;
}

// Level 1 ********************************************************************************************************** //

/** \brief Set the indicator in the appropriate index for whether the second sub-element operator should be used for a
 *         differentiation related operation. */
static void set_ind_spec_diff
	(bool* ind_spec,    ///< The array of indices for the special operators.
	 const int s_type,  ///< \ref Element::s_type.
	 const int op_val_d ///< The differentiation index of the current operator.
	);

/** \brief Set the indicator in the appropriate index for whether the second sub-element operator should be used for a
 *         face/edge computational element related operation. */
static void set_ind_spec_ce
	(bool* ind_spec,     ///< The array of indices for the special operators.
	 const int e_type,   ///< \ref Element::type.
	 const int op_val_ce ///< The computational element index of the current operator.
	);

static void set_op_MMd
	(const struct const_Multiarray_Matrix_d* op_MMd[DMAX], const struct Operators_TP* ops_tp, const int op_val_d,
	 const int op_val_ce, const int e_type)
{
	const int s_type = compute_super_from_elem_type(e_type);

	assert(DMAX == 3);
	assert((op_val_d == OP_INVALID_IND) || (op_val_ce == OP_INVALID_IND));
	assert((s_type == ST_TP) || (s_type == ST_WEDGE));

	int* ind_ops_tp = ( s_type == ST_TP ? (int[]) { 0, 0, 0, } : (int[]) { 0, OP_INVALID_IND, 1, } );

	bool ind_spec[] = { false, false, false, };
	set_ind_spec_d(ind_spec,s_type,op_val_d);

	for (int i = 0; i < DMAX; ++i) {
		const int ind_op_tp = ind_ops_tp[i];
		if (ind_op_tp != OP_INVALID_IND) {
			if (!ind_spec[i])
				op_MMd[i] = ops_tp->op[ind_op_tp][0];
			else {
				op_MMd[i] = ops_tp->op[ind_op_tp][1];
			}
		} else {
			op_MMd[i] = NULL;
		}
	}
}

// Level 1 ********************************************************************************************************** //

static void set_ind_spec_diff (bool* ind_spec, const int s_type, const int op_val_d)
{
	switch (s_type) {
	case ST_TP:
		if (op_val_d != OP_INVALID_IND) {
			assert(op_val_d >= 0);
			if (op_val_d <= 2)
				ind_spec[op_val_d] = true;
			else
				EXIT_ERROR("Unsupported: %d\n",op_val_d);
		}
		break;
	case ST_WEDGE:
		if (op_val_d != OP_INVALID_IND) {
			assert(op_val_d >= 0);
			if (op_val_d < 2)
				ind_spec[0] = true;
			else if (op_val_d == 2)
				ind_spec[2] = true;
			else
				EXIT_ERROR("Unsupported: %d\n",op_val_d);
		}
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",s_type);
		break;
	}
}

static void set_ind_spec_ce (bool* ind_spec, const int e_type, const int op_val_ce)
{
	bool* ind_spec_l = NULL;
	int n_spec = -1;

/// \todo make external functions here?
	switch (e_type) {
	case QUAD:
		n_spec = 2;
		switch (ce_o) {
		case 'v':
			ind_spec_l = (bool[]) { false, false, };
			break;
		case 'f': // fallthrough
		case 'e':
			switch (op_val_ce) {
			case 0:         // vertices: {0,6}.
			case 1:         // vertices: {2,8}.
			case 4: case 5: // vertices: {0,3}, {3,6}.
			case 6: case 7: // vertices: {2,5}, {5,8}.
				ind_spec_l = (bool[]) { true, false, };
				break;
			case 2:           // face(s) with vertices: {0,2}.
			case 3:           // face(s) with vertices: {6,8}.
			case 8:  case 9:  // face(s) with vertices: {0,1}, {1,2}.
			case 10: case 11: // face(s) with vertices: {6,7}, {7,8}.
				ind_spec_l = (bool[]) { false, true, };
				break;
			default:
				EXIT_ERROR("Unsupported: %d.\n",op_val_ce);
				break;
			}
			break;
		default:
			EXIT_ERROR("Unsupported: %c.\n",ce_o);
			break;
		}
		break;
	case HEX:
		n_spec = 3;
		switch (ce_o) {
		case 'v':
			ind_spec_l = (bool[]) { false, false, false, };
			break;
		case 'f':
			switch (op_val_ce) {
			case 0:                             // { 0, 6,18,24}.
			case 1:                             // { 2, 8,20,26}.
			case 6:  case 7:  case 8:  case 9:  // { 0, 3, 9,12}, { 3, 6,12,15}, { 9,12,18,21}, {12,15,21,24}.
			case 10: case 11: case 12: case 13: // { 2, 5,11,14}, { 5, 8,14,17}, {11,14,20,23}, {14,17,23,26}.
				ind_spec_l = (bool[]) { true, false, false, };
				break;
// Continue comment change below
EXIT_UNSUPPORTED;
			case 2:                             // bottom face (complete).
			case 3:                             // top    face (complete).
			case 14: case 15: case 16: case 17: // bottom face (refined - isotropic).
			case 18: case 19: case 20: case 21: // top    face (refined - isotropic).
				ind_spec_l = (bool[]) { false, true, false, };
				break;
			case 4:                             // back  face (complete).
			case 5:                             // front face (complete).
			case 22: case 23: case 24: case 25: // back  face (refined - isotropic).
			case 26: case 27: case 28: case 29: // front face (refined - isotropic).
				ind_spec_l = (bool[]) { false, false, true, };
				break;
			default:
				EXIT_ERROR("Unsupported: %d.\n",op_val_ce);
				break;
			}
			break;
		default:
			EXIT_ERROR("Unsupported: %c.\n",ce_o);
			break;
		}
		break;
	case WEDGE:
		n_spec = 3;
		switch (ce_o) {
		case 'v':
			ind_spec_l = (bool[]) { false, false, false, };
			break;
		case 'f':
			switch (op_val_ce) {
			case 0:                             // face opposite vertices {0,3} (complete).
			case 1:                             // face opposite vertices {1,4} (complete).
			case 2:                             // face opposite vertices {2,5} (complete).
			case 5:  case 6:  case 7:  case 8:  // face opposite vertices {0,3} (refined - isotropic).
			case 9:  case 10: case 11: case 12: // face opposite vertices {1,4} (refined - isotropic).
			case 13: case 14: case 15: case 16: // face opposite vertices {2,5} (refined - isotropic).
				ind_spec_l = (bool[]) { true, false, false, };
				break;
			case 3:                             // face with vertices {0,1,2} (complete).
			case 4:                             // face with vertices {3,4,5} (complete).
			case 17: case 18: case 19: case 20: // face with vertices {0,1,2} (refined - isotropic).
			case 21: case 22: case 23: case 24: // face with vertices {3,4,5} (refined - isotropic).
				ind_spec_l = (bool[]) { false, false, true, };
				break;
			default:
				EXIT_ERROR("Unsupported: %d.\n",op_val_ce);
				break;
			}
			break;
		default:
			EXIT_ERROR("Unsupported: %c.\n",ce_o);
			break;
		}
		break;
	default:
		EXIT_ERROR("Unsupported: %d.\n",e_type);
		break;
	}

	for (int i = 0; i < n_spec; ++i)
		ind_spec[i] = ind_spec_l[i];
}
