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
#include "definitions_h_ref.h"

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

/** \brief Set the input `op_MMd[DMAX]` to point to the appropriate operators from `ops_tp`.
 *
 *  `spec_indices` can be interpreted as the direction(s) (n-Cube: r,s,t; wedge: rs,NULL,t) in which to use the
 *  'special' operator (the second operator for the sub-element).
 */
static void set_ops_MMd
	(const struct const_Multiarray_Matrix_d* op_MMd[DMAX], ///< The array of \ref Multiarray_Matrix_d operators.
	 bool* spec_indices,                                   ///< The array of indices for the special operators.
	 const struct Operators_TP* ops_tp,                    ///< \ref Operators_TP.
	 const int* op_values,                                 ///< The index values of the current operator.
	 const struct Operator_Info* op_info                   ///< \ref Sub_Operator_Info.
	);

/// \brief Set the the input `op_Md[DMAX]` to point to the appropriate matrices in `op_MMd[DMAX]`.
static void set_ops_Md
	(const struct const_Matrix_d* op_Md[DMAX],             ///< The array of \ref Matrix_d operators.
	 const struct const_Multiarray_Matrix_d* op_MMd[DMAX], ///< The array of \ref Multiarray_Matrix_d operators.
	 const int* op_values_i,                               ///< The index values of the current operator.
	 const struct Operator_Info* op_info,                  ///< \ref Sub_Operator_Info.
	 const bool* spec_indices                              ///< Defined for \ref set_ops_MMd.
	);

static void set_sub_operator_info
	(struct Sub_Operator_Info* sub_op_info, const struct Operator_Info* op_info, const struct Operators_TP* ops_tp,
	 const int ind_values)
{
	const struct Op_IO* op_io = op_info->op_io;
	assert(op_io[OP_IND_I].ce == 'v');

	const int* op_values = get_row_const_Matrix_i(ind_values,op_info->values_op);
	if (check_op_info_loss(op_values))
		EXIT_ERROR("Ensure that all is working as expected.\n");
		// Invalid condition in \ref set_ops_Md for j == OP_IND_H.

	// See comments in \ref element_operators_tp.h
	assert((op_values[OP_IND_D] == OP_INVALID_IND) ||
	       (op_values[OP_IND_CE] == OP_INVALID_IND && op_values[OP_IND_CE+1] == OP_INVALID_IND));

	const struct const_Multiarray_Matrix_d* op_MMd[] = { NULL, NULL, NULL };
	bool spec_indices[] = { false, false, false, };
	set_ops_MMd(op_MMd,spec_indices,ops_tp,op_values,op_info);

	const struct const_Matrix_d* op_Md[] = { NULL, NULL, NULL };
	set_ops_Md(op_Md,op_MMd,op_values,op_info,spec_indices);



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
static void set_spec_indices_diff
	(bool* spec_indices, ///< The array of indices for the special operator.
	 const int op_val_d, ///< The differentiation index of the current operator.
	 const int s_type    ///< \ref Element::s_type.
	);

/** \brief Set the indicator in the appropriate index for whether the second sub-element operator should be used for a
 *         face/edge computational element related operation. */
static void set_spec_indices_ce
	(bool* spec_indices,                 ///< The array of indices for the special operator.
	 const int op_val_ce,                ///< The computational element index of the current operator.
	 const struct Operator_Info* op_info ///< \ref Sub_Operator_Info.
	);

/// \brief Set the sub-operator indices for the input element type: ce_o, h_i, h_o.
static void set_sub_op_values
	(struct Matrix_i* sub_op_values,         ///< The sub-operator values for each of the tensor-product operators.
	 const struct const_Vector_i* op_values, ///< The operator values for the standard operator.
	 const struct Operator_Info* op_info,    ///< \ref Operator_Info.
	 const bool* spec_indices                ///< Defined for \ref set_ops_MMd.
	);

static void set_ops_MMd
	(const struct const_Multiarray_Matrix_d* op_MMd[DMAX], bool* spec_indices, const struct Operators_TP* ops_tp,
	 const int* op_values, const struct Operator_Info* op_info)
{
	const int op_val_d    = op_values[OP_IND_D];
	const int op_val_ce_o = op_values[OP_IND_CE+OP_IND_O];
	const int s_type      = op_info->element->s_type;

	assert(DMAX == 3);
	assert((op_val_d == OP_INVALID_IND) || (op_val_ce_o == OP_INVALID_IND));
	assert((s_type == ST_TP) || (s_type == ST_WEDGE));

	int* ind_ops_tp = ( s_type == ST_TP ? (int[]) { 0, 0, 0, } : (int[]) { 0, OP_INVALID_IND, 1, } );

	set_spec_indices_diff(spec_indices,op_val_d,s_type);
	set_spec_indices_ce(spec_indices,op_val_ce_o,op_info);

	for (int i = 0; i < DMAX; ++i) {
		const int ind_op_tp = ind_ops_tp[i];
		if (ind_op_tp != OP_INVALID_IND) {
			if (!spec_indices[i])
				op_MMd[i] = ops_tp->op[ind_op_tp][0];
			else
				op_MMd[i] = ops_tp->op[ind_op_tp][1];
		} else {
			op_MMd[i] = NULL;
		}
	}
}

static void set_ops_Md
	(const struct const_Matrix_d* op_Md[DMAX], const struct const_Multiarray_Matrix_d* op_MMd[DMAX],
	 const int* op_values_i, const struct Operator_Info* op_info, const bool* spec_indices)
{
	assert(op_values_i[OP_IND_H+OP_IND_I] == 0); // Not valid for fine to coarse (i.e. if info_loss == true)

/// \todo Move this to function above.
	const struct const_Vector_i* op_values = constructor_copy_const_Vector_i_i(OP_ORDER_MAX,op_values_i); // tbd

	struct Matrix_i* sub_op_values = constructor_empty_Matrix_i('R',DMAX,OP_ORDER_MAX); // tbd
	set_to_value_Matrix_i(sub_op_values,OP_INVALID_IND);
	const int set_inds[] = { OP_IND_P+OP_IND_I, OP_IND_P+OP_IND_O, OP_IND_H+OP_IND_I };
	for (int i = 0; i < (int)(sizeof(set_inds)/sizeof(*set_inds)); ++i)
		set_col_to_val_Matrix_i(i,sub_op_values,op_values->data[i]);

	set_sub_op_values(sub_op_values,op_values,op_info,spec_indices);

	EXIT_ADD_SUPPORT;
UNUSED(op_Md);
UNUSED(op_MMd);
}

// Level 1 ********************************************************************************************************** //

/// \brief Set the sub-operator indices for the QUAD element type.
static void set_sub_op_values_quad
	(struct Matrix_i* sub_op_values,         ///< Defined for \ref set_sub_op_values.
	 const struct const_Vector_i* op_values, ///< Defined for \ref set_sub_op_values.
	 const struct Operator_Info* op_info,    ///< Defined for \ref set_sub_op_values.
	 const bool* spec_indices                ///< Defined for \ref set_sub_op_values.
	);

/// \brief Set the sub-operator indices for the HEX element type.
static void set_sub_op_values_hex
	(struct Matrix_i* sub_op_values,         ///< Defined for \ref set_sub_op_values.
	 const struct const_Vector_i* op_values, ///< Defined for \ref set_sub_op_values.
	 const struct Operator_Info* op_info,    ///< Defined for \ref set_sub_op_values.
	 const bool* spec_indices                ///< Defined for \ref set_sub_op_values.
	);

/// \brief Set the sub-operator indices for the WEDGE element type.
static void set_sub_op_values_wedge
	(struct Matrix_i* sub_op_values,         ///< Defined for \ref set_sub_op_values.
	 const struct const_Vector_i* op_values, ///< Defined for \ref set_sub_op_values.
	 const struct Operator_Info* op_info,    ///< Defined for \ref set_sub_op_values.
	 const bool* spec_indices                ///< Defined for \ref set_sub_op_values.
	);

static void set_spec_indices_diff (bool* spec_indices, const int op_val_d, const int s_type)
{
	if (op_val_d == OP_INVALID_IND)
		return;

	assert(op_val_d >= 0);
	assert(op_val_d < DMAX);

	int ind_spec = -1;
	switch (s_type) {
	case ST_TP:
		ind_spec = op_val_d;
		break;
	case ST_WEDGE:
		if (op_val_d < 2)
			ind_spec = 0;
		else if (op_val_d == 2)
			ind_spec = 2;
		else
			EXIT_ERROR("Unsupported: %d\n",op_val_d);
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",s_type);
		break;
	}

	spec_indices[ind_spec] = true;
}

static void set_spec_indices_ce (bool* spec_indices, const int op_val_ce_o, const struct Operator_Info* op_info)
{
	if (op_val_ce_o == OP_INVALID_IND)
		return;

	const int e_type = op_info->element->type;
	const char ce_o  = op_info->op_io[OP_IND_O].ce;

/// \todo make external functions here?
	switch (e_type) {
	case QUAD:
		switch (ce_o) {
		case 'v':
			// Do nothing.
			break;
		case 'f': // fallthrough
		case 'e':
			switch (op_val_ce_o) {
			case H_QUAD1_F0:
			case H_QUAD1_F1:
			case H_QUAD4_F00: case H_QUAD4_F01:
			case H_QUAD4_F10: case H_QUAD4_F11:
				spec_indices[0] = true;
				break;
			case H_QUAD1_F2:
			case H_QUAD1_F3:
			case H_QUAD4_F20: case H_QUAD4_F21:
			case H_QUAD4_F30: case H_QUAD4_F31:
				spec_indices[1] = true;
				break;
			default:
				EXIT_ERROR("Unsupported: %d.\n",op_val_ce_o);
				break;
			}
			break;
		default:
			EXIT_ERROR("Unsupported: %c.\n",ce_o);
			break;
		}
		break;
	case HEX:
		switch (ce_o) {
		case 'v':
			// Do nothing.
			break;
		case 'f':
			switch (op_val_ce_o) {
			case H_HEX1_F0:
			case H_HEX1_F1:
			case H_HEX8_F00: case H_HEX8_F01: case H_HEX8_F02: case H_HEX8_F03:
			case H_HEX8_F10: case H_HEX8_F11: case H_HEX8_F12: case H_HEX8_F13:
				spec_indices[0] = true;
				break;
			case H_HEX1_F2:
			case H_HEX1_F3:
			case H_HEX8_F20: case H_HEX8_F21: case H_HEX8_F22: case H_HEX8_F23:
			case H_HEX8_F30: case H_HEX8_F31: case H_HEX8_F32: case H_HEX8_F33:
				spec_indices[1] = true;
				break;
			case H_HEX1_F4:
			case H_HEX1_F5:
			case H_HEX8_F40: case H_HEX8_F41: case H_HEX8_F42: case H_HEX8_F43:
			case H_HEX8_F50: case H_HEX8_F51: case H_HEX8_F52: case H_HEX8_F53:
				spec_indices[2] = true;
				break;
			default:
				EXIT_ERROR("Unsupported: %d.\n",op_val_ce_o);
				break;
			}
			break;
		default:
			EXIT_ERROR("Unsupported: %c.\n",ce_o);
			break;
		}
		break;
	case WEDGE:
		switch (ce_o) {
		case 'v':
			// Do nothing.
			break;
		case 'f':
			switch (op_val_ce_o) {
			case H_WEDGE1_F0:
			case H_WEDGE1_F1:
			case H_WEDGE1_F2:
			case H_WEDGE8_F00: case H_WEDGE8_F01: case H_WEDGE8_F02: case H_WEDGE8_F03:
			case H_WEDGE8_F10: case H_WEDGE8_F11: case H_WEDGE8_F12: case H_WEDGE8_F13:
			case H_WEDGE8_F20: case H_WEDGE8_F21: case H_WEDGE8_F22: case H_WEDGE8_F23:
				spec_indices[0] = true;
				break;
			case H_WEDGE1_F3:
			case H_WEDGE1_F4:
			case H_WEDGE8_F30: case H_WEDGE8_F31: case H_WEDGE8_F32: case H_WEDGE8_F33:
			case H_WEDGE8_F40: case H_WEDGE8_F41: case H_WEDGE8_F42: case H_WEDGE8_F43:
				spec_indices[2] = true;
				break;
			default:
				EXIT_ERROR("Unsupported: %d.\n",op_val_ce_o);
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
}

static void set_sub_op_values
	(struct Matrix_i* sub_op_values, const struct const_Vector_i* op_values, const struct Operator_Info* op_info,
	 const bool* spec_indices)
{
	switch (op_info->element->type) {
		case QUAD:  set_sub_op_values_quad(sub_op_values,op_values,op_info,spec_indices);  break;
		case HEX:   set_sub_op_values_hex(sub_op_values,op_values,op_info,spec_indices);   break;
		case WEDGE: set_sub_op_values_wedge(sub_op_values,op_values,op_info,spec_indices); break;
		default:    EXIT_ERROR("Unsupported: %d\n",op_info->element->type);   break;
	}
}

// Level 2 ********************************************************************************************************** //

static void set_sub_op_values_quad
	(struct Matrix_i* sub_op_values, const struct const_Vector_i* op_values, const struct Operator_Info* op_info,
	 const bool* spec_indices)
{
	const char ce_o = op_info->op_io[OP_IND_O].ce;

	int* s_vs[2] = { get_row_Matrix_i(0,sub_op_values),
	                 get_row_Matrix_i(1,sub_op_values), };

	// h_o
	const int d = op_info->element->d;
	const int h_o = OP_IND_H+OP_IND_O;
	const int op_val_h_o = op_values->data[h_o];
	switch (ce_o) {
	case 'v':
		switch (op_val_h_o) {
		case H_QUAD1_V0: s_vs[0][h_o] = H_LINE1_V0; s_vs[1][h_o] = H_LINE1_V0; break;

		case H_QUAD4_V0: s_vs[0][h_o] = H_LINE2_V0; s_vs[1][h_o] = H_LINE2_V0; break;
		case H_QUAD4_V1: s_vs[0][h_o] = H_LINE2_V1; s_vs[1][h_o] = H_LINE2_V0; break;
		case H_QUAD4_V2: s_vs[0][h_o] = H_LINE2_V0; s_vs[1][h_o] = H_LINE2_V1; break;
		case H_QUAD4_V3: s_vs[0][h_o] = H_LINE2_V1; s_vs[1][h_o] = H_LINE2_V1; break;

		default: EXIT_ERROR("Unsupported: %d.\n",op_val_h_o); break;
		}

		// differentiation index (if applicable)
		for (int i = 0; i < d; ++i) {
			if (spec_indices[i])
				s_vs[i][OP_IND_D] = 0;
		}
	case 'f': // fallthrough
	case 'e': {
		const int c_o = OP_IND_CE+OP_IND_O;
		const int op_val_ce_o = op_values->data[c_o];
		assert((0 <= op_val_ce_o) && (op_val_ce_o <= 3));

		int ind_f = op_val_ce_o / 2;
		s_vs[ind_f][h_o] = H_POINT1_V0;

		const int ind_v = ( ind_f == 0 ? 1 : 0);
		s_vs[ind_v][h_o] = op_val_h_o;

		// computational element index
		for (int i = 0; i < d; ++i) {
			if (spec_indices[i])
				s_vs[i][c_o] = op_val_ce_o % 2;
		}
	} default:
		EXIT_ERROR("Unsupported: %d\n",op_info->element->type);
		break;
	}
}
