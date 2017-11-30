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

#include "const_cast.h"
#include "element.h"
#include "element_operators.h"
#include "multiarray_operator.h"
#include "operator.h"
#include "simulation.h"

// Static function declarations ************************************************************************************* //

/// Container for information needed to set the sub-operators.
struct Sub_Operator_Info {
	const int d;      ///< The dimension of the input.
	const int s_type; ///< The \ref Element::s_type for which the operator is being constructed.

	const struct const_Matrix_d* sub_op_Md[DMAX]; ///< The array of matrix sub-operators.
};

/// \brief Constructor standard operators from the tensor-product sub-operators.
static void construct_operators_std
	(const struct Multiarray_Operator* op ///< The \ref Multiarray_Operator\* input.
	);

/// \brief Set the pointers in \ref Sub_Operator_Info::sub_op_Md to the appropriate sub-operator matrices.
static void set_sub_operator_info
	(struct Sub_Operator_Info* sub_op_info, ///< \ref Sub_Operator_Info.
	 const struct Operator_Info* op_info,   ///< \ref Sub_Operator_Info.
	 const struct Operators_TP* ops_tp,     ///< \ref Operators_TP.
	 const int ind_values                   ///< The index of the row of \ref Operator_Info::values_op to use.
	);

/** \brief Set the row permutation indices for the standard operator based on the input and desired output index
 *         ordering. */
static void set_row_permutation_indices
	(ptrdiff_t* perm,            ///< Set to the permutation indices.
	 const int ordering_i[DMAX], ///< The input ordering.
	 const int ordering_o[DMAX], ///< The output ordering.
	 const int*const n_rows      ///< The number of rows in each of the DMAX operator directions.
	);

// Interface functions ********************************************************************************************** //

void set_operators_tp
	(struct Operators_TP* ops_tp,
	 const struct Multiarray_Operator* op_00, const struct Multiarray_Operator* op_01,
	 const struct Multiarray_Operator* op_10, const struct Multiarray_Operator* op_11)
{
	assert(op_00 != NULL);
	if (op_01 == NULL)
		op_01 = op_00;
	if (op_10 == NULL)
		op_10 = op_00;
	if (op_11 == NULL)
		op_11 = op_10;

	ops_tp->op[0][0] = op_00;
	ops_tp->op[0][1] = op_01;
	ops_tp->op[1][0] = op_10;
	ops_tp->op[1][1] = op_11;
}

const struct Multiarray_Operator* constructor_operators_tp
	(const char*const name_type, const char*const name_in, const char*const name_out, const char*const name_range,
	 const struct const_Element* element, const struct Simulation* sim, const struct Operators_TP* ops_tp)
{
	struct Operator_Info* op_info =
		constructor_Operator_Info(name_type,name_in,name_out,name_range,element,sim); // destructed

	const struct Op_IO* op_io = op_info->op_io;

	const char ce_i = op_io[OP_IND_I].ce;

	assert(ce_i == 'v'); // May potentially be made flexible in future.
	const int d = element->d;

	const struct Multiarray_Operator* op = constructor_empty_Multiarray_Operator_V(op_info->extents_op); // returned

	const ptrdiff_t row_max = op_info->values_op->ext_0;
	for (ptrdiff_t row = 0; row < row_max; ++row) {
		struct Sub_Operator_Info sub_op_info =
			{ .d         = element->d,
			  .s_type    = element->s_type,
			  .sub_op_Md = { NULL, NULL, NULL, }, };

		set_sub_operator_info(&sub_op_info,op_info,ops_tp,(int)row);

		const int order_op = (int)op_info->extents_op->ext_0;
		const int* op_values = get_row_const_Matrix_i(row,op_info->values_op);
		const struct const_Vector_i* indices_op =
			constructor_indices_Vector_i(order_op,op_values,NULL); // destructed
		const ptrdiff_t ind_op = compute_index_sub_container_pi(op->order,0,op->extents,indices_op->data);
		destructor_const_Vector_i(indices_op);

		// Set the tensor-product operators
		const_constructor_move_Multiarray_Matrix_d(
			&op->data[ind_op]->ops_tp,constructor_empty_Multiarray_Matrix_d(false,1,(ptrdiff_t[]){d})); // keep

		// Sub-operators are always owned by sub-elements.
		const_cast_b(&op->data[ind_op]->ops_tp->owns_data,false);
		for (int i = 0; i < d; ++i)
			const_constructor_move_const_Matrix_d(&op->data[ind_op]->ops_tp->data[i],sub_op_info.sub_op_Md[i]);
	}
	construct_operators_std(op);

	destructor_Operator_Info(op_info);

	return op;
}

const struct const_Matrix_d* constructor_op_std (const struct const_Multiarray_Matrix_d* ops_tp)
{
	/** If this function is found to be slow while profiling, the following changes may be implemented:
	 *  - replace the dense blas function calls with block sparse blas mm_d function calls when applying the op_s
	 *    and op_t matrices.
	 *  - remove the redundant permutation between rs and rst for hex operators (i.e. {1,0,2} -> {2,0,1} directly).
	 */
	assert(ops_tp->order == 1);
	const int d_op = (int)compute_size(ops_tp->order,ops_tp->extents);
	assert(2 <= d_op);
	assert(d_op <= DMAX);

	int n_rows_sub[DMAX] = { 0, 0, 0, },
	    n_cols_sub[DMAX] = { 0, 0, 0, };
	set_ops_tp_n_rows_cols(n_rows_sub,n_cols_sub,ops_tp);

	int n_blocks = -1;
	const struct const_Matrix_d* sub_op = NULL;
	struct Vector_i* n_rows_op =
		constructor_copy_Vector_i_i(DMAX,(int[]){n_cols_sub[0],n_cols_sub[1],n_cols_sub[2]}); // destructed

	// r-direction
	int ind_sub_op = 0;
	sub_op = ops_tp->data[ind_sub_op];
	assert(sub_op != NULL);

	n_rows_op->data[ind_sub_op] = n_rows_sub[ind_sub_op];
	n_blocks = n_rows_op->data[1]*n_rows_op->data[2];
	const struct const_Matrix_d* op_r = constructor_block_diagonal_const_Matrix_d(sub_op,n_blocks); // destructed

	// s-direction
	++ind_sub_op;
	sub_op = ops_tp->data[ind_sub_op];

	const struct const_Matrix_d* op_rs = NULL;
	if (sub_op) {
		ptrdiff_t perm_r[prod_Vector_i(n_rows_op)];
		set_row_permutation_indices(perm_r,(int[]){0,1,2},(int[]){1,0,2},n_rows_op->data);

		permute_Matrix_d((struct Matrix_d*)op_r,perm_r);

		n_rows_op->data[ind_sub_op] = n_rows_sub[ind_sub_op];
		n_blocks = n_rows_op->data[0]*n_rows_op->data[2];
		// Note: Potentially very sparse.
		const struct const_Matrix_d* op_s =
			constructor_block_diagonal_const_Matrix_d(sub_op,n_blocks); // destructed

		op_rs = constructor_mm_NN1R_const_Matrix_d(op_s,op_r); // destructed
		destructor_const_Matrix_d(op_s);
		destructor_const_Matrix_d(op_r);

		ptrdiff_t perm_rs[prod_Vector_i(n_rows_op)];
		set_row_permutation_indices(perm_rs,(int[]){1,0,2},(int[]){0,1,2},n_rows_op->data);

		permute_Matrix_d((struct Matrix_d*)op_rs,perm_rs);
	} else {
		assert(d_op == 3);
		n_rows_op->data[ind_sub_op] = n_rows_sub[ind_sub_op];
		op_rs = op_r;
	}

	// t-direction
	const struct const_Matrix_d* op_rst = NULL;
	if (d_op > 2) {
		++ind_sub_op;
		sub_op = ops_tp->data[ind_sub_op];
		assert(sub_op != NULL);

		ptrdiff_t perm_rs[prod_Vector_i(n_rows_op)];
		set_row_permutation_indices(perm_rs,(int[]){0,1,2},(int[]){2,0,1},n_rows_op->data);

		permute_Matrix_d((struct Matrix_d*)op_rs,perm_rs);

		n_rows_op->data[ind_sub_op] = n_rows_sub[ind_sub_op];
		n_blocks = n_rows_op->data[0]*n_rows_op->data[1];
		// Note: Potentially very sparse.
		const struct const_Matrix_d* op_t =
			constructor_block_diagonal_const_Matrix_d(sub_op,n_blocks); // destructed

		op_rst = constructor_mm_NN1R_const_Matrix_d(op_t,op_rs); // destructed
		destructor_const_Matrix_d(op_t);
		destructor_const_Matrix_d(op_rs);

		ptrdiff_t perm_rst[prod_Vector_i(n_rows_op)];
		set_row_permutation_indices(perm_rst,(int[]){2,0,1},(int[]){0,1,2},n_rows_op->data);

		permute_Matrix_d((struct Matrix_d*)op_rst,perm_rst);
	} else {
		assert(d_op == 2);
		op_rst = op_rs;
	}

	destructor_Vector_i(n_rows_op);
	return op_rst;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Set the input `op_MO[DMAX]` to point to the appropriate operators from `ops_tp`.
 *
 *  `spec_indices` can be interpreted as the direction(s) (n-Cube: r,s,t; wedge: rs,NULL,t) in which to use the
 *  'special' operator (the second operator for the sub-element).
 */
static void set_ops_MO
	(const struct Multiarray_Operator* op_MO[DMAX], ///< The array of \ref Multiarray_Matrix_d operators.
	 bool* spec_indices,                            ///< The array of indices for the special operators.
	 const struct Operators_TP* ops_tp,             ///< \ref Operators_TP.
	 const struct const_Vector_i* op_values,        ///< The index values of the current operator.
	 const struct Operator_Info* op_info            ///< \ref Sub_Operator_Info.
	);

/// \brief Set the the input `op_Md[DMAX]` to point to the appropriate matrices in `op_MO[DMAX]`.
static void set_ops_Md
	(const struct const_Matrix_d* op_Md[DMAX],      ///< The array of \ref Matrix_d operators.
	 const struct Multiarray_Operator* op_MO[DMAX], ///< The array of \ref Multiarray_Matrix_d operators.
	 const struct const_Vector_i* op_values,        ///< The index values of the current operator.
	 const struct Operator_Info* op_info,           ///< \ref Sub_Operator_Info.
	 const bool* spec_indices                       ///< Defined for \ref set_ops_MO.
	);

static void set_sub_operator_info
	(struct Sub_Operator_Info* sub_op_info, const struct Operator_Info* op_info, const struct Operators_TP* ops_tp,
	 const int ind_values)
{
	const struct Op_IO* op_io = op_info->op_io;
	assert(op_io[OP_IND_I].ce == 'v');

	const struct const_Matrix_i* values_op = op_info->values_op;
	const struct const_Vector_i* op_values =
		constructor_move_const_Vector_Matrix_row_i(ind_values,values_op,false); // destructed

	if (op_should_use_L2(op_values->data,op_io))
		EXIT_ERROR("Ensure that all is working as expected.\n");
		// Invalid condition in \ref set_ops_Md for j == OP_IND_H.

	// See comments in \ref element_operators_tp.h
	assert((op_values->data[OP_IND_D] == OP_INVALID_IND) ||
	       (op_values->data[OP_IND_CE] == OP_INVALID_IND && op_values->data[OP_IND_CE+1] == OP_INVALID_IND));

	const struct Multiarray_Operator* op_MO[] = { NULL, NULL, NULL };
	bool spec_indices[] = { false, false, false, };
	set_ops_MO(op_MO,spec_indices,ops_tp,op_values,op_info);
	set_ops_Md(sub_op_info->sub_op_Md,op_MO,op_values,op_info,spec_indices);

	destructor_const_Vector_i(op_values);
}

static void construct_operators_std (const struct Multiarray_Operator* op)
{
	const ptrdiff_t size = compute_size(op->order,op->extents);
	for (int i = 0; i < size; ++i) {
		const struct Operator* op_c = op->data[i];
		if (op_c->ops_tp)
			const_constructor_move_const_Matrix_d(&op_c->op_std,constructor_op_std(op_c->ops_tp)); // keep
	}
}

static void set_row_permutation_indices
	(ptrdiff_t* perm, const int ordering_i[DMAX], const int ordering_o[DMAX], const int*const n_rows)
{
	const int i_0 = ordering_i[0],
	          i_1 = ordering_i[1],
	          i_2 = ordering_i[2];
	const int o_0 = ordering_o[0],
	          o_1 = ordering_o[1],
	          o_2 = ordering_o[2];

	int ind[DMAX];
	for (ind[o_2] = 0; ind[o_2] < n_rows[o_2]; ++ind[o_2]) {
	for (ind[o_1] = 0; ind[o_1] < n_rows[o_1]; ++ind[o_1]) {
	for (ind[o_0] = 0; ind[o_0] < n_rows[o_0]; ++ind[o_0]) {
		*perm++ = ind[i_0]+n_rows[i_0]*(ind[i_1]+n_rows[i_1]*(ind[i_2]));
	}}}
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
	 const bool* spec_indices                ///< Defined for \ref set_ops_MO.
	);

static void set_ops_MO
	(const struct Multiarray_Operator* op_MO[DMAX], bool* spec_indices, const struct Operators_TP* ops_tp,
	 const struct const_Vector_i* op_values, const struct Operator_Info* op_info)
{
	const int op_val_d    = op_values->data[OP_IND_D];
	const int op_val_ce_o = op_values->data[OP_IND_CE+OP_IND_O];
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
				op_MO[i] = ops_tp->op[ind_op_tp][0];
			else
				op_MO[i] = ops_tp->op[ind_op_tp][1];
		} else {
			op_MO[i] = NULL;
		}
	}
}

static void set_ops_Md
	(const struct const_Matrix_d* op_Md[DMAX], const struct Multiarray_Operator* op_MO[DMAX],
	 const struct const_Vector_i* op_values, const struct Operator_Info* op_info, const bool* spec_indices)
{
	assert(op_values->data[OP_IND_H+OP_IND_I] == 0); // Not valid for fine to coarse (i.e. if info_loss == true)

	const int d = op_info->element->d;
	struct Matrix_i* sub_op_values = constructor_empty_Matrix_i('R',d,OP_ORDER_MAX); // destructed
	set_to_value_Matrix_i(sub_op_values,OP_INVALID_IND);
	const int set_inds[] = { OP_IND_P+OP_IND_I, OP_IND_P+OP_IND_O, OP_IND_H+OP_IND_I };
	for (int i = 0; i < (int)(sizeof(set_inds)/sizeof(*set_inds)); ++i) {
		const int col = set_inds[i];
		set_col_to_val_Matrix_i(col,sub_op_values,op_values->data[col]);
	}

	set_sub_op_values(sub_op_values,op_values,op_info,spec_indices);

	for (int i = 0; i < d; ++i) {
		const int* op_values = get_row_Matrix_i(i,sub_op_values);
		const struct const_Vector_i* inds_op = constructor_indices_Vector_i(-1,op_values,NULL); // destructed
		const int ind_sub_op = (int)compute_index_sub_container_pi(op_MO[i]->order,0,op_MO[i]->extents,inds_op->data);
		destructor_const_Vector_i(inds_op);

		op_Md[i] = op_MO[i]->data[ind_sub_op]->op_std;
	}
	destructor_Matrix_i(sub_op_values);
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

/** \brief Set the missing dimension indices based on the input index and the dimension.
 *  \return A pointer to an array of size d-1 holding the missing dimension indices. */
static const int* set_missing_dimension_indices
	(const int dim, ///< The dimension index which is already accounted for.
	 const int d    ///< The dimension.
	);

static void set_sub_op_values_quad
	(struct Matrix_i* sub_op_values, const struct const_Vector_i* op_values, const struct Operator_Info* op_info,
	 const bool* spec_indices)
{
	const char ce_o = op_info->op_io[OP_IND_O].ce;

	int* s_vs[2] = { get_row_Matrix_i(0,sub_op_values),
	                 get_row_Matrix_i(1,sub_op_values), };

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
		break;
	case 'f': // fallthrough
	case 'e': {
		const int c_o = OP_IND_CE+OP_IND_O;
		const int op_val_ce_o = op_values->data[c_o];

		assert((0 <= op_val_ce_o) && (op_val_ce_o < op_info->element->n_f));

		const int ind_f = op_val_ce_o / 2;
		s_vs[ind_f][h_o] = H_POINT1_V0;
		const int* ind_v = set_missing_dimension_indices(ind_f,d);
		s_vs[ind_v[0]][h_o] = op_val_h_o;

		// computational element index
		for (int i = 0; i < d; ++i) {
			if (spec_indices[i])
				s_vs[i][c_o] = op_val_ce_o % 2;
		}
		break;
	} default:
		EXIT_ERROR("Unsupported: %d\n",op_info->element->type);
		break;
	}
}

static void set_sub_op_values_hex
	(struct Matrix_i* sub_op_values, const struct const_Vector_i* op_values, const struct Operator_Info* op_info,
	 const bool* spec_indices)
{
	const char ce_o = op_info->op_io[OP_IND_O].ce;

	int* s_vs[3] = { get_row_Matrix_i(0,sub_op_values),
	                 get_row_Matrix_i(1,sub_op_values),
	                 get_row_Matrix_i(2,sub_op_values), };

	const int d = op_info->element->d;
	const int h_o = OP_IND_H+OP_IND_O;
	const int op_val_h_o = op_values->data[h_o];
	switch (ce_o) {
	case 'v':
		switch (op_val_h_o) {
		case H_HEX1_V0: s_vs[0][h_o] = H_LINE1_V0; s_vs[1][h_o] = H_LINE1_V0; s_vs[2][h_o] = H_LINE1_V0; break;

		case H_HEX8_V0: s_vs[0][h_o] = H_LINE2_V0; s_vs[1][h_o] = H_LINE2_V0; s_vs[2][h_o] = H_LINE2_V0; break;
		case H_HEX8_V1: s_vs[0][h_o] = H_LINE2_V1; s_vs[1][h_o] = H_LINE2_V0; s_vs[2][h_o] = H_LINE2_V0; break;
		case H_HEX8_V2: s_vs[0][h_o] = H_LINE2_V0; s_vs[1][h_o] = H_LINE2_V1; s_vs[2][h_o] = H_LINE2_V0; break;
		case H_HEX8_V3: s_vs[0][h_o] = H_LINE2_V1; s_vs[1][h_o] = H_LINE2_V1; s_vs[2][h_o] = H_LINE2_V0; break;
		case H_HEX8_V4: s_vs[0][h_o] = H_LINE2_V0; s_vs[1][h_o] = H_LINE2_V0; s_vs[2][h_o] = H_LINE2_V1; break;
		case H_HEX8_V5: s_vs[0][h_o] = H_LINE2_V1; s_vs[1][h_o] = H_LINE2_V0; s_vs[2][h_o] = H_LINE2_V1; break;
		case H_HEX8_V6: s_vs[0][h_o] = H_LINE2_V0; s_vs[1][h_o] = H_LINE2_V1; s_vs[2][h_o] = H_LINE2_V1; break;
		case H_HEX8_V7: s_vs[0][h_o] = H_LINE2_V1; s_vs[1][h_o] = H_LINE2_V1; s_vs[2][h_o] = H_LINE2_V1; break;

		default: EXIT_ERROR("Unsupported: %d.\n",op_val_h_o); break;
		}

		// differentiation index (if applicable)
		for (int i = 0; i < d; ++i) {
			if (spec_indices[i])
				s_vs[i][OP_IND_D] = 0;
		}
		break;
	case 'f': {
		const int c_o = OP_IND_CE+OP_IND_O;
		const int op_val_ce_o = op_values->data[c_o];

		assert((0 <= op_val_ce_o) && (op_val_ce_o < op_info->element->n_f));

		const int ind_f = op_val_ce_o / 2;
		s_vs[ind_f][h_o] = H_POINT1_V0;

		const int* ind_v = set_missing_dimension_indices(ind_f,d);

/// \todo call set_sub_op_values_quad 'v' here H_QUAD4_V0:H_QUAD4_V3 here.
		switch (op_val_h_o) {
			case H_QUAD1_V0: s_vs[ind_v[0]][h_o] = H_LINE1_V0; s_vs[ind_v[1]][h_o] = H_LINE1_V0; break;

			case H_QUAD4_V0: s_vs[ind_v[0]][h_o] = H_LINE2_V0; s_vs[ind_v[1]][h_o] = H_LINE2_V0; break;
			case H_QUAD4_V1: s_vs[ind_v[0]][h_o] = H_LINE2_V1; s_vs[ind_v[1]][h_o] = H_LINE2_V0; break;
			case H_QUAD4_V2: s_vs[ind_v[0]][h_o] = H_LINE2_V0; s_vs[ind_v[1]][h_o] = H_LINE2_V1; break;
			case H_QUAD4_V3: s_vs[ind_v[0]][h_o] = H_LINE2_V1; s_vs[ind_v[1]][h_o] = H_LINE2_V1; break;

			default: EXIT_ERROR("Unsupported: %d.\n",op_val_h_o); break;
		}

		for (int i = 0; i < d; ++i) {
			if (spec_indices[i])
				s_vs[i][c_o] = op_val_ce_o % 2;
		}
		break;
	} case 'e': {
		assert(op_val_h_o == 0); // Need to add support for h-refinement edge operators if necessary.

		const int c_o = OP_IND_CE+OP_IND_O;
		const int op_val_ce_o = op_values->data[c_o];
		assert((0 <= op_val_ce_o) && (op_val_ce_o < op_info->element->n_e));

		const int ind_v = op_val_ce_o / 4;
		s_vs[ind_v][h_o] = H_LINE1_V0;

		const int* ind_f = set_missing_dimension_indices(ind_v,d);

		s_vs[ind_f[0]][h_o] = H_POINT1_V0;
		s_vs[ind_f[1]][h_o] = H_POINT1_V0;

		// computational element index
		assert(spec_indices[ind_f[0]] == true);
		assert(spec_indices[ind_f[1]] == true);

		s_vs[ind_f[0]][c_o] = op_val_ce_o % 2;
		s_vs[ind_f[1]][c_o] = (op_val_ce_o % 4) / 2;
		break;
	} default:
		EXIT_ERROR("Unsupported: %d\n",op_info->element->type);
		break;
	}
}

static void set_sub_op_values_wedge
	(struct Matrix_i* sub_op_values, const struct const_Vector_i* op_values, const struct Operator_Info* op_info,
	 const bool* spec_indices)
{
	EXIT_ADD_SUPPORT;
UNUSED(sub_op_values);
UNUSED(op_values);
UNUSED(op_info);
UNUSED(spec_indices);
}

// Level 3 ********************************************************************************************************** //

static const int* set_missing_dimension_indices (const int dim, const int d)
{
	static int missed_inds[2] = { -1, -1, };
	switch (d) {
	case 2:
		switch (dim) {
			case 0: missed_inds[0] = 1; break;
			case 1: missed_inds[0] = 0; break;
			default: EXIT_ERROR("Unsupported: %d\n",dim); break;
		}
		break;
	case 3:
		switch (dim) {
			case 0: missed_inds[0] = 1; missed_inds[1] = 2; break;
			case 1: missed_inds[0] = 0; missed_inds[1] = 2; break;
			case 2: missed_inds[0] = 0; missed_inds[1] = 1; break;
			default: EXIT_ERROR("Unsupported: %d\n",dim); break;
		}
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",d);
		break;
	}

	return missed_inds;
}
