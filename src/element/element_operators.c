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

#include "element_operators.h"

#include <assert.h>
#include <string.h>
#include <ctype.h>
#include "gsl/gsl_math.h"

#include "macros.h"
#include "definitions_bases.h"
#include "definitions_core.h"
#include "definitions_element_operators.h"
#include "definitions_elements.h"
#include "definitions_h_ref.h"
#include "definitions_test_case.h"

#include "multiarray.h"
#include "matrix.h"
#include "vector.h"

#include "element.h"

#include "bases.h"
#include "const_cast.h"
#include "multiarray_operator.h"
#include "nodes.h"
#include "nodes_correspondence.h"
#include "nodes_operators.h"
#include "operator.h"
#include "simulation.h"

// Static function declarations ************************************************************************************* //

/** \brief Version of \ref set_operator_fptr setting standard operator(s).
 *  \note Standard operators include all permutations of coef/value to coef/value with possible differentiation.
 *
 *  Currently, the nodes are recomputed upon each entry into this function. It would be more efficient to compute the
 *  required nodes rules in advance and pass a multiarray of nodes to this function. This would also allow for more
 *  efficient computation of the nodes themselves.
 */
static void set_operator_std
	(ptrdiff_t*const ind_values,           ///< See brief.
	 const struct Multiarray_Operator* op, ///< See brief.
	 const struct Operator_Info* op_info,  ///< See brief.
	 const struct Simulation* sim          ///< See brief.
	);

/** \brief Version of \ref set_operator_fptr setting solver operator(s).
 *  \ref set_operator_std may contain relevant comments. */
static void set_operator_solver
	(ptrdiff_t*const ind_values,           ///< Defined for \ref set_operator_fptr.
	 const struct Multiarray_Operator* op, ///< Defined for \ref set_operator_fptr.
	 const struct Operator_Info* op_info,  ///< Defined for \ref set_operator_fptr.
	 const struct Simulation* sim          ///< Defined for \ref set_operator_fptr.
	);

/// \brief Set the node correspondence for the supported arrangements.
static void set_node_correspondence
	(ptrdiff_t ind_values,                       ///< Defined for \ref set_operator_fptr.
	 const struct const_Multiarray_Vector_i* nc, ///< Container for the node correspondence data.
	 const struct Operator_Info* op_info,        ///< Defined for \ref set_operator_fptr.
	 const struct Simulation* sim                ///< Defined for \ref set_operator_fptr.
	);

/// \brief Set the cubature weights for the supported arrangements.
static void set_weights
	(ptrdiff_t ind_values,                      ///< Defined for \ref set_operator_fptr.
	 const struct const_Multiarray_Vector_d* w, ///< Container for the cubature weight data.
	 const struct Operator_Info* op_info,       ///< Defined for \ref set_operator_fptr.
	 const struct Simulation* sim               ///< Defined for \ref set_operator_fptr.
	);

/** \brief Convert the `char*` input to the appropriate definition of OP_T_*.
 *  \return See brief. */
static int convert_to_type
	(const char* name_type ///< Defined for \ref constructor_Operator_Info.
	);

/** \brief Check if the computed operators should be transposed.
 *  \return `true` if yes; `false` otherwise. */
static int check_for_transpose
	(const char* name_type, ///< Defined for \ref constructor_Operator_Info.
	 const int op_type      ///< \ref Operator_Info::op_type.
	);

/** \brief Convert the `char*` input to the appropriate definition of OP_R_D_*.
 *  \return See brief. */
static int convert_to_range_d
	(const char* name_type ///< Defined for \ref constructor_Operator_Info.
	);

/** \brief Convert the `char*` inputs to the appropriate definition of OP_R_CE_*.
 *  \return See brief. */
static int convert_to_range_ce
	(const char ce_i, ///< `name_in[0]`  which is defined for \ref constructor_Operator_Info.
	 const char ce_o  ///< `name_out[0]` which is defined for \ref constructor_Operator_Info.
	);

/** \brief Convert the `char*` input to the appropriate definition of OP_R_*.
 *  \return See brief. */
static int convert_to_range
	(const char type_range,      ///< The type of range parameter. Options: 'h', 'p'.
	 const char*const name_range ///< Defined for \ref constructor_Operator_Info.
	);

/// \brief Set up \ref Operator_Info extents_* members.
static void set_up_extents
	(struct Operator_Info* op_info ///< \ref Operator_Info.
	);

/// \brief Set up \ref Operator_Info::values_op.
static void set_up_values_op
	(struct Operator_Info* op_info ///< \ref Operator_Info.
	);

/// \brief Transpose the computed operators.
static void transpose_operators
	(struct Multiarray_Operator* op ///< Multiarray of operators.
	);

/** \brief Return the names of the basis types to use for the constructor in \ref constructor_operators_bt.
 *  \return See brief. */
static const int* get_basis_types
	(const char* name_type,       ///< Defined for \ref constructor_operators_bt.
	 const struct Simulation* sim ///< \ref Simulation.
	);

/** \brief Return the index of the operator currently being set.
 *  \return See brief. */
static ptrdiff_t get_ind_op
	(const struct Operator_Info* op_info, ///< \ref Operator_Info.
	 const int* op_values,                ///< One line of values of \ref Operator_Info::values_op.
	 const struct Multiarray_Operator* op ///< Defined for \ref set_operator_fptr.
	);

/// \brief Set the values of \ref Op_IO for the current operator.
static void set_current_op_io
	(const struct Operator_Info* op_info, ///< \ref Operator_Info.
	 const int* op_values                 ///< One line of values of \ref Operator_Info::values_op.
	);

/** \brief Constructor for the coefficient to value operator, selecting the operator basis based on `kind`.
 *  \return See brief. */
static const struct const_Matrix_d* constructor_cv
	(const struct const_Matrix_d* rst, ///< The rst coordinates.
	 const struct Op_IO* op_io,        ///< \ref Op_IO.
	 const struct Simulation* sim,     ///< \ref Simulation.
	 const int basis_type_i            ///< May be set to hold the input basis type.
	);

// Interface functions ********************************************************************************************** //

const struct Multiarray_Operator* constructor_operators
	(const char*const name_type, const char*const name_in, const char*const name_out, const char*const name_range,
	 const struct const_Element* element, const struct Simulation* sim)
{
	struct Operator_Info* op_info =
		constructor_Operator_Info(name_type,name_in,name_out,name_range,element,sim); // destructed
	assert((op_info->range_d == OP_R_D_0) || (op_info->range_d == OP_R_D_ALL));

	const struct Multiarray_Operator* op = constructor_empty_Multiarray_Operator_V(op_info->extents_op); // returned

	const ptrdiff_t row_max = op_info->values_op->ext_0;
	for (ptrdiff_t row = 0; row < row_max; ) // row is incremented when setting the operators.
		op_info->set_operator(&row,op,op_info,sim);

	if (op_info->transpose)
		transpose_operators((struct Multiarray_Operator*)op);
	destructor_Operator_Info(op_info);

	return op;
}

const struct const_Multiarray_Vector_i* constructor_operators_nc
	(const char*const name_in, const char*const name_out, const char*const name_range,
	 const struct const_Element* element, const struct Simulation* sim)
{
	// This will fail for PYR/WEDGE. Ensure that all is working as expected (that operators are built for both face
	// element types).
	assert(get_number_of_face_elements(element) == 1);

	struct Operator_Info* op_info =
		constructor_Operator_Info("UNUSED0",name_in,name_out,name_range,element,sim); // destructed

	const struct const_Vector_i* e_o = op_info->extents_op;

	const struct const_Element* f_element = element->face_element[0];
	assert(f_element != NULL);

	// Add an additional index for the possible permutations
	const ptrdiff_t len_e = e_o->ext_0+1;
	struct Vector_i* e_o_p1 = constructor_empty_Vector_i(len_e); // destructed
	e_o_p1->data[0] = (int)get_n_perm_corr(f_element->d,f_element->s_type);
	for (int i = 0; i < len_e-1; i++)
		e_o_p1->data[i+1] = e_o->data[i];

	const struct const_Multiarray_Vector_i* nc =
		constructor_empty_const_Multiarray_Vector_i_V(false,(struct const_Vector_i*)e_o_p1); // returned
	destructor_Vector_i(e_o_p1);

	const ptrdiff_t row_max = op_info->values_op->ext_0;
	for (ptrdiff_t row = 0; row < row_max; ++row)
		set_node_correspondence(row,nc,op_info,sim);
	destructor_Operator_Info(op_info);

	return nc;
}

const struct const_Multiarray_Vector_d* constructor_operators_w
	(const char*const name_in, const char*const name_out, const char*const name_range, const int p_ref[2],
	 const struct const_Element* element, const struct Simulation* sim)
{
UNUSED(p_ref);
	struct Operator_Info* op_info =
		constructor_Operator_Info("UNUSED0",name_in,name_out,name_range,element,sim); // destructed

	const struct const_Multiarray_Vector_d* w =
		constructor_empty_const_Multiarray_Vector_d_V(false,op_info->extents_op); // returned

	const ptrdiff_t row_max = op_info->values_op->ext_0;
	for (ptrdiff_t row = 0; row < row_max; ++row)
		set_weights(row,w,op_info,sim);
	destructor_Operator_Info(op_info);

	return w;
}

const struct Multiarray_Operator* constructor_operators_bt
	(const char*const name_type, const char*const name_in, const char*const name_out, const char*const name_range,
	 const struct const_Element* element, const struct Simulation* sim)
{
	struct Operator_Info* op_info =
		constructor_Operator_Info(name_type,name_in,name_out,name_range,element,sim); // destructed
	assert(op_info->range_d == OP_R_D_0);

	const struct Multiarray_Operator* op = constructor_empty_Multiarray_Operator_V(op_info->extents_op); // returned

	const int* basis_type = get_basis_types(name_type,sim);

	const ptrdiff_t row_max = op_info->values_op->ext_0;
	for (ptrdiff_t row = 0; row < row_max; ++row) {
		const int* op_values = get_row_const_Matrix_i(row,op_info->values_op);
		set_current_op_io(op_info,op_values);

		const struct Op_IO* op_io = op_info->op_io;
		const struct const_Nodes* nodes = constructor_const_Nodes_h(OP_IND_I,op_io,op_info->element,sim); // dest.

		const struct const_Matrix_d* cv_i = constructor_cv(nodes->rst,op_io,sim,basis_type[0]), // destructed
		                           * cv_o = constructor_cv(nodes->rst,op_io,sim,basis_type[1]); // destructed
		destructor_const_Nodes(nodes);


		const struct const_Matrix_d* op_io0 = constructor_sgesv_const_Matrix_d(cv_o,cv_i); // moved
		destructor_const_Matrix_d(cv_i);
		destructor_const_Matrix_d(cv_o);

		ptrdiff_t ind_op = get_ind_op(op_info,op_values,op);
		const_constructor_move_const_Matrix_d(&op->data[ind_op]->op_std,op_io0);
	}
	destructor_Operator_Info(op_info);

	return op;
}

const struct const_Multiarray_Vector_d* constructor_operators_ones_coef
	(const char*const name_io, const struct const_Element* element, const struct Simulation* sim)
{
	struct Operator_Info* op_info =
		constructor_Operator_Info("vv0",name_io,name_io,"H_1_P_PM0",element,sim); // destructed

	const struct Multiarray_Operator* op = constructor_empty_Multiarray_Operator_V(op_info->extents_op); // destructed

	const ptrdiff_t row_max = op_info->values_op->ext_0;
	for (ptrdiff_t row = 0; row < row_max; ) // row is incremented when setting the operators.
		op_info->set_operator(&row,op,op_info,sim);

	const struct const_Multiarray_Vector_d* ones_coef =
		constructor_empty_const_Multiarray_Vector_d_V(false,op_info->extents_op); // returned
	destructor_Operator_Info(op_info);

	const ptrdiff_t size = compute_size(op->order,op->extents);
	for (ptrdiff_t i = 0; i < size; ++i) {
		const struct const_Matrix_d*const vv0 = op->data[i]->op_std;
		if (vv0 == NULL) {
			const_constructor_move_const_Vector_d(&ones_coef->data[i],NULL);
			continue;
		}

		struct Vector_d*const ones = constructor_empty_Vector_d(vv0->ext_0); // destructed
		set_to_value_Vector_d(ones,1.0);

		const_constructor_move_const_Vector_d(&ones_coef->data[i],
			constructor_sgesv_const_Vector_d(vv0,(struct const_Vector_d*)ones)); // keep
		destructor_Vector_d(ones);
	}
	destructor_Multiarray_Operator(op);

	return ones_coef;
}

const struct Multiarray_Operator* constructor_operators_tens3
	(const struct Multiarray_Operator*const op_l, const struct Multiarray_Operator*const op_r)
{
	/* If this function is ever to be used for operators other than those having `ce_i = ce_o = 'v'`, changes will
	 * likely be necessary. */
	const int order               = op_r->order;
	const ptrdiff_t*const extents = op_r->extents;
	const ptrdiff_t size          = compute_size(order,extents);

	const struct Multiarray_Operator* op = constructor_empty_Multiarray_Operator(order,extents); // returned

	ptrdiff_t counter[order];
	for (ptrdiff_t i = 0; i < order; ++i)
		counter[i] = 0;

	for (ptrdiff_t sz = 0; sz < size; sz++) {
		const ptrdiff_t ind_op_r = compute_index_sub_container(order,0,extents,counter);
		const struct const_Matrix_d* op_r_M_ptr = op_r->data[ind_op_r]->op_std;
		if (op_r_M_ptr == NULL) {
			increment_counter_MaO(order,extents,counter);
			continue;
		}

		const int order_l               = op_l->order;
		const ptrdiff_t*const extents_l = op_l->extents;
		assert((order_l == order) || (order_l == order-1));

		const int ind_c = order-order_l;
		const ptrdiff_t ind_op_l = compute_index_sub_container(order_l,0,extents_l,&counter[ind_c]);

		const struct const_Matrix_d* op_l_M_ptr = op_l->data[ind_op_l]->op_std;
		assert(op_l_M_ptr != NULL);

		const struct const_Matrix_d*const op_l_M = constructor_copy_const_Matrix_d(op_l_M_ptr); // destructed
		const struct const_Matrix_d*const op_r_M = constructor_copy_const_Matrix_d(op_r_M_ptr); // destructed
		transpose_Matrix_d((struct Matrix_d*)op_l_M,true);
		transpose_Matrix_d((struct Matrix_d*)op_r_M,true);
		assert((op_l_M->layout == 'C') && (op_r_M->layout == 'C'));

		const ptrdiff_t ext_0_r = op_r_M->ext_0,
		                ext_1_r = op_r_M->ext_1,
		                ext_0_l = op_l_M->ext_0,
		                ext_1_l = op_l_M->ext_1;
		assert(ext_0_r == op_l_M->ext_0);

		struct const_Vector_d op_l_diag = { .ext_0 = ext_0_l, .owns_data = false, .data = NULL, };

		struct Matrix_d* op_lr_local = constructor_zero_Matrix_d('C',ext_0_r,ext_1_r);        // destructed
		struct Matrix_d* op_lr       = constructor_empty_Matrix_d('C',ext_0_r,ext_1_r*ext_1_l); // moved
		for (int j = 0; j < ext_1_l; ++j) {
			const_cast_d1(&op_l_diag.data,get_col_const_Matrix_d(j,op_l_M));
			mm_diag_d('L',1.0,0.0,op_r_M,&op_l_diag,op_lr_local,false);
			set_block_Matrix_d(op_lr,0,j*ext_1_r,
			                   (struct const_Matrix_d*)op_lr_local,0,0,op_lr_local->ext_0,op_lr_local->ext_1,'i');
		}
		destructor_Matrix_d(op_lr_local);
		destructor_const_Matrix_d(op_l_M);
		destructor_const_Matrix_d(op_r_M);

		const_constructor_move_const_Matrix_d(&op->data[ind_op_r]->op_std,(struct const_Matrix_d*)op_lr); // rtrnd
		increment_counter_MaO(order,extents,counter);
	}

	return op;
}

struct Operator_Info* constructor_Operator_Info
	(const char*const name_type, const char*const name_in, const char*const name_out, const char*const name_range,
	 const struct const_Element* element, const struct Simulation* sim)
{
	struct Operator_Info* op_info = malloc(sizeof *op_info); // returned

	op_info->element = element;

	const_cast_i(&op_info->op_type,convert_to_type(name_type));
	const_cast_b(&op_info->transpose,check_for_transpose(name_type,op_info->op_type));

	const_cast_c(&op_info->op_io[OP_IND_I].ce,  name_in[0]);
	const_cast_c(&op_info->op_io[OP_IND_I].kind,name_in[1]);
	const_cast_c(&op_info->op_io[OP_IND_I].sc,  name_in[2]);

	const_cast_c(&op_info->op_io[OP_IND_O].ce,  name_out[0]);
	const_cast_c(&op_info->op_io[OP_IND_O].kind,name_out[1]);
	const_cast_c(&op_info->op_io[OP_IND_O].sc,  name_out[2]);

	const int ranges[] =
		{ convert_to_range_d(name_type),
		  convert_to_range_ce(name_in[0],name_out[0]),
		  convert_to_range('h',name_range),
		  convert_to_range('p',name_range), };

	const_cast_i(&op_info->range_d, ranges[0]);
	const_cast_i(&op_info->range_ce,ranges[1]);
	const_cast_i(&op_info->range_h, ranges[2]);
	const_cast_i(&op_info->range_p, ranges[3]);

	const_cast_i_n(op_info->p_ref,sim->p_ref,2);
	set_up_extents(op_info);
	set_up_values_op(op_info);

	switch (op_info->op_type) {
	case OP_T_CV: case OP_T_CC: case OP_T_VV: case OP_T_VC:
		op_info->set_operator = set_operator_std;
		break;
	case OP_T_TW:
		op_info->set_operator = set_operator_solver;
		break;
	case OP_T_UNUSED:
		; // Do nothing
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}

	return op_info;
}

void destructor_Operator_Info (struct Operator_Info* op_info)
{
	destructor_const_Vector_i(op_info->extents_op);
	destructor_const_Matrix_i(op_info->values_op);
	free(op_info);
}

int compute_p_basis (const struct Op_IO* op_io, const struct Simulation* sim)
{
	const int s_type     = op_io->s_type,
	          p_op       = op_io->p_op,
	          nodes_kind = op_io->kind,
	          op_sc      = op_io->sc;
	switch (nodes_kind) {
	case 's': // solution
		return p_op+sim->p_s_v_p;
		break;
	case 'f': // flux
		return p_op+sim->p_s_f_p;
		break;
	case 'r': // g'r'adient
		return p_op+sim->p_sg_v_p;
		break;
	case 'g': // geometry
		if (op_sc == 's')
			return 1;

		const_cast_c(&op_io->kind,'s');
		const int p_s = compute_p_basis(op_io,sim);
		const_cast_c(&op_io->kind,'g');

		if (strcmp(sim->geom_rep,"isoparametric") == 0)
			return GSL_MAX(p_s,1);
		else if (strcmp(sim->geom_rep,"superparametric") == 0)
			return p_s+1;
		else if (strcmp(sim->geom_rep,"superparametric_p_le_1") == 0)
			return ( p_s < 2 ? p_s+1 : p_s );
		else if (strcmp(sim->geom_rep,"superparametric2") == 0)
			return p_s+2;
		else if (strstr(sim->geom_rep,"fixed"))
			EXIT_ADD_SUPPORT; // Find number in geom_rep (use something similar to 'convert_to_range_d').
		else
			EXIT_ERROR("Unsupported: %s\n",sim->geom_rep);
		break;
	case 'm': { // metric
		const_cast_c(&op_io->kind,'g');
		const int p_g = compute_p_basis(op_io,sim);
		const_cast_c(&op_io->kind,'m');

		/// \todo Revisit this based on the results of the geometric conservation testing.
		switch (s_type) {
		case ST_TP: // fallthrough
		case ST_SI:
			return p_g;
			break;
		default:
			EXIT_ERROR("Unsupported: %d\n",s_type);
			break;
		}
		break;
	} case 'p': // plotting
		return GSL_MAX(1,p_op);
		break;
	case 't': { // test space
		const_cast_c(&op_io->kind,'s');
		const int p_s = compute_p_basis(op_io,sim);
		const_cast_c(&op_io->kind,'t');

		const int ind_p_t = ( op_sc == 's' ? 0 : 1 );
		return p_s + sim->p_t_p[ind_p_t];
		break;
	} case 'v': // vertex
		assert(p_op == 1 || p_op == 2);
		return p_op;
		break;
	case 'c': // fallthrough
	default:
		EXIT_ERROR("Unsupported: %c\n",nodes_kind);
		break;
	}
	EXIT_UNSUPPORTED;
	return -1;
}

const struct const_Vector_i* constructor_indices_Vector_i
	(const int ext_0_expected, const int* op_values, const bool*const indices_skip)
{
	struct Vector_i* indices = ( ext_0_expected == -1 ? constructor_empty_Vector_i(OP_ORDER_MAX)
	                                                  : constructor_empty_Vector_i(ext_0_expected) ); // returned
	int* data = indices->data;

	int ind = 0;
	for (int i = 0; i < OP_ORDER_MAX; ++i) {
		if (!(indices_skip && indices_skip[i]) && op_values[i] != OP_INVALID_IND)
			data[ind++] = op_values[i];
	}

	if (ext_0_expected == -1)
		indices->ext_0 = ind;
	else
		assert(ind == ext_0_expected);

	return (const struct const_Vector_i*)indices;
}

bool op_should_use_L2 (const int*const op_values, const struct Op_IO* op_io)
{
	const char ce_i = op_io[OP_IND_I].ce,
	           ce_o = op_io[OP_IND_O].ce;
	// Don't use L2 for operators having projected nodes.
	if (!((ce_i == ce_o) || (ce_i == 'v'))) // ((vv || ff || ee) || (vf || ve))
		return false;

	if ((op_values[OP_IND_H+OP_IND_I] > op_values[OP_IND_H+OP_IND_O]) ||
	    (op_values[OP_IND_P+OP_IND_I] > op_values[OP_IND_P+OP_IND_O])) {
		switch (op_io[OP_IND_O].kind) {
		case 'c':
		case 'v':
		case 'g':
			// Do nothing.
			break;
		case 's': // fallthrough
		case 'f': // fallthrough
		case 'r': // fallthrough
		case 'p': // fallthrough
			return true;
			break;
		default:
			EXIT_ERROR("Unsupported: %c\n",op_io[OP_IND_O].kind);
			break;
		}
	}
	return false;
}

int compute_super_type_op (const char ce, const int h_op, const struct const_Element* element)
{
	const int sub_e_type = compute_elem_type_sub_ce(element->type,ce,h_op);
	return compute_super_from_elem_type(sub_e_type);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Return the number of operators to be set depending on \ref Operator_Info::range_d.
 *  \return See brief. */
static int get_n_op
	(const int range_d, ///< \ref Operator_Info::range_d.
	 const int d        ///< The dimension of the element for which operators are being constructed.
	);

/** \brief Constructor for a \ref Multiarray_Matrix_T\* of order 1 with `owns_data = false` to hold locally computed
 *         operators.
 *  \return See brief. */
static const struct const_Multiarray_Matrix_d* constructor_op_MMd
	(const bool owns_data, ///< Defined in \ref Multiarray_Matrix_T.
	 const ptrdiff_t ext_0 ///< The size of the single extent.
	);

/// \brief Version of \ref set_operator_std for interpolation operators.
static void set_operator_std_interp
	(const ptrdiff_t ind_values,                    ///< See brief.
	 const struct Operator_Info* op_info,           ///< See brief.
	 const struct Simulation* sim,                  ///< See brief.
	 const struct const_Multiarray_Matrix_d* op_ioN ///< Will hold the matrices of the standard operators.
	);

/// \brief Version of \ref set_operator_std for L2 projection operators.
static void set_operator_std_L2
	(const ptrdiff_t ind_values,                    ///< See brief.
	 const struct Operator_Info* op_info,           ///< See brief.
	 const struct Simulation* sim,                  ///< See brief.
	 const struct const_Multiarray_Matrix_d* op_ioN ///< Will hold the matrices of the standard operators.
	);

/** \brief Get the value for the maximum number of h-refinements.
 *  \return See brief. */
static int get_n_ref_max
	(const struct Operator_Info* op_info ///< \ref Operator_Info.
	);

/// \brief Compute range (min, max) for the index of `var_type`.
static void compute_range
	(int x_mm[2],                         ///< Set to the min and max values.
	 const struct Operator_Info* op_info, ///< \ref Operator_Info.
	 const char var_type,                 ///< Variable type. Options: 'd', 'c'e, 'h', 'p'.
	 const char var_io                    ///< Whether the variable is of 'i'nput or 'o'utput type.
	);

/// \brief Compute range (min, max) for `p_o` based on the values of `p_i` and `op_info->p_range`.
static void compute_range_p_o
	(int p_o_mm[2],                       ///< Set to the min and max values.
	 const struct Operator_Info* op_info, ///< \ref Operator_Info.
	 const int p_i                        ///< The input order.
	);

/** \brief Return the index of the node correspondence sub multiarray currently being set.
 *  \return See brief. */
static ptrdiff_t get_ind_nc
	(const struct Operator_Info* op_info,       ///< \ref Operator_Info.
	 const int* op_values,                      ///< One line of values of \ref Operator_Info::values_op.
	 const struct const_Multiarray_Vector_i* nc ///< Defined for \ref set_node_correspondence.
	);

/** \brief Return the index of the cubature weight currently being set.
 *  \return See brief. */
static ptrdiff_t get_ind_w
	(const struct Operator_Info* op_info,      ///< \ref Operator_Info.
	 const int* op_values,                     ///< One line of values of \ref Operator_Info::values_op.
	 const struct const_Multiarray_Vector_d* w ///< Defined for \ref set_weights.
	);

static void set_operator_std
	(ptrdiff_t*const ind_values, const struct Multiarray_Operator* op, const struct Operator_Info* op_info,
	 const struct Simulation* sim)
{
	const int* op_values = get_row_const_Matrix_i(*ind_values,op_info->values_op);
	const int n_op       = get_n_op(op_info->range_d,op_info->element->d);
	const struct const_Multiarray_Matrix_d* op_ioN = constructor_op_MMd(false,n_op); // destructed

	if (!op_should_use_L2(op_values,op_info->op_io))
		set_operator_std_interp(*ind_values,op_info,sim,op_ioN);
	else
		set_operator_std_L2(*ind_values,op_info,sim,op_ioN);

	ptrdiff_t ind_op = get_ind_op(op_info,op_values,op);
	for (int i = 0; i < n_op; ++i)
		const_constructor_move_const_Matrix_d(&op->data[ind_op++]->op_std,op_ioN->data[i]);
	*ind_values += n_op;

	destructor_const_Multiarray_Matrix_d(op_ioN);
}

static void set_operator_solver
	(ptrdiff_t*const ind_values, const struct Multiarray_Operator* op, const struct Operator_Info* op_info,
	 const struct Simulation* sim)
{
	const int op_type = op_info->op_type;
	assert(op_type == OP_T_TW);    // To be made flexible.

	const int* op_values = get_row_const_Matrix_i(*ind_values,op_info->values_op);

	const_cast_i(&op_info->op_type,OP_T_CV);
	set_operator_std(ind_values,op,op_info,sim);
	const_cast_i(&op_info->op_type,op_type);

	const struct const_Element* element = op_info->element;

	const struct Op_IO* op_io = op_info->op_io;
	const struct const_Nodes* nodes_o = constructor_const_Nodes_h(OP_IND_O,op_io,element,sim); // destructed
	assert(nodes_o->has_weights);

	const int n_op   = get_n_op(op_info->range_d,element->d);
	ptrdiff_t ind_op = get_ind_op(op_info,op_values,op);

	for (int i = 0; i < n_op; ++i) {
		struct Matrix_d* op_std = (struct Matrix_d*) op->data[ind_op+i]->op_std;
		transpose_Matrix_d(op_std,false);
		scale_Matrix_d_by_Vector_d('R',1.0,op_std,nodes_o->w,false);
	}
	destructor_const_Nodes(nodes_o);

	if (sim->collocated) {
		EXIT_ERROR("Ensure that all is working as expected.\n");
		// Note: Test space is higher degree than solution so efficiency advantages would not be present.
		if (sim->method == METHOD_DPG)
			return;

		const struct const_Nodes* nodes_i = constructor_const_Nodes_h(OP_IND_I,op_io,element,sim); // destructed
		assert(nodes_i->has_weights);

		for (int i = 0; i < n_op; ++i) {
			struct Matrix_d* op_std = (struct Matrix_d*) op->data[ind_op+i]->op_std;
			scale_Matrix_d_by_Vector_d('L',1.0,op_std,nodes_i->w,true);
		}
		destructor_const_Nodes(nodes_i);
	}
}

static void set_node_correspondence
	(ptrdiff_t ind_values, const struct const_Multiarray_Vector_i* nc, const struct Operator_Info* op_info,
	 const struct Simulation* sim)
{
	const int* op_values = get_row_const_Matrix_i(ind_values,op_info->values_op);

	set_current_op_io(op_info,op_values);
	const struct Op_IO* op_io = op_info->op_io;

	const struct const_Element* element = op_info->element;
	const struct const_Multiarray_Vector_i* nc_curr =
		constructor_nodes_face_corr_op(&op_io[OP_IND_O],element,sim); // destructed/moved

	const ptrdiff_t n_nc = nc->extents[0];
	ptrdiff_t ind_nc = get_ind_nc(op_info,op_values,nc);
	for (ptrdiff_t i = 0; i < n_nc; ++i)
		const_constructor_move_const_Vector_i(&nc->data[ind_nc++],nc_curr->data[i]);

	const_cast_b(&nc_curr->owns_data,false);
	destructor_const_Multiarray_Vector_i(nc_curr);
}

static void set_weights
	(ptrdiff_t ind_values, const struct const_Multiarray_Vector_d* w, const struct Operator_Info* op_info,
	 const struct Simulation* sim)
{
	const int* op_values = get_row_const_Matrix_i(ind_values,op_info->values_op);

	set_current_op_io(op_info,op_values);
	const struct Op_IO* op_io = op_info->op_io;

	const struct const_Element* element = op_info->element;
	const struct const_Vector_d* w_curr = constructor_weights(&op_io[OP_IND_O],element,sim); // moved

	ptrdiff_t ind_w = get_ind_w(op_info,op_values,w);
	const_constructor_move_const_Vector_d(&w->data[ind_w++],w_curr);
}

static int convert_to_type (const char* name_type)
{
	if (strstr(name_type,"cv"))
		return OP_T_CV;
	else if (strstr(name_type,"cc"))
		return OP_T_CC;
	else if (strstr(name_type,"vv"))
		return OP_T_VV;
	else if (strstr(name_type,"vc"))
		return OP_T_VC;
	else if (strstr(name_type,"tw"))
		return OP_T_TW;
	else if (strstr(name_type,"UNUSED"))
		return OP_T_UNUSED;
	else
		EXIT_ERROR("(%s)\n",name_type);
	return -1;
}

int check_for_transpose (const char* name_type, const int op_type)
{
	switch (op_type) {
	case OP_T_CV: // fallthrough
	case OP_T_CC: // fallthrough
	case OP_T_VV: // fallthrough
	case OP_T_VC:
		if (name_type[2] == 't')
			return true;
		break;
	case OP_T_TW:     // fallthrough
	case OP_T_UNUSED:
		break; // Do nothing.
	default:
		EXIT_ERROR("Unsupported: %d.\n",op_type);
		break;
	}
	return false;
}

static int convert_to_range_d (const char* name_type)
{
	char* char_end = NULL;
	long order_d = -1;
	while (*name_type) {
		if (isdigit(*name_type)) {
			order_d = strtol(name_type,&char_end,10);
			break;
		}
		++name_type;
	}

	switch (order_d) {
		case 0:  return OP_R_D_0;   break;
		case 1:  return OP_R_D_ALL; break;
		default: EXIT_UNSUPPORTED;  break;
	}
	EXIT_ERROR("Did not find the operator range.");
	return -1;
}

static int convert_to_range_ce (const char ce_i, const char ce_o)
{
	switch (ce_i) {
	case 'v':
		switch (ce_o) {
			case 'v': return OP_R_CE_VV; break;
			case 'f': return OP_R_CE_VF; break;
			case 'e': return OP_R_CE_VE; break;
			default:  EXIT_UNSUPPORTED;  break;
		}
		break;
	case 'f':
		switch (ce_o) {
			case 'v': return OP_R_CE_FV; break;
			case 'f': return OP_R_CE_FF; break;
			case 'e': return OP_R_CE_FE; break;
			default:  EXIT_UNSUPPORTED;  break;
		}
		break;
	case 'e':
		switch (ce_o) {
			case 'v': return OP_R_CE_EV; break;
			case 'f': return OP_R_CE_EF; break;
			case 'e': return OP_R_CE_EE; break;
			default:  EXIT_UNSUPPORTED;  break;
		}
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}
	EXIT_ERROR("Did not find the operator range.");
	return -1;
}

static int convert_to_range (const char type_range, const char*const name_range)
{
	switch (type_range) {
	case 'h':
		if (strstr(name_range,"H_1"))
			return OP_R_H_1;
		else if (strstr(name_range,"H_CF"))
			return OP_R_H_CF;
		else if (strstr(name_range,"H_FC"))
			return OP_R_H_FC;
		else if (strstr(name_range,"H_ALL"))
			return OP_R_H_ALL;
		else
			EXIT_UNSUPPORTED;
		break;
	case 'p':
		if (strstr(name_range,"P_1P"))
			if (strstr(name_range,"P_1PPM1"))
				return OP_R_P_1PPM1;
			else
				return OP_R_P_1P;
		else if (strstr(name_range,"P_1"))
			if (strstr(name_range,"P_12"))
				return OP_R_P_12;
			else
				return OP_R_P_1;
		else if (strstr(name_range,"P_P1"))
			return OP_R_P_P1;
		else if (strstr(name_range,"P_PM0"))
			return OP_R_P_PM0;
		else if (strstr(name_range,"P_PM1"))
			return OP_R_P_PM1;
		else if (strstr(name_range,"P_ALL"))
			return OP_R_P_ALL;
		else
			EXIT_ERROR("Unsupported: %s\n",name_range);
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}
	EXIT_ERROR("Did not find the operator range.");
	return -1;
}

static void set_up_extents (struct Operator_Info* op_info)
{
	const struct const_Element* element = op_info->element;

	struct Vector_i* extents_op = constructor_default_Vector_i(); // keep

	switch (op_info->range_d) {
		case OP_R_D_0:
			// Do nothing
			break;
		case OP_R_D_ALL:
			push_back_Vector_i(extents_op,element->d,false,false);
			break;
		default:
			EXIT_UNSUPPORTED;
			break;
	}

	switch (op_info->range_ce) {
		case OP_R_CE_VV: // fallthrough
		case OP_R_CE_EE:
			// Do nothing
			break;
		case OP_R_CE_FF: {
			const int n_fe = get_number_of_face_elements(op_info->element);
			switch (n_fe) {
			case 1: // fallthrough
			case 2:
				push_back_Vector_i(extents_op,n_fe,false,false);
				push_back_Vector_i(extents_op,n_fe,false,false);
				break;
			default:
				EXIT_ERROR("Unsupported: %d\n",n_fe);
				break;
			}
			break;
		}
		case OP_R_CE_VF:
		case OP_R_CE_FV:
			push_back_Vector_i(extents_op,element->n_f,false,false);
			break;
		case OP_R_CE_VE:
		case OP_R_CE_EV:
			push_back_Vector_i(extents_op,element->n_e,false,false);
			break;
		case OP_R_CE_FE:
		case OP_R_CE_EF:
			// Note: Potential special cases for wedge/pyr elements with faces having varying number of edges.
			EXIT_ADD_SUPPORT;
			break;
		default:
			EXIT_UNSUPPORTED;
			break;
	}

	switch (op_info->range_h) { // h_o, h_i
		case OP_R_H_1:
			push_back_Vector_i(extents_op,1,false,false);
			push_back_Vector_i(extents_op,1,false,false);
			break;
		case OP_R_H_CF: {
			const int n_ref_max = get_n_ref_max(op_info);
			push_back_Vector_i(extents_op,n_ref_max,false,false);
			push_back_Vector_i(extents_op,1,false,false);
			break;
		} case OP_R_H_FC: {
			const int n_ref_max = get_n_ref_max(op_info);
			push_back_Vector_i(extents_op,1,false,false);
			push_back_Vector_i(extents_op,n_ref_max,false,false);
			break;
		} case OP_R_H_ALL: {
			const int n_ref_max = get_n_ref_max(op_info);
			push_back_Vector_i(extents_op,n_ref_max,false,false);
			push_back_Vector_i(extents_op,n_ref_max,false,false);
			break;
		} default:
			EXIT_UNSUPPORTED;
			break;
	}

	switch (op_info->range_p) { // p_o, p_i
		case OP_R_P_1:
			push_back_Vector_i(extents_op,2,false,false);
			push_back_Vector_i(extents_op,2,false,false);
			break;
		case OP_R_P_12:
			push_back_Vector_i(extents_op,3,false,false);
			push_back_Vector_i(extents_op,3,false,false);
			break;
		case OP_R_P_1P:    // fallthrough
		case OP_R_P_1PPM1:
			push_back_Vector_i(extents_op,op_info->p_ref[1]+1,false,false);
			push_back_Vector_i(extents_op,2,false,false);
			break;
		case OP_R_P_P1:
			push_back_Vector_i(extents_op,2,false,false);
			push_back_Vector_i(extents_op,op_info->p_ref[1]+1,false,false);
			break;
		case OP_R_P_PM0:   // fallthrough
		case OP_R_P_PM1:   // fallthrough
		case OP_R_P_ALL:
			push_back_Vector_i(extents_op,op_info->p_ref[1]+1,false,false);
			push_back_Vector_i(extents_op,op_info->p_ref[1]+1,false,false);
			break;
		default:
			EXIT_UNSUPPORTED;
			break;
	}

	op_info->extents_op  = (const struct const_Vector_i*) extents_op;
}

static void set_up_values_op (struct Operator_Info* op_info)
{
	const ptrdiff_t size = prod_const_Vector_i(op_info->extents_op);
	assert(size != 0);

	struct Matrix_i* values = constructor_empty_Matrix_i('R',size,OP_ORDER_MAX); // keep

	int d_mm[2],
	    ce_o_mm[2],
	    ce_i_mm[2],
	    h_o_mm[2],
	    h_i_mm[2],
	    p_o_mm[2],
	    p_i_mm[2];

	compute_range(d_mm,op_info,'d','i');
	compute_range(ce_o_mm,op_info,'c','o');
	compute_range(ce_i_mm,op_info,'c','i');
	compute_range(h_o_mm,op_info,'h','o');
	compute_range(h_i_mm,op_info,'h','i');
	compute_range(p_i_mm,op_info,'p','i');

	int row = 0;
	for (int p_i = p_i_mm[0]; p_i < p_i_mm[1]; ++p_i) {
		compute_range_p_o(p_o_mm,op_info,p_i);
		for (int p_o = p_o_mm[0]; p_o < p_o_mm[1]; ++p_o) {
		for (int h_i = h_i_mm[0]; h_i < h_i_mm[1]; ++h_i) {
		for (int h_o = h_o_mm[0]; h_o < h_o_mm[1]; ++h_o) {
			if (!(h_o == 0 || h_i == 0))
				continue;
			for (int ce_i = ce_i_mm[0]; ce_i < ce_i_mm[1]; ++ce_i) {
			for (int ce_o = ce_o_mm[0]; ce_o < ce_o_mm[1]; ++ce_o) {
			for (int d = d_mm[0]; d < d_mm[1]; ++d) {
				set_row_Matrix_i(row,values,(int[]){d,ce_o,ce_i,h_o,h_i,p_o,p_i});
				++row;
			}}}
		}}}
	}
	values->ext_0 = row;

	op_info->values_op = (const struct const_Matrix_i*) values;
}

static void transpose_operators (struct Multiarray_Operator* op)
{
	const ptrdiff_t size = compute_size(op->order,op->extents);
	for (ptrdiff_t i = 0; i < size; ++i) {
		struct Matrix_d* op_std = (struct Matrix_d*) op->data[i]->op_std;
		if (!op_std)
			continue;

		transpose_Matrix_d(op_std,false);
	}
}

static const int* get_basis_types (const char* name_type, const struct Simulation* sim)
{
	static int basis_type[2] = { -1, -1, };
	if (strstr(name_type,"SB")) {
		basis_type[0] = get_basis_i_from_s(sim->basis_sol);
		basis_type[1] = BASIS_BEZIER;
	} else if (strstr(name_type,"BS")) {
		basis_type[0] = BASIS_BEZIER;
		basis_type[1] = get_basis_i_from_s(sim->basis_sol);
	} else {
		EXIT_ADD_SUPPORT;
	}
	return basis_type;
}

static ptrdiff_t get_ind_op
	(const struct Operator_Info* op_info, const int* op_values, const struct Multiarray_Operator* op)
{
	const int order_op = (int)op_info->extents_op->ext_0;
	const struct const_Vector_i* indices_op = constructor_indices_Vector_i(order_op,op_values,NULL); // destructed

	ptrdiff_t ind_op = compute_index_sub_container_pi(op->order,0,op->extents,indices_op->data);
	destructor_const_Vector_i(indices_op);

	return ind_op;
}

static void set_current_op_io (const struct Operator_Info* op_info, const int* op_values)
{
	const struct Op_IO* op_io = op_info->op_io;

	const int*const ce_ptr = &op_values[OP_IND_CE],
	         *const h_ptr = &op_values[OP_IND_H],
	         *const p_ptr = &op_values[OP_IND_P];
	for (int i = 0; i < 2; ++i) {
		const_cast_i(&op_io[i].ce_op,ce_ptr[i]);
		const_cast_i(&op_io[i].h_op,h_ptr[i]);
		const_cast_i(&op_io[i].p_op,p_ptr[i]);
		const_cast_i(&op_io[i].s_type,compute_super_type_op(op_io[i].ce,op_io[i].h_op,op_info->element));
	}
}

static const struct const_Matrix_d* constructor_cv
	(const struct const_Matrix_d* rst, const struct Op_IO* op_io, const struct Simulation* sim,
	 const int basis_type_i)
{
	int basis_type = -1;
	switch (basis_type_i) {
	case BASIS_ORTHO: case BASIS_LAGRANGE: case BASIS_BEZIER:
		basis_type = basis_type_i;
		break;
	default:
		switch (op_io->kind) {
		case 'g': // fallthrough
		case 'm':
			basis_type = get_basis_i_from_s(sim->basis_geom);
			break;
		case 's': // fallthrough
		case 'f': // fallthrough
		case 'r': // fallthrough
		case 't':
			basis_type = get_basis_i_from_s(sim->basis_sol);
			break;
		case 'p': // fallthrough
		case 'v':
			basis_type = get_basis_i_from_s("lagrange");
			break;
		default:
			EXIT_ERROR("Unsupported: %c.\n",op_io->kind);
			break;
		}
		break;
	}

	const struct const_Matrix_d* cv = NULL;
	if (basis_type == BASIS_LAGRANGE) {
		cv = constructor_identity_const_Matrix_d('R',rst->ext_0); // returned
	} else {
		const int s_type = op_io->s_type;
		const int p_basis = compute_p_basis(op_io,sim);

		constructor_basis_fptr constructor_basis = get_constructor_basis_by_super_type_i(s_type,basis_type);
		cv = constructor_basis(p_basis,rst); // returned
	}

	return cv;
}

// Level 1 ********************************************************************************************************** //

/** \brief Constructor for the \*c0_\*\*\*_\*\*\* operator from input to output computational element.
 *  \return The required operator if required; `NULL` otherwise. */
static const struct const_Multiarray_Matrix_d* constructor_op_Xc_ce_io
	(const ptrdiff_t ind_values,          ///< Defined for \ref set_operator_std_L2.
	 const struct Operator_Info* op_info, ///< Defined for \ref set_operator_std_L2.
	 const struct Simulation* sim         ///< Defined for \ref set_operator_std_L2.
	);

/** \brief Construct the mass-type \ref const_Matrix_T\* with the given input operators.
 *  \return See brief. */
static const struct const_Matrix_d* constructor_mass
	(const struct const_Matrix_d* cv0_l, ///< The 'l'eft  cv matrix (not transposed).
	 const struct const_Matrix_d* cv0_r, ///< The 'r'ight cv matrix.
	 const struct const_Vector_d* w,     ///< The cubature weights.
	 const double volume_ratio           ///< See \ref compute_volume_ratio.
	);

/** \brief Get the `char` representation of the output of \ref Operator_Info::op_type.
 *  \return 'c'oefficient or 'v'alue. */
static char get_op_type_out
	(const int op_type ///< \ref Operator_Info::op_type.
	);

/** \brief Compute the volume relative to the complete reference element.
 *  \return See brief.
 *
 *  When nodes are given for a sub-region of the reference element, this ratio should be used to scale any integrals
 *  computed as if the complete reference element was being used.
 */
static double compute_volume_ratio
	(const int e_type, ///< \ref Element::type.
	 const int ind_h   ///< The h-refinement index.
	);

static int get_n_op (const int range_d, const int d)
{
	return ( range_d == OP_R_D_0 ? 1 : d );
}

static const struct const_Multiarray_Matrix_d* constructor_op_MMd (const bool owns_data, const ptrdiff_t ext_0)
{
	struct Multiarray_Matrix_d* op_MMd = constructor_empty_Multiarray_Matrix_d(false,1,&ext_0); // returned
	op_MMd->owns_data = owns_data;

	return (const struct const_Multiarray_Matrix_d*) op_MMd;
}

static void set_operator_std_interp
	(const ptrdiff_t ind_values, const struct Operator_Info* op_info, const struct Simulation* sim,
	 const struct const_Multiarray_Matrix_d* op_ioN)
{
	const int* op_values = get_row_const_Matrix_i(ind_values,op_info->values_op);
	const int n_op = (int)compute_size(op_ioN->order,op_ioN->extents);

	const struct const_Element* element = op_info->element;

	set_current_op_io(op_info,op_values);
	const struct Op_IO* op_io = op_info->op_io;

	const int p_i = compute_p_basis(&op_io[OP_IND_I],sim);

	const int s_type = op_io[OP_IND_I].s_type;
	constructor_basis_fptr constructor_basis = get_constructor_basis_by_super_type(s_type,"orthonormal");

	/* Compute `op_cvNr`: [operator, coefficients to values, differentiation order N, reference basis].
	 * op_cvNr == basis(p_i,rst_o). */
	const struct const_Nodes* nodes_io = constructor_const_Nodes_h(OP_IND_O,op_io,element,sim); // destructed

	const struct const_Multiarray_Matrix_d* op_cvNr = NULL;
	if (op_info->range_d == OP_R_D_0) {
		const struct const_Matrix_d* cv0r = constructor_basis(p_i,nodes_io->rst); // moved

		op_cvNr = constructor_op_MMd(true,n_op); // destructed
		const_constructor_move_const_Matrix_d(&op_cvNr->data[0],cv0r); // destructed
	} else if (op_info->range_d == OP_R_D_ALL) {
		constructor_grad_basis_fptr constructor_grad_basis =
			get_constructor_grad_basis_by_super_type(s_type,"orthonormal");
		op_cvNr = constructor_grad_basis(p_i,nodes_io->rst); // destructed
	} else {
		EXIT_ERROR("Unsupported: %d\n",op_info->range_d);
	}
	destructor_const_Nodes(nodes_io);

	/* Compute `op_cvN`: [operator, coefficients to values, differentiation order N].
	 * op_cvN == op_cvNr*T_i (Corollary 2.2, \cite Zwanenburg2016). */
	const struct const_Nodes* nodes_ii = constructor_const_Nodes_h(OP_IND_I,op_io,element,sim); // destructed

	const struct const_Matrix_d* cv0r_ii = constructor_basis(p_i,nodes_ii->rst);                 // destructed
	const struct const_Matrix_d* cv0_ii  = constructor_cv(nodes_ii->rst,&op_io[OP_IND_I],sim,0); // destructed
	destructor_const_Nodes(nodes_ii);

	const struct const_Matrix_d* T_i = constructor_sgesv_const_Matrix_d(cv0r_ii,cv0_ii); // destructed

	destructor_const_Matrix_d(cv0r_ii);

	const struct const_Multiarray_Matrix_d* op_cvN = constructor_op_MMd(true,n_op);  // destructed
	for (int i = 0; i < n_op; ++i) {
		const struct const_Matrix_d* cv0 = constructor_mm_NN1R_const_Matrix_d(op_cvNr->data[i],T_i); // moved
		const_constructor_move_const_Matrix_d(&op_cvN->data[i],cv0); // destructed
	}
	destructor_const_Matrix_d(T_i);
	destructor_const_Multiarray_Matrix_d(op_cvNr);

	const int op_type = op_info->op_type;

	/* Compute `op_coN`: [operator, coefficients to output (coefs/values), differentiation order N].
	 * op_coN == inv_cv0_oo*op_cvN   <=>   op_coN = sgesv(cv0_oo,op_cvN). */
	const struct const_Multiarray_Matrix_d* op_coN = constructor_op_MMd(true,n_op);  // destructed
	if (op_type == OP_T_CC || op_type == OP_T_VC) {
		const struct Op_IO op_io_oo[2] = { op_io[OP_IND_O], op_io[OP_IND_O], };
		const struct const_Nodes* nodes_oo = constructor_const_Nodes_h(OP_IND_O,op_io_oo,element,sim); // destructed

		const struct const_Matrix_d* cv0_oo = constructor_cv(nodes_oo->rst,&op_io[OP_IND_O],sim,0); // destructed
		destructor_const_Nodes(nodes_oo);

		for (int i = 0; i < n_op; ++i) {
			const struct const_Matrix_d* cv0 = op_cvN->data[i];
			const struct const_Matrix_d* co0 = constructor_sgesv_const_Matrix_d(cv0_oo,cv0); // moved
			const_constructor_move_const_Matrix_d(&op_coN->data[i],co0); // destructed
		}
		destructor_const_Matrix_d(cv0_oo);
	} else {
		const_cast_b(&op_cvN->owns_data,false);
		for (int i = 0; i < n_op; ++i)
			const_constructor_move_const_Matrix_d(&op_coN->data[i],op_cvN->data[i]); // destructed
	}
	destructor_const_Multiarray_Matrix_d(op_cvN);

	/* Compute `op_ioN`: [operator, input (coefs/values) to output (coefs/values), differentiation order N].
	 * op_ioN == op_coN*inv_cv0_ii. */
	if (op_type == OP_T_VV || op_type == OP_T_VC) {
		const struct const_Matrix_d* inv_cv0_ii = constructor_inverse_const_Matrix_d(cv0_ii); // destructed
		for (int i = 0; i < n_op; ++i) {
			const struct const_Matrix_d* co0 = op_coN->data[i];
			const struct const_Matrix_d* io0 = constructor_mm_NN1R_const_Matrix_d(co0,inv_cv0_ii); // moved
			const_constructor_move_const_Matrix_d(&op_ioN->data[i],io0); // returned
		}
		destructor_const_Matrix_d(inv_cv0_ii);
	} else {
		const_cast_b(&op_coN->owns_data,false);
		for (int i = 0; i < n_op; ++i)
			const_constructor_move_const_Matrix_d(&op_ioN->data[i],op_coN->data[i]); // returned
	}
	destructor_const_Matrix_d(cv0_ii);
	destructor_const_Multiarray_Matrix_d(op_coN);
}

static void set_operator_std_L2
	(const ptrdiff_t ind_values, const struct Operator_Info* op_info, const struct Simulation* sim,
	 const struct const_Multiarray_Matrix_d* op_ioN)
{
	const struct Op_IO* op_io = op_info->op_io;

	const char ce_i = op_io[OP_IND_I].ce,
	           ce_o = op_io[OP_IND_O].ce;
	// Add support if required.
	assert((ce_i == 'v' && ce_o == 'v') || (ce_i == 'v' && ce_o == 'f') || (ce_i == 'v' && ce_o == 'e') ||
	       (ce_i == 'f' && ce_o == 'f'));

	const int n_op = (int)compute_size(op_ioN->order,op_ioN->extents);
	assert(n_op == 1); //< Should not have any differentiation operators here.

	const int* op_values = get_row_const_Matrix_i(ind_values,op_info->values_op);

	const int p_o = op_values[OP_IND_P+OP_IND_O],
	          p_i = op_values[OP_IND_P+OP_IND_I];

	assert((OP_IND_O == 0) && (OP_IND_I == 1)); // Otherwise, swap op_io_*[0] and op_io_*[1] below.
	const struct const_Element* e_op = op_info->element;

	const int ind_ce = op_values[OP_IND_CE+OP_IND_O];

	const struct const_Element* el = NULL;
	switch (ce_o) {
		case 'v': el = e_op; break;
		case 'f': el = get_element_by_face(e_op,ind_ce); break;
		case 'e': el = get_element_by_edge(e_op,ind_ce); break;
		default:  EXIT_ERROR("Unsupported: %c\n",ce_o); break;
	}

	const int p_op_c   = GSL_MAX(p_o,p_i),
	          s_type_o = el->s_type;

	constructor_basis_fptr constructor_basis = get_constructor_basis_by_super_type(s_type_o,"orthonormal");

	/* Note: The value of "sc" for the input must be used for both coarse and fine nodes such that the same cubature
	 *       order is selected for both sets of operators. */
	const char sc = op_io[OP_IND_I].sc;

	// 'C'oarse
	struct Op_IO op_io_C[2] = {
		{.ce = 'v', .kind = 'c', .sc = sc, .h_op = 0, .p_op = p_op_c, .s_type = s_type_o, },
		{.ce = 'v', .s_type = s_type_o, },
	};
	const struct const_Nodes* nodes_Cc = constructor_const_Nodes_h(OP_IND_O,op_io_C,el,sim); // destructed

	const_cast_i(&op_io_C[OP_IND_O].p_op,p_o);
	const_cast_c(&op_io_C[0].kind,op_io[OP_IND_O].kind);
	const struct const_Nodes* nodes_Cb = constructor_const_Nodes_h(OP_IND_O,op_io_C,el,sim); // destructed

	const int p_Cb = compute_p_basis(&op_io_C[OP_IND_O],sim);
	const struct const_Matrix_d* cv0r_Cb_Cb = constructor_basis(p_Cb,nodes_Cb->rst);                  // dest.
	const struct const_Matrix_d* cv0_Cb_Cb  = constructor_cv(nodes_Cb->rst,&op_io_C[OP_IND_O],sim,0); // dest.
	const struct const_Matrix_d* T_C        = constructor_sgesv_const_Matrix_d(cv0r_Cb_Cb,cv0_Cb_Cb); // dest.
	destructor_const_Matrix_d(cv0r_Cb_Cb);
	destructor_const_Nodes(nodes_Cb);

	const struct const_Matrix_d* cv0r_Cb_Cc = constructor_basis(p_Cb,nodes_Cc->rst);              // destructed
	const struct const_Matrix_d* cv0_Cb_Cc  = constructor_mm_NN1R_const_Matrix_d(cv0r_Cb_Cc,T_C); // destructed
	destructor_const_Matrix_d(cv0r_Cb_Cc);
	destructor_const_Matrix_d(T_C);

	// Can possibly implement a simplification as mass_C is always identity in the reference basis.
	const struct const_Matrix_d* mass_C = constructor_mass(cv0_Cb_Cc,cv0_Cb_Cc,nodes_Cc->w,1.0); // destructed
	destructor_const_Nodes(nodes_Cc);


	// 'F'ine
	const int h_i = op_values[OP_IND_H+OP_IND_I];
	struct Op_IO op_io_F[2] = {
		{.ce = 'v', .kind = 'c', .sc = sc, .h_op = h_i, .p_op = p_op_c, .s_type = s_type_o, },
		{.ce = 'v', .s_type = s_type_o, },
	};

	const struct const_Nodes* nodes_Fc = constructor_const_Nodes_h(OP_IND_O,op_io_F,el,sim); // destructed

	const_cast_i(&op_io_F[OP_IND_O].p_op,p_i);
	const_cast_c(&op_io_F[0].kind,op_io[OP_IND_I].kind);
	const struct const_Nodes* nodes_Fb = constructor_const_Nodes_h(OP_IND_O,op_io_F,el,sim); // destructed

	const int p_Fb = compute_p_basis(&op_io_F[OP_IND_O],sim);
	const struct const_Matrix_d* cv0r_Fb_Fb = constructor_basis(p_Fb,nodes_Fb->rst);                  // dest.
	const struct const_Matrix_d* cv0_Fb_Fb  = constructor_cv(nodes_Fb->rst,&op_io_F[OP_IND_O],sim,0); // dest.
	const struct const_Matrix_d* T_F        = constructor_sgesv_const_Matrix_d(cv0r_Fb_Fb,cv0_Fb_Fb); // dest.
	destructor_const_Matrix_d(cv0r_Fb_Fb);
	destructor_const_Matrix_d(cv0_Fb_Fb);
	destructor_const_Nodes(nodes_Fb);

	const struct const_Matrix_d* cv0r_Fb_Fc = constructor_basis(p_Fb,nodes_Fc->rst);              // destructed
	const struct const_Matrix_d* cv0_Fb_Fc  = constructor_mm_NN1R_const_Matrix_d(cv0r_Fb_Fc,T_F); // destructed
	destructor_const_Matrix_d(cv0r_Fb_Fc);
	destructor_const_Matrix_d(T_F);

	const double v_ratio = compute_volume_ratio(el->type,h_i);
	const struct const_Matrix_d* mass_F = constructor_mass(cv0_Cb_Cc,cv0_Fb_Fc,nodes_Fc->w,v_ratio); // destructed
	destructor_const_Matrix_d(cv0_Cb_Cc);
	destructor_const_Matrix_d(cv0_Fb_Fc);
	destructor_const_Nodes(nodes_Fc);

	const struct const_Matrix_d* op_cc = constructor_sgesv_const_Matrix_d(mass_C,mass_F); // destructed
	destructor_const_Matrix_d(mass_C);
	destructor_const_Matrix_d(mass_F);

	const struct const_Matrix_d* op_ic = NULL;
	if (ce_i != ce_o) {
		const struct const_Multiarray_Matrix_d* op_Xc_ce_io =
			constructor_op_Xc_ce_io(ind_values,op_info,sim); // destructed

		op_ic = constructor_mm_NN1R_const_Matrix_d(op_cc,op_Xc_ce_io->data[0]); // destructed/returned
		destructor_const_Matrix_d(op_cc);
		destructor_const_Multiarray_Matrix_d(op_Xc_ce_io);
	} else {
		op_ic = op_cc;
	}

	const struct const_Matrix_d* op_io_M = NULL;
	if (get_op_type_out(op_info->op_type) == 'c') {
		op_io_M = op_ic;
	} else {
		op_io_M = constructor_mm_NN1R_const_Matrix_d(cv0_Cb_Cb,op_ic); // moved
		destructor_const_Matrix_d(op_ic);
	}
	destructor_const_Matrix_d(cv0_Cb_Cb);

	const_constructor_move_const_Matrix_d(&op_ioN->data[0],op_io_M); // returned
}

static int get_n_ref_max (const struct Operator_Info* op_info)
{
	// There are currently no h-adaptation operators for other `range_ce` values.
	if (op_info->range_ce == OP_R_CE_VV)
		return op_info->element->n_ref_max_v;
	else if (op_info->range_ce == OP_R_CE_VF || op_info->range_ce == OP_R_CE_FF)
		return op_info->element->n_ref_max_f;
	EXIT_UNSUPPORTED;
}

static void compute_range (int x_mm[2], const struct Operator_Info* op_info, const char var_type, const char var_io)
{
	assert((var_io == 'i' || var_io == 'o'));

	switch (var_type) {
	case 'd':
		switch (op_info->range_d) {
		case OP_R_D_0:
			x_mm[0] = OP_INVALID_IND;
			x_mm[1] = OP_INVALID_IND+1;
			break;
		case OP_R_D_ALL:
			x_mm[0] = 0;
			x_mm[1] = op_info->element->d;
			break;
		default:
			EXIT_UNSUPPORTED;
			break;
		}
		break;
	case 'c': // ce
		switch (op_info->range_ce) {
		case OP_R_CE_VV: // fallthrough
		case OP_R_CE_EE:
			x_mm[0] = OP_INVALID_IND;
			x_mm[1] = OP_INVALID_IND+1;
			break;
		case OP_R_CE_FF:
			x_mm[0] = 0;
			x_mm[1] = get_number_of_face_elements(op_info->element);
			break;
		case OP_R_CE_VF:
			if (var_io == 'i') {
				x_mm[0] = OP_INVALID_IND;
				x_mm[1] = OP_INVALID_IND+1;
			} else if (var_io == 'o') {
				x_mm[0] = 0;
				x_mm[1] = op_info->element->n_f;
			}
			break;
		case OP_R_CE_FV:
			if (var_io == 'i') {
				x_mm[0] = 0;
				x_mm[1] = op_info->element->n_f;
			} else if (var_io == 'o') {
				x_mm[0] = OP_INVALID_IND;
				x_mm[1] = OP_INVALID_IND+1;
			}
			break;
		case OP_R_CE_VE:
			if (var_io == 'i') {
				x_mm[0] = OP_INVALID_IND;
				x_mm[1] = OP_INVALID_IND+1;
			} else if (var_io == 'o') {
				x_mm[0] = 0;
				x_mm[1] = op_info->element->n_e;
			}
			break;
		case OP_R_CE_EV:
			if (var_io == 'i') {
				x_mm[0] = 0;
				x_mm[1] = op_info->element->n_e;
			} else if (var_io == 'o') {
				x_mm[0] = OP_INVALID_IND;
				x_mm[1] = OP_INVALID_IND+1;
			}
			break;
		case OP_R_CE_FE:
		case OP_R_CE_EF:
			EXIT_ADD_SUPPORT;
			break;
		default:
			EXIT_UNSUPPORTED;
			break;
		}
		break;
	case 'h':
		x_mm[0] = 0;
		switch (op_info->range_h) {
		case OP_R_H_1:
			x_mm[1] = 1;
			break;
		case OP_R_H_CF:
			if (var_io == 'i')
				x_mm[1] = 1;
			else if (var_io == 'o')
				x_mm[1] = get_n_ref_max(op_info);
			break;
		case OP_R_H_FC:
			if (var_io == 'i')
				x_mm[1] = get_n_ref_max(op_info);
			else if (var_io == 'o')
				x_mm[1] = 1;
			break;
		case OP_R_H_ALL:
			x_mm[1] = get_n_ref_max(op_info);
			break;
		default:
			EXIT_UNSUPPORTED;
			break;
		}
		break;
	case 'p':
		// In general, p_o depends on p_i and requires more information than that provided here.
		assert(var_io == 'i');

		switch (op_info->range_p) {
		case OP_R_P_1:  // fallthrough
		case OP_R_P_1P: // fallthrough
		case OP_R_P_1PPM1:
			x_mm[0] = 1;
			x_mm[1] = 1+1;
			break;
		case OP_R_P_12:
			x_mm[0] = 1;
			x_mm[1] = 2+1;
			break;
		case OP_R_P_P1:  // fallthrough
		case OP_R_P_PM0: // fallthrough
		case OP_R_P_PM1: // fallthrough
		case OP_R_P_ALL:
			x_mm[0] = op_info->p_ref[0];
			x_mm[1] = op_info->p_ref[1]+1;
			break;
		default:
			EXIT_UNSUPPORTED;
			break;
		}
		break;
	default:
		EXIT_ERROR("Unsupported: %c\n",var_type);
		break;
	}
}

static void compute_range_p_o (int p_o_mm[2], const struct Operator_Info* op_info, const int p_i)
{
	const int*const p_ref = op_info->p_ref;
	switch (op_info->range_p) {
	case OP_R_P_1: // fallthrough
	case OP_R_P_P1:
		p_o_mm[0] = 1;
		p_o_mm[1] = 1+1;
		break;
	case OP_R_P_12:
		p_o_mm[0] = 1;
		p_o_mm[1] = 2+1;
		break;
	case OP_R_P_PM0:
		p_o_mm[0] = p_i;
		p_o_mm[1] = p_i+1;
		break;
	case OP_R_P_PM1:
		p_o_mm[0] = GSL_MAX(p_i-1,p_ref[0]);
		p_o_mm[1] = GSL_MIN(p_i+1,p_ref[1])+1;
		break;
	case OP_R_P_1PPM1:
		p_o_mm[0] = GSL_MAX(p_ref[0]-1,p_ref[0]);
		p_o_mm[1] = GSL_MIN(p_ref[1]+1,p_ref[1])+1;
		break;
	case OP_R_P_1P: // fallthrough
	case OP_R_P_ALL:
		p_o_mm[0] = p_ref[0];
		p_o_mm[1] = p_ref[1]+1;
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}
	assert(p_o_mm[0] >= 0);
	assert(p_o_mm[1] >= p_o_mm[0]);
}

static ptrdiff_t get_ind_nc
	(const struct Operator_Info* op_info, const int* op_values, const struct const_Multiarray_Vector_i* nc)
{
	const int order_op = (int)op_info->extents_op->ext_0;
	const struct const_Vector_i* indices_op = constructor_indices_Vector_i(order_op,op_values,NULL); // destructed

	ptrdiff_t ind_nc = compute_index_sub_container_pi(nc->order,1,nc->extents,indices_op->data);
	destructor_const_Vector_i(indices_op);

	return ind_nc;
}

static ptrdiff_t get_ind_w
	(const struct Operator_Info* op_info, const int* op_values, const struct const_Multiarray_Vector_d* w)
{
	const int order_op = (int)op_info->extents_op->ext_0;
	const struct const_Vector_i* indices_op = constructor_indices_Vector_i(order_op,op_values,NULL); // destructed

	ptrdiff_t ind_w = compute_index_sub_container_pi(w->order,0,w->extents,indices_op->data);
	destructor_const_Vector_i(indices_op);

	return ind_w;
}

// Level 2 ********************************************************************************************************** //

/** \brief Get the type of coarse to fine projection which is resulting in the usage of the L2 operator.
 *  \return 'h'-projection, 'p'-projection. */
static char get_op_L2_type
	(const int*const op_values ///< Values for the operator indices.
	);

/// \brief Update \ref Operator_Info::op_type such that the output is of type 'c'oefficient.
static void update_op_type_Xc
	(const int*const op_type ///< \ref Operator_Info::op_type.
	);

static const struct const_Multiarray_Matrix_d* constructor_op_Xc_ce_io
	(const ptrdiff_t ind_values, const struct Operator_Info* op_info, const struct Simulation* sim)
{
	const int n_op = 1;

	int* op_values = (int*) get_row_const_Matrix_i(ind_values,op_info->values_op);
	const struct Op_IO* op_io = op_info->op_io;

	const char type_L2 = get_op_L2_type(op_values);

	// Should only have a fine to coarse projection from different computational elements of type 'p'.
	assert(type_L2 == 'p');

	const int op_type = op_info->op_type;
	update_op_type_Xc(&op_info->op_type);

	const int p_i = op_values[OP_IND_P+OP_IND_I],
	          p_o = op_values[OP_IND_P+OP_IND_O];
	op_values[OP_IND_P+OP_IND_O] = p_i;

	const char kind_i = op_io[OP_IND_I].kind,
	           kind_o = op_io[OP_IND_O].kind;
	const_cast_c(&op_io[OP_IND_O].kind,kind_i);

	const char sc_i = op_io[OP_IND_I].sc,
	           sc_o = op_io[OP_IND_O].sc;
	const_cast_c(&op_io[OP_IND_O].sc,sc_i);

	const struct const_Multiarray_Matrix_d* op_Xc_ce_io = constructor_op_MMd(true,n_op); // returned
	set_operator_std_interp(ind_values,op_info,sim,op_Xc_ce_io);

	const_cast_i(&op_info->op_type,op_type);
	op_values[OP_IND_P+OP_IND_O] = p_o;
	const_cast_c(&op_io[OP_IND_O].kind,kind_o);
	const_cast_c(&op_io[OP_IND_O].sc,sc_o);

	return op_Xc_ce_io;
}

static const struct const_Matrix_d* constructor_mass
	(const struct const_Matrix_d* cv0_l, const struct const_Matrix_d* cv0_r, const struct const_Vector_d* w,
	 const double volume_ratio)
{
	const struct const_Matrix_d* mass_r = constructor_mm_diag_const_Matrix_d(1.0,cv0_r,w,'L',false); // destructed

	const struct const_Matrix_d* mass =
		constructor_mm_const_Matrix_d('T','N',volume_ratio,cv0_l,mass_r,'R'); // returned
	destructor_const_Matrix_d(mass_r);

	return mass;
}

static char get_op_type_out (const int op_type)
{
	switch (op_type) {
	case OP_T_CV: case OP_T_VV:
		return 'v';
		break;
	case OP_T_CC: case OP_T_VC:
		return 'c';
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",op_type);
		break;
	}
}

static double compute_volume_ratio (const int e_type, const int ind_h)
{
	if (ind_h == 0)
		return 1.0;

	switch (e_type) {
	case LINE:
		switch (ind_h) {
		case H_LINE2_V0: case H_LINE2_V1:
			return 0.5;
			break;
		default:
			EXIT_ERROR("Unsupported: %d\n",ind_h);
			break;
		}
		break;
	case TRI:
		switch (ind_h) {
		case H_TRI4_V0: case H_TRI4_V1: case H_TRI4_V2: case H_TRI4_V3:
			return 0.25;
			break;
		default:
			EXIT_ERROR("Unsupported: %d\n",ind_h);
			break;
		}
		break;
	case TET:
		EXIT_ADD_SUPPORT;
		switch (ind_h) {
//		case H_TET8_V0: case H_TET8_V1: case H_TET8_V2: case H_TET8_V3:
//		case H_TET8_V4: case H_TET8_V5: case H_TET8_V6: case H_TET8_V7:
//			return 0.125;
//			break;
		default:
			EXIT_ERROR("Unsupported: %d\n",ind_h);
			break;
		}
		break;
	case PYR:
		switch (ind_h) {
		case H_PYR10_V0: case H_PYR10_V1: case H_PYR10_V2: case H_PYR10_V3:
		case H_PYR10_V4: case H_PYR10_V5: case H_PYR10_V6: case H_PYR10_V7:
		case H_PYR10_V8: case H_PYR10_V9:
			return 0.125;
			break;
		default:
			EXIT_ERROR("Unsupported: %d\n",ind_h);
			break;
		}
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",e_type);
		break;
	}
}

// Level 3 ********************************************************************************************************** //

static char get_op_L2_type (const int*const op_values)
{
	if (op_values[OP_IND_H+OP_IND_I] > op_values[OP_IND_H+OP_IND_O])
		return 'h';
	else if (op_values[OP_IND_P+OP_IND_I] > op_values[OP_IND_P+OP_IND_O])
		return 'p';
	else
		EXIT_ERROR("Should not be here.\n");
}

static void update_op_type_Xc (const int*const op_type)
{
	int op_type_mod = -1;
	switch (*op_type) {
		case OP_T_CV: case OP_T_CC:
			op_type_mod = OP_T_CC;
			break;
		case OP_T_VV: case OP_T_VC:
			op_type_mod = OP_T_VC;
			break;
		default:
			EXIT_ERROR("Unsupported: %d\n",*op_type);
			break;
	}
	const_cast_i(op_type,op_type_mod);
}
