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
/// \todo Delete unused functions.

#include "element_operators.h"

#include <assert.h>
#include <string.h>

#include "macros.h"
#include "definitions_core.h"
#include "definitions_cubature.h"
#include "definitions_element_operators.h"
#include "definitions_elements.h"

#include "multiarray.h"
#include "vector.h"

#include "simulation.h"
#include "element.h"
#include "const_cast.h"
#include "bases.h"
#include "cubature.h"

// Static function declarations ************************************************************************************* //

///\{ \name Invalid operator index
#define OP_INVALID_IND -999
///\}

/** \brief Constructor for the \ref Operator_Info\* having the given inputs.
 *  \return Standard. */
struct Operator_Info* constructor_Operator_Info
	(const char*const name_type,         ///< Defined for \ref constructor_operators.
	 const char*const name_in,           ///< Defined for \ref constructor_operators.
	 const char*const name_out,          ///< Defined for \ref constructor_operators.
	 const char*const name_range,        ///< Defined for \ref constructor_operators.
	 const int p_ref[2],                 ///< Defined for \ref constructor_operators.
	 const struct const_Element* element ///< Defined for \ref constructor_operators.
	);

/// \brief Destructor for a \ref Operator_Info\*.
void destructor_Operator_Info
	(struct Operator_Info* op_ranges ///< Standard.
	);

/** \brief Set the current standard operator(s).
 *  \note Standard operators include all permutations of coef/value to coef/value with possible differentiation. */
static void set_operator_std
	(const int ind_values,                       /**< The index of the first row of \ref Operator_Info::op_values
	                                              *   to use. */
	 const struct const_Multiarray_Matrix_d* op, ///< The multiarray of operators.
	 struct Operator_Info* op_info               ///< \ref Operator_Info.
	);


/** \brief Constructor for the \ref Multiarray_Cubature::data.
 *  \return Standard. */
/*const struct const_Cubature** constructor_cub_data_array
	(const struct Simulation* sim,        ///< \ref Simulation.
	 const struct const_Element* element, ///< \ref const_Element.
	 const struct Operator_Info* op_info  ///< \ref Operator_Info.
	);

/// \brief Destructor for a \ref Multiarray_Cubature\* container.
static void destructor_Multiarray_Cubature
	(struct Multiarray_Cubature* a ///< Standard.
	);*/

// Interface functions ********************************************************************************************** //

/*
const struct const_Multiarray_Cubature* constructor_const_Multiarray_Cubature
	(const struct Simulation* sim, const struct const_Element* element, const struct Operator_Info* op_info)
{
	// Compute \ref Cubature data.
	const struct const_Cubature** cub_data = constructor_cub_data_array(sim,element,op_info); // keep

	// Construct the Multiarray.
	struct Multiarray_Cubature* cub_MA = malloc(sizeof *cub_MA); // returned

	cub_MA->order     = order;
	cub_MA->extents   = extents;
	cub_MA->owns_data = true;
	cub_MA->data      = cub_data;

	return (const struct const_Multiarray_Cubature*) cub_MA;
}

void destructor_const_Multiarray_Cubature (const struct const_Multiarray_Cubature*const a)
{
	destructor_Multiarray_Cubature((struct Multiarray_Cubature*)a);
}*/

const struct const_Multiarray_Matrix_d* constructor_operators
	(const char*const name_type, const char*const name_in, const char*const name_out, const char*const name_range,
	 const int p_ref[2], const struct const_Element* element, const struct Simulation* sim)
{
	struct Operator_Info* op_info =
		constructor_Operator_Info(name_type,name_in,name_out,name_range,p_ref,element); // destructed

	const ptrdiff_t order = op_info->extents_op->ext_0;
	const struct const_Multiarray_Matrix_d* op =
		constructor_empty_const_Multiarray_Matrix_d_V(false,op_info->extents_op); // returned

	ptrdiff_t* extents = op->extents;

	const ptrdiff_t row_max = op_info->values_op->ext_0;
	for (ptrdiff_t row = 0; row < row_max; ++row) {
		switch (op_info->op_type) {
		case OP_T_CV:
		case OP_T_CC:
		case OP_T_VV:
		case OP_T_VC:
// Currently only Volume -> Volume, Face, Edge operators.
			set_operator_std(row,op,op_info);
			break;
		default:
			EXIT_UNSUPPORTED;
			break;
		}
	}

	destructor_Operator_Info(op_info);

	return (const struct const_Multiarray_Matrix_d*) op;
}



/** \brief Check if information can be lost while performing the operation.
 *  \return `true` if: `p_i > p_o` or `h_i > h_o`; `false` otherwise. */
bool check_op_info_loss
	(const int*const op_values)
{
	if ((op_values[OP_IND_H+1] > op_values[OP_IND_H]) ||
	    (op_values[OP_IND_P+1] > op_values[OP_IND_P]))
		return true;
	return false;
}




// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Convert the `char*` input to the appropriate definition of OP_R_D_*.
 *  \return See brief. */
int convert_to_range_d
	(char* name_type ///< Defined for \ref constructor_Operator_Info.
	);

/** \brief Convert the `char*` inputs to the appropriate definition of OP_R_CE_*.
 *  \return See brief. */
int convert_to_range_ce
	(const char ce_i, ///< `name_in[0]`  which is defined for \ref constructor_Operator_Info.
	 const char ce_o  ///< `name_out[0]` which is defined for \ref constructor_Operator_Info.
	);

/** \brief Convert the `char*` input to the appropriate definition of OP_R_*.
 *  \return See brief. */
int convert_to_range
	(const char type_range       ///< The type of range parameter. Options: 'h', 'p'.
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

/** \brief Constructor for a \ref Vector_i\* of indices for the current operator.
 *  \return See brief. */
static struct Vector_i* constructor_indices_Vector_i
	(const int order_expected,     ///< The expected order.
	 int* op_values,               ///< The operator values.
	 const bool*const indices_skip ///< Indices to skip (if not NULL).
	);

struct Operator_Info* constructor_Operator_Info
	(const char*const name_type, const char*const name_in, const char*const name_out, const char*const name_range,
	 const int p_ref[2], const struct const_Element* element)
{
	const int op_type = convert_to_type(name_type);

	const int ranges[] =
		{ convert_to_range_d(name_type),
	          convert_to_range_ce(name_in[0],name_in[1]),
	          convert_to_range('h',name_range),
	          convert_to_range('p',name_range), };

	struct Operator_Info* op_info = malloc(sizeof *op_info); // returned

	op_info->element = element;

	const_cast_i(&op_info->op_type,op_type);

	const_cast_c(&op_info->op_io[OP_IND_IN].ce,  name_in[0]);
	const_cast_c(&op_info->op_io[OP_IND_IN].kind,name_in[1]);
	const_cast_c(&op_info->op_io[OP_IND_IN].sc,  name_in[2]);

	const_cast_c(&op_info->op_io[OP_IND_OUT].ce,  name_out[0]);
	const_cast_c(&op_info->op_io[OP_IND_OUT].kind,name_out[1]);
	const_cast_c(&op_info->op_io[OP_IND_OUT].sc,  name_out[2]);

	const_cast_i(&op_info->range_d,ranges[0]);
	const_cast_i(&op_info->range_f,ranges[1]);
	const_cast_i(&op_info->range_h,ranges[2]);
	const_cast_i(&op_info->range_p,ranges[3]);

	const_cast_i1(op_info->p_ref,p_ref,2);

	set_up_extents(op_info);
	set_up_values_op(op_info);

	return op_info;
}

void destructor_Operator_Info (struct Operator_Info* op_info)
{
	destructor_Vector_i(op_info->extents_cub);
	destructor_Vector_i(op_info->extents_op);
	destructor_Matrix_i(op_info->values_op);
	free(op_info);
}

static void set_operator_std
	(const int ind_values, const struct const_Multiarray_Matrix_d* op, const struct Operator_Info* op_info)
{
	const bool info_loss = check_op_info_loss(op_values);

	if (!info_loss) {
		const int*const h_ptr = &op_values[OP_IND_H],
		         *const p_ptr = &op_values[OP_IND_P];

		// Construct T_i
		const int p_i         = p_ptr[OP_IND_I];
		const int node_type_i = // Maybe move node_type computation into cubature function. Pass op_info instead.
		const struct const_Cubature*const cub_i =
			constructor_const_Cubature_h(d,p_x[1],node_type_i,s_type,h_ptr[1]); // destructed
	}



	const struct const_Element* element = op_info->element;

	const int d      = element->d,
	          s_type = element->s_type;

/*	const int* op_values = get_row_const_Matrix_i(ind_values,op_info->values_op);

	const int order_op  = op_info->extents_op->ext_0,
	          order_cub = op_info->extents_cub->ext_0;

	struct Vector_i* indices_op = constructor_indices_Vector_i(order_op,op_values,NULL); // destructed

	const bool skip_cub[] = {true,false,false,true,false,true};
	struct Vector_i* indices_cub = constructor_indices_Vector_i(order_cub,op_values,skip_cub); // destructed
*/

// separate function here
//const ptrdiff_t ind_cub_array = compute_index_sub_container(order_cub,1,cub_array->extents,indices_cub);
//const struct const_Cubature*const cub = cub_array->data[ind_cub_array];

//	destructor_Vector_i(indices_op);
//	destructor_Vector_i(indices_cub);



//grad_basis_fptr constructor_grad_basis = get_grad_basis_by_super_type(element->s_type,"ortho");
	if (!info_loss) {
		const int*const p_x[2] = { compute_p_cub(p_ptr[0],op_info->cub_type,s_type,sim),
		                           compute_p_cub(p_ptr[1],op_info->cub_type,s_type,sim), };
// NOTE: nodes out and basis out must coincide when this operator is available. I.e. you cannot project a non-basis set
// of nodes to coefficients as you would never have the cc, vc operator.
// Thus, only two sets of cubature nodes are required.

// Get cubature nodes possibly restricted to subset of element.
		const int basis_node_type = get_node_type(basis comp elem (v,f,e), basis type (g,m,s,...),s_type);
		const struct const_Cubature*const cub_i =
			constructor_const_Cubature_h(d,p_x[1],basis_node_type,s_type,h_ptr[1]); // destructed
		const int cub_node_type = get_node_type(cub comp elem (v,f,e), cub type (g,m,s,...),s_type);
		const struct const_Cubature*const cub_o =
			constructor_const_Cubature_h(d,p_x[0],cub_node_type,s_type,h_ptr[0]); // destructed

// Construct cv
		basis_fptr constructor_basis_ref = get_basis_by_super_type(s_type,"ortho");
		basis_fptr constructor_basis     = get_basis_by_super_type(s_type,"");
		const struct const_Matrix_d* phi_ref_cv_ii = constructor_basis_ref(p_x[1],d,cub_i->rst); // tbd
// external function to compute T below.
		const struct const_Matrix_d* phi_cv_ii = NULL;
		if (ortho)
			phi_cv_ii = constructor_copy_const_Matrix_d(phi_ref_cv_ii); // tbd
		else if (lagrange)
			phi_cv_ii = constructor_identity_const_Matrix_d(phi_ref_cv_ii->ext_0); // tbd
		else if (bezier)
			phi_cv_ii = constructor_basis_bezier(p_x[1],d,cub_i->rst); // tbd
		else
			EXIT_UNSUPPORTED;

		const struct const_Matrix_d* T_i = constructor_sgesv_const_Matrix_d(phi_ref_cv_ii,phi_cv_ii); // tbd

if (D_OP_R_0)
		const struct const_Matrix_d* phi_ref_cv_io = constructor_basis_ref(p_x[1],d,cub_o->rst); // tbd
		const struct const_Matrix_d* phi_cv_io     = constructor_mm_const_Matrix_d(phi_ref_cv_io,T_i); // tbd

		const struct const_Matrix_d* output = NULL;
		if (op_type == OP_T_CV)
			op_out = phi_cv_io;

		if (op_type == OP_T_CC) {
			const struct const_Matrix_d* phi_ref_oo = constructor_basis_ref(p_x[0],d,cub_o->rst);
		}
else
	// compute

// Update for other operators

// index of first operator: const ptrdiff_t ind_op_MA = compute_index_sub_container(OP_ORDER_MAX,0,extents,indices_op);
// if: op_values[0] == OP_INVALID_IND, set from Matrix_d.
// else: loop over dim set from Multiarray_Matrix_d.

	} else {
		EXIT_ADD_SUPPORT; // L2 projection operators.
	}
}


// Level 1 ********************************************************************************************************** //

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

static void set_up_extents (struct Operator_Info* op_info)
{
/// \todo Check whether extents_cub is actually needed. Delete if not.
	struct Vector_i* extents_cub = constructor_default_Vector_i(); // keep
	struct Vector_i* extents_op  = constructor_default_Vector_i(); // keep

	switch (op_info->range_d) {
		case OP_R_D_0:
			/* Do nothing */
			break;
		case OP_R_D_ALL:
			push_back_Vector_i(extents_op,element->d,false,false);
			break;
		default:
			EXIT_UNSUPPORTED;
			break;
	}

	switch (op_info->range_ce) {
//			push_back_Vector_i(extents_cub,element->n_f,false,false);
		case OP_R_CE_VV:
		case OP_R_CE_FF:
		case OP_R_CE_EE:
			/* Do nothing */
			break;
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
			// Note: Potential special cases for wedge/pyr elements with faces having varying number of
			//       edges.
			EXIT_ADD_SUPPORT;
			break;
		default:
			EXIT_UNSUPPORTED;
			break;
	}

	switch (op_info->range_h) { // h_o, h_i
		// There are currently no h-adaptation operators for other `range_ce` values.
		assert((op_info->range_ce == OP_R_CE_VV) ||
		       (op_info->range_ce == OP_R_CE_VF));

		int n_ref_max = -1;
		if (op_info->range_ce == OP_R_CE_VV)
			n_ref_max = op_info->element->n_ref_max;
		else if (op_info->range_ce == OP_R_CE_VF)
			n_ref_max = op_info->element->n_ref_f_max;

		case OP_R_H_1:
			push_back_Vector_i(extents_cub,1,false,false);
			push_back_Vector_i(extents_op,1,false,false);
			push_back_Vector_i(extents_op,1,false,false);
			break;
		case OP_R_H_CF:
			push_back_Vector_i(extents_cub,n_ref_max,false,false);
			push_back_Vector_i(extents_op,n_ref_max,false,false);
			push_back_Vector_i(extents_op,1,false,false);
			break;
		case OP_R_H_FC:
			push_back_Vector_i(extents_cub,n_ref_max,false,false);
			push_back_Vector_i(extents_op,1,false,false);
			push_back_Vector_i(extents_op,n_ref_max,false,false);
			break;
		default:
			EXIT_UNSUPPORTED;
			break;
	}

	switch (op_info->range_p) { // p_o, p_i
		case OP_R_P_1:
			push_back_Vector_i(extents_cub,1,false,false);
			push_back_Vector_i(extents_op,1,false,false);
			push_back_Vector_i(extents_op,1,false,false);
			break;
		case OP_R_P_PM0: /* fallthrough */
		case OP_R_P_PM1: /* fallthrough */
		case OP_R_P_ALL:
			push_back_Vector_i(extents_cub,op_info->p_ref[1]+1,false,false);
			push_back_Vector_i(extents_op,op_info->p_ref[1]+1,false,false);
			push_back_Vector_i(extents_op,op_info->p_ref[1]+1,false,false);
			break;
		default:
			EXIT_UNSUPPORTED;
			break;
	}

	op_info->extents_cub = (const struct const_Vector_i*) extents_cub;
	op_info->extents_op  = (const struct const_Vector_i*) extents_op;
}

static void set_up_values_op (struct Operator_Info* op_info)
{
	const ptrdiff_t size = prod_Vector_i(op_info->extents);

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
		for (int ce_i = ce_i_mm[0]; ce_i < ce_i_mm[1]; ++ce_i) {
		for (int ce_o = ce_o_mm[0]; ce_o < ce_o_mm[1]; ++ce_o) {
		for (int d = d_mm[0]; d < d_mm[1]; ++d) {
			set_row_Matrix_i(row,values,(int[]){d,ce_o,ce_i,h_o,h_i,p_o,p_i});
			++row;
		}}}}}}
	}
	values->ext_0 = row;

	op_info->values_op = (const struct const_Matrix_i*) values;
}

// Level 2 ********************************************************************************************************** //

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
	case 'c': // ce
		switch (op_info->range_ce) {
		case OP_R_CE_VV:
		case OP_R_CE_FF:
		case OP_R_CE_EE:
			x_mm[0] = OP_INVALID_IND;
			x_mm[1] = OP_INVALID_IND+1;
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
	case 'h':
		// There are currently no h-adaptation operators for other `range_ce` values.
		assert((op_info->range_ce == OP_R_CE_VV) ||
		       (op_info->range_ce == OP_R_CE_VF));

		int n_ref_max = -1;
		if (op_info->range_ce == OP_R_CE_VV)
			n_ref_max = op_info->element->n_ref_max;
		else if (op_info->range_ce == OP_R_CE_VF)
			n_ref_max = op_info->element->n_ref_f_max;

		x_mm[0] = 0;
		switch (op_info->range_h) {
		case OP_R_H_1:
			x_mm[1] = 1;
			break;
		case OP_R_H_CF:
			if (var_io == 'i')
				x_mm[1] = 1;
			else if (var_io == 'o')
				x_mm[1] = n_ref_max;
			break;
		case OP_R_H_FC:
			if (var_io == 'i')
				x_mm[1] = n_ref_max;
			else if (var_io == 'o')
				x_mm[1] = 1;
			break;
		default:
			EXIT_UNSUPPORTED;
			break;
		}
	case 'p':
		// In general, p_o depends on p_i and requires more information than that provided here.
		assert(var_io == 'i');

		switch (op_info->range_p) {
		case OP_R_P_1:
			x_mm[0] = 1;
			x_mm[1] = 2;
			break;
		case OP_R_P_PM0:
		case OP_R_P_PM1: /* fallthrough */
		case OP_R_P_ALL:
			x_mm[0] = op_info->p_ref[0];
			x_mm[1] = op_info->p_ref[1];
			break;
		default:
			EXIT_UNSUPPORTED;
			break;
		}
	default:
		EXIT_UNSUPPORTED;
		break;
	}
}

static void compute_range_p_o (int p_o_mm[2], const struct Operator_Info* op_info, const int p_i)
{
	switch (op_info->range_p) {
	case OP_R_P_1:
		p_o_mm[0] = 1;
		p_o_mm[1] = 2;
		break;
	case OP_R_P_PM0:
		p_o_mm[0] = p_i;
		p_o_mm[1] = p_i+1;
		break;
	case OP_R_P_PM1:
		p_o_mm[0] = max(p_i-1,op_info->p_ref[0]);
		p_o_mm[1] = min(p_i+1,op_info->p_ref[1])+1;
		break;
	case OP_R_P_ALL:
		x_mm[0] = op_info->p_ref[0];
		x_mm[1] = op_info->p_ref[1]+1;
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}
}









/** \brief Compute the node type associated with the input cubature type.
 *  \return See brief. */
static int compute_node_type
	(const struct Operator_Info* op_info, ///< \ref Operator_Info.
	 const int s_type,                    ///< \ref Element::s_type.
	 const struct Simulation *sim         ///< \ref Simulation.
	);

/// \brief Compute the minimum and maximum orders for the cubature nodes of the given type.
static void compute_p_range
	(int*const p_min,                    ///< The minimum order in the range.
	 int*const p_max,                    ///< The maximum order in the range.
	 const struct Operator_Info* op_info ///< \ref Operator_Info.
	);

/** \brief Compute the cubature node order associated with the current reference order for the give cub_type.
 *  \return See brief. */
static int compute_p_cub
	(const int p_ref,             ///< The input refenrece order.
	 const int cub_type,          ///< \ref Operator_Info::cub_type.
	 const int s_type,            ///< \ref Element::s_type.
	 const struct Simulation *sim ///< \ref Simulation.
	);

/** \brief Compute the index of computational element associated with the cubature type.
 *  \return See brief. */
static int compute_cub_ce
	(const int cub_type ///< The cubature type.
	);

static void destructor_Multiarray_Cubature (struct Multiarray_Cubature* a)
{
	assert(a != NULL);

	if (a->owns_data) {
		const ptrdiff_t size = compute_size(a->order,a->extents);
		for (ptrdiff_t i = 0; i < size; i++)
			destructor_const_Cubature(a->data[i]);
		free(a->data);
	}
	free(a->extents);
	free(a);
}

const struct const_Cubature** constructor_cub_data_array
	(const struct Simulation* sim, const struct const_Element* element, const struct Operator_Info* op_info)
{
	const int order = op_info->order_cub;
	const ptrdiff_t*const extents = op_info->extents_cub;

	const ptrdiff_t size = compute_size(order,extents);
	const struct const_Cubature** cub_data = malloc(size * sizeof *cub_data); // returned

	const int d = element->d;
	cubature_fptr constructor_Cubature = get_cubature_by_super_type(element->s_type);

	const int comp_elem = compute_cub_ce(op_info->cub_type);
	if (comp_elem == CUB_CE_V) {
		assert(order == 2);

		const int node_type = compute_node_type(op_info,element->s_type,sim);

		int p_min = -1,
		    p_max = -1;
		compute_p_range(&p_min,&p_max,op_info);
		for (ptrdiff_t ind = 0, p = p_min; p <= p_max; ++p) {
			compute_p_cub(p,op_info->cub_type,element->s_type,sim);
			for (ptrdiff_t h = 0; h < extents[0]; ++h) {
				if (h > 0)
					EXIT_ADD_SUPPORT; // sub-elements
				cub_data[ind++] = constructor_Cubature(d,p,node_type); // keep
			}
		}
	} else if (comp_elem == CUB_CE_F) {
		assert(order == 3);
		EXIT_ADD_SUPPORT;
	} else {
		EXIT_UNSUPPORTED;
	}

	return cub_data;
}

static struct Vector_i* constructor_indices_Vector_i
	(const int order_expected, int* op_values, const bool*const indices_skip)
{
	struct Vector_i* indices = constructor_empty_Vector_i(order_expected); // returned
	int* data = indices->data;

	int ind = 0;
	for (int i = 0; i < OP_ORDER_MAX; ++i) {
		if (!(indices_skip && indices_skip[i]) && op_values[i] != OP_INVALID_IND)
			data[ind++] = op_values[i];
	}
	assert(ind == order_expected);

	return indices;
}


int convert_to_range_d (char* name_type)
{
	long order_d = -1;
	while (*name_type) {
		if (isdigit(*name_type)) {
			order_d = strtol(name_type,&name_type,10);
			break;
		} else {
			++name_type;
		}
	}

	switch (order_d) {
		case 0:  return OP_R_D_0;   break;
		case 1:  return OP_R_D_ALL; break;
		default: EXIT_UNSUPPORTED;  break;
	}
	EXIT_ERROR("Did not find the operator range.");
	return -1;
}

int convert_to_range_ce (const char ce_i, const char ce_o)
{
	switch (ce_i) {
	case 'v':
		switch (ce_o) {
			case 'v': return OP_R_CE_VV;   break;
			case 'f': return OP_R_CE_VF;   break;
			case 'e': return OP_R_CE_VE;   break;
			default:  EXIT_UNSUPPORTED; break;
		}
		break;
	case 'f':
		switch (ce_o) {
			case 'v': return OP_R_CE_FV;   break;
			case 'f': return OP_R_CE_FF;   break;
			case 'e': return OP_R_CE_FE;   break;
			default:  EXIT_UNSUPPORTED; break;
		}
		break;
	case 'e':
		switch (ce_o) {
			case 'v': return OP_R_CE_EV;   break;
			case 'f': return OP_R_CE_EF;   break;
			case 'e': return OP_R_CE_EE;   break;
			default:  EXIT_UNSUPPORTED; break;
		}
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}
	EXIT_ERROR("Did not find the operator range.");
	return -1;
}

int convert_to_range (const char type_range, const char*const name_range)
{
	switch (type_range) {
	case 'h':
		if (strstr(name_range,"H_1"))
			return OP_R_H_1;
		else if (strstr(name_range,"H_CF"))
			return OP_R_H_CF;
		else if (strstr(name_range,"H_FC"))
			return OP_R_H_FC;
		else
			EXIT_UNSUPPORTED;
		break;
	case 'p':
		if (strstr(name_range,"P_1"))
			return OP_R_P_1;
		else if (strstr(name_range,"P_PM0"))
			return OP_R_P_PM0;
		else if (strstr(name_range,"P_PM1"))
			return OP_R_P_PM1;
		else if (strstr(name_range,"P_ALL"))
			return OP_R_P_ALL;
		else
			EXIT_UNSUPPORTED;
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}
	EXIT_ERROR("Did not find the operator range.");
	return -1;
}

/** \brief Compute the index of straight/curved indicator associated with the cubature type.
 *  \return See brief. */
static int compute_cub_sc
	(const int cub_type ///< The cubature type.
	);

/** \brief Compute the index of entity associated with the cubature type.
 *  \return See brief. */
static int compute_cub_ent
	(const int cub_type ///< The cubature type.
	);

static int compute_node_type (const struct Operator_Info* op_info, const int s_type, const struct Simulation *sim)
{
	const int cub_ent = compute_cub_ent(op_info->cub_type);
	switch (cub_ent) {
	case CUB_ENT_G:
		if (s_type == ST_TP || s_type == ST_PYR)
			return CUB_GLL;
		else if (s_type == ST_SI)
			return CUB_AO;
		else
			EXIT_UNSUPPORTED;
		break;
UNUSED(sim);
	default:
		EXIT_UNSUPPORTED;
		break;
	}
}

static void compute_p_range
	(int*const p_min, int*const p_max, const struct Operator_Info* op_info)
{
	switch (op_info->range_p) {
		case OP_R_P_1:
			*p_min = 1;
			*p_max = 1;
			break;
		case OP_R_P_PM0: /* fallthrough */
		case OP_R_P_PM1: /* fallthrough */
		case OP_R_P_ALL: /* fallthrough */
			*p_min = op_info->p_ref[0];
			*p_max = op_info->p_ref[1];
			break;
		default:
			EXIT_UNSUPPORTED;
			break;
	}
}

static int compute_p_cub (const int p_ref, const int cub_type, const int s_type, const struct Simulation *sim)
{
	const int cub_sc  = compute_cub_sc(cub_type),
	          cub_ce  = compute_cub_ce(cub_type),
	          cub_ent = compute_cub_ent(cub_type);
	switch (cub_ent) {
	case CUB_ENT_S:
		return p_ref;
		break;
	case CUB_ENT_G: {
		if (cub_sc == CUB_SC_S)
			return 1;

		if (strcmp(sim->geom_rep,"isoparametric") == 0)
			return p_ref;
		else if (strcmp(sim->geom_rep,"superparametric") == 0)
			return p_ref+1;
		else if (strstr(sim->geom_rep,"fixed"))
			EXIT_ADD_SUPPORT;
		else
			EXIT_UNSUPPORTED;
		break;
	} case CUB_ENT_M: {
//		const int p_g = compute_p_cub (p_ref,cub_sc+cub_ce+CUB_ENG_G,s_type,sim);
UNUSED(cub_ce);
UNUSED(s_type);

		break;
	} default:
		EXIT_UNSUPPORTED;
		break;
	}
	EXIT_UNSUPPORTED;
	return -1;
}

static int compute_cub_ce (const int cub_type)
{
	const int cub_sc       = compute_cub_sc(cub_type),
	          cub_ce_inter = cub_type - cub_sc;

	// note the integer division.
	return (cub_ce_inter / CUB_CE_MULT)*CUB_CE_MULT;
}



static int compute_cub_sc (const int cub_type)
{
	// note the integer division.
	return (cub_type / CUB_SC_MULT)*CUB_SC_MULT;
}

static int compute_cub_ent (const int cub_type)
{
	const int cub_sc = compute_cub_sc(cub_type),
	          cub_ce = compute_cub_ce(cub_type);
	return cub_type - cub_sc - cub_ce;
}
