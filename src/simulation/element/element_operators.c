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

///\{ \name Operator related parameters
#define OP_ORDER_MAX 6 // d(1) + f(1) + h(2) + p(2)
#define OP_IND_H     2
#define OP_IND_P     4
///\}

///\{ \name Invalid operator index
#define OP_INVALID_IND -999
///\}

/** \brief Constructor for the \ref Multiarray_Cubature::data.
 *  \return Standard. */
const struct const_Cubature** constructor_cub_data_array
	(const struct Simulation* sim,        ///< \ref Simulation.
	 const struct const_Element* element, ///< \ref const_Element.
	 const struct Operator_Info* op_info  ///< \ref Operator_Info.
	);

/// \brief Destructor for a \ref Multiarray_Cubature\* container.
static void destructor_Multiarray_Cubature
	(struct Multiarray_Cubature* a ///< Standard.
	);

/// \brief Set up \ref Operator_Info order_* and extents_* members.
static void set_up_order_and_extents
	(struct Operator_Info* op_info ///< \ref Operator_Info.
	);

/// \brief Set up \ref Operator_Info::values_op.
static void set_up_basis_values
	(struct Operator_Info* op_info ///< \ref Operator_Info.
	);

// Interface functions ********************************************************************************************** //

struct Operator_Info* constructor_Operator_Info
	(const int range_d, const int range_f, const int range_p, const int range_h, const int cub_type,
	 const int p_ref[2], const struct const_Element* element)
{
	struct Operator_Info* op_info = malloc(sizeof *op_info); // returned

	op_info->element = element;

	const_cast_i(&op_info->range_d,range_d);
	const_cast_i(&op_info->range_f,range_f);
	const_cast_i(&op_info->range_h,range_h);
	const_cast_i(&op_info->range_p,range_p);

	const_cast_i(&op_info->cub_type,cub_type);
	const_cast_i1(op_info->p_ref,p_ref,2);

	set_up_order_and_extents(op_info);
	set_up_basis_values(op_info);

	return op_info;
}

void destructor_Operator_Info (struct Operator_Info* op_info)
{
	free(op_info->extents_cub);
	free(op_info->extents_op);
	destructor_Matrix_i(op_info->values_op);
	free(op_info);
}

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
}

struct Multiarray_Matrix_d* constructor_operators_Multiarray_Matrix_d_V
	()
const int op_type,    ///< Operator type (CV,CC,VV,VC,DG_Weak_VV,...)
const int diff_order, ///< Differentiation order.
const struct Simulation* sim,
const struct const_Element* element,
const struct Operator_Info* op_info,
//const struct const_Multiarray_Cubature* cub_array
{
	const int order = op_info->order_op;
	const ptrdiff_t*const extents = op_info->extents_op;

	struct Multiarray_Matrix_d* op = constructor_empty_Multiarray_Matrix_d_V(false,order,extents); // returned

	const ptrdiff_t size = compute_size(order,extents);

	basis_fptr      constructor_basis      = get_basis_by_super_type     (element->s_type,"ortho");
	grad_basis_fptr constructor_grad_basis = get_grad_basis_by_super_type(element->s_type,"ortho");

	const ptrdiff_t row_max = op_info->values_op->ext_0;
	for (ptrdiff_t row = 0; row < row_max; ++row) {
		int* op_values = get_row_Matrix_i(row,op_info->values_op);

// separate function here and below
		int ind = 0;
		ptrdiff_t indices_op[order];
		for (int i = 0; i < OP_ORDER_MAX; ++i) {
			if (op_values[i] != OP_INVALID_IND)
				indices_op[ind++] = op_values[i];
		}
		assert(ind == order);

		const int order_cub = op_info->order_cub;

		ind = 0;
		const int indices_skip[OP_ORDER_MAX] = {1,0,0,1,0}; // cubature only uses, f, h, p_out
		ptrdiff_t indices_cub[order_cub];
		for (int i = 0; i < OP_ORDER_MAX; ++i) {
			if (!indices_skip[i] && op_values[i] != OP_INVALID_IND)
				indices_cub[ind++] = op_values[i];
		}
		assert(ind == order_cub);
//		const ptrdiff_t row = compute_index_sub_container(OP_ORDER_MAX,0,extents,(ptrdiff_t[]){d,f,h,p_o,p_i});

// separate function here
		const ptrdiff_t ind_basis = indices_op[order-2];
//		const ptrdiff_t ind_cub_array = compute_index_sub_container(order_cub,1,cub_array->extents,indices_cub);
//		const struct const_Cubature*const cub = cub_array->data[ind_cub_array];
		if (op_values[0] == OP_INVALID_IND) {
			const struct const_Matrix_d* op_ref = constructor_basis(ind_basis,d,cub->rst); // tbd
			const struct const_Matrix_d* op =
				constructor_mm_const_Matrix_d('N','N',1.0,0.0,op_ref,op_t,'R'); // tbd

		} else {
			const struct const_Multiarray_Matrix_d* op_ref = constructor_grad_basis(ind_basis,d,cub->rst); // tbd

// Increment `row` here for d > 1.
		}

//		op->data[row] =

	}

	return op;
}

/** \brief Check if information can be lost while performing the operation.
 *  \return `true` if: `p_i > p_o` or `h_i > h_o`; `false` otherwise. */
bool check_op_info_loss
	(const int*const op_indices)
{
	if ((op_indices[OP_IND_H+1] > op_indices[OP_IND_H]) ||
	    (op_indices[OP_IND_P+1] > op_indices[OP_IND_P]))
		return true;
	return false;
}

// For cv, cc, vv, vc operators
void compute_the_operator_std
()
const int*const op_indices;
const int op_type;
const struct const_Element* element;
{
	const bool info_loss = check_op_info_loss(op_indices);

	const int* h_ptr = &op_indices[OP_IND_H];
	const int* p_ptr = &op_indices[OP_IND_P];

	if (!info_loss) {
		const int*const p_x[2] = { compute_p_cub(p_ptr[0],op_info->cub_type,element->s_type,sim),
		                           compute_p_cub(p_ptr[1],op_info->cub_type,element->s_type,sim), };
// Need new function in cubature.h
		const struct const_Cubature*const cub_o =
			constructor_const_Cubature_hp(element->s_type,h_ptr[0],p_x[0]); // destructed
		const struct const_Cubature*const cub_i =
			constructor_const_Cubature_hp(element->s_type,h_ptr[1],p_x[1]); // destructed

// Construct cv
	// basis_ref == ortho (make function pointer)
		const struct const_Matrix_d* phi_ref_cv_ii = constructor_basis_ref(p_x[1],element->d,cub_i->rst); // tbd
// external function to computed T below.
		const struct const_Matrix_d* phi_cv_ii = NULL;
		if (ortho)
			phi_cv_ii = constructor_copy_const_Matrix_d(phi_ref_cv_ii); // tbd
		else if (lagrange)
			phi_cv_ii = constructor_identity_const_Matrix_d(phi_ref_cv_ii->ext_0); // tbd
		else if (bezier)
			phi_cv_ii = constructor_basis(p_x[1],element->d,cub_i->rst); // tbd
		else
			EXIT_UNSUPPORTED;

		const struct const_Matrix_d* T_i = constructor_sgesv_const_Matrix_d(phi_ref_cv_ii,phi_cv_ii); // tbd

		const struct const_Matrix_d* phi_ref_cv_io = constructor_basis_ref(p_x[1],element->d,cub_o->rst); // tbd
		const struct const_Matrix_d* phi_cv_io     = constructor_mm_const_Matrix_d(phi_ref_cv_io,T_i); // tbd
	// compute

// Update for other operators

	} else {
		EXIT_ADD_SUPPORT; // L2 projection operators.
	}
}


// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

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

/// \brief Compute range (min, max) for the index of `var_type`.
static void compute_range
	(int x_mm[2],                         ///< Set to the min and max values.
	 const struct Operator_Info* op_info, ///< \ref Operator_Info.
	 const char var_type                  ///< Variable type.
	);

/// \brief Compute range (min, max) for `p_o` based on the values of `p_i` and `op_info->p_range`.
static void compute_range_p_o
	(int p_o_mm[2],                       ///< Set to the min and max values.
	 const struct Operator_Info* op_info, ///< \ref Operator_Info.
	 const int p_i                        ///< The input order.
	);

static void set_up_order_and_extents (struct Operator_Info* op_info)
{
	int order_cub = 0;
	int order_op = 0;
	if (op_info->range_d != RANGE_D_0)
		++order_op;

	if (op_info->range_f != RANGE_F_0) {
		++order_cub;
		++order_op;
	}
	++order_cub; // range_h
	++order_op;
	++order_cub; // range_p
	++order_op;
	++order_op;

	ptrdiff_t* extents_cub = malloc(order_cub * sizeof *extents_cub); // keep
	ptrdiff_t* extents_op  = malloc(order_op  * sizeof *extents_op);  // keep
	int ind_cub = 0,
	    ind_op = 0;
	switch (op_info->range_d) {
		case RANGE_D_0:
			/* Do nothing */
			break;
		case RANGE_D_ALL:
			extents_op[ind_op++] = element->d;
			break;
		default:
			EXIT_UNSUPPORTED;
			break;
	}

	switch (op_info->range_f) {
		case RANGE_F_0:
			/* Do nothing */
			break;
		case RANGE_F_ALL:
			extents_cub[ind_cub++] = element->n_f;
			extents_op[ind_op++]   = element->n_f;
			break;
		default:
			EXIT_UNSUPPORTED;
			break;
	}

	switch (op_info->range_h) {
		case RANGE_H_1:
			extents_cub[ind_cub++] = 1;
			extents_op[ind_op++]   = 1;
			break;
		case RANGE_H_ALL:
			extents_cub[ind_cub++] = element->n_ref_max;
			extents_op[ind_op++]   = element->n_ref_max;
			break;
		default:
			EXIT_UNSUPPORTED;
			break;
	}

	switch (op_info->range_p) {
		case RANGE_P_1:
			extents_cub[ind_cub++] = 1;
			extents_op[ind_op++]   = 1;
			extents_op[ind_op++]   = 1;
			break;
		case RANGE_P_PM0: /* fallthrough */
		case RANGE_P_PM1: /* fallthrough */
		case RANGE_P_ALL:
			extents_cub[ind_cub++] = op_info->p_ref[1]+1;
			extents_op[ind_op++]   = op_info->p_ref[1]+1;
			extents_op[ind_op++]   = op_info->p_ref[1]+1;
			break;
		default:
			EXIT_UNSUPPORTED;
			break;
	}

	op_info->order_cub   = order_cub;
	op_info->extents_cub = extents_cub;
	op_info->order_op    = order_op;
	op_info->extents_op  = extents_op;
}

static void set_up_basis_values (struct Operator_Info* op_info)
{
	const int order = op_info->order_op;
	const ptrdiff_t*const extents = op_info->extents_op;

	const ptrdiff_t size = compute_size(order,extents);

	struct Matrix_i* values = constructor_empty_Matrix_i('R',size,OP_ORDER_MAX); // keep

	int d_mm[2],
	    f_mm[2],
	    h_mm[2],
	    p_i_mm[2],
	    p_o_mm[2];

	compute_range(d_mm,op_info,'d');
	compute_range(f_mm,op_info,'f');
	compute_range(h_mm,op_info,'h');
	compute_range(p_i_mm,op_info,'p');

	int row = 0;
	for (int p_i = p_i_mm[0]; p_i <= p_i_mm[1]; ++p_i) {
		compute_range_p_o(p_o_mm,op_info,p_i);
		for (int p_o = p_o_mm[0]; p_o <= p_o_mm[1]; ++p_o) {
		for (int h = h_mm[0]; h <= h_mm[1]; ++h) {
		for (int f = f_mm[0]; f <= f_mm[1]; ++f) {
		for (int d = d_mm[0]; d <= d_mm[1]; ++d) {
			set_row_Matrix_i(row,values,(int[]){d,f,h,p_o,p_i});
			++row;
		}}}}
	}
	values->ext_0 = row;

	op_info->values_op = values;
}

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

// Level 1 ********************************************************************************************************** //

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
		case RANGE_P_1:
			*p_min = 1;
			*p_max = 1;
			break;
		case RANGE_P_PM0: /* fallthrough */
		case RANGE_P_PM1: /* fallthrough */
		case RANGE_P_ALL: /* fallthrough */
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

static void compute_range (int x_mm[2], const struct Operator_Info* op_info, const char var_type)
{
	switch (var_type) {
	case 'd':
		switch (op_info->range_d) {
		case RANGE_D_0:
			x_mm[0] = OP_INVALID_IND;
			x_mm[1] = OP_INVALID_IND;
			break;
		case RANGE_D_ALL:
			x_mm[0] = 0;
			x_mm[1] = op_info->element->d;
			break;
		default:
			EXIT_UNSUPPORTED;
			break;
		}
	case 'f':
		switch (op_info->range_f) {
		case RANGE_F_0:
			x_mm[0] = OP_INVALID_IND;
			x_mm[1] = OP_INVALID_IND;
			break;
		case RANGE_F_ALL:
			x_mm[0] = 0;
			x_mm[1] = op_info->element->n_f;
			break;
		default:
			EXIT_UNSUPPORTED;
			break;
		}
	case 'h':
		switch (op_info->range_h) {
		case RANGE_H_1:
			x_mm[0] = 0;
			x_mm[1] = 0;
			break;
		case RANGE_H_ALL:
			x_mm[0] = 0;
			x_mm[1] = op_info->element->n_ref_max;
			break;
		default:
			EXIT_UNSUPPORTED;
			break;
		}
	case 'p':
		switch (op_info->range_p) {
		case RANGE_P_1:
			x_mm[0] = 1;
			x_mm[1] = 1;
			break;
		case RANGE_P_PM0: /* fallthrough */
		case RANGE_P_PM1: /* fallthrough */
		case RANGE_P_ALL:
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
	case RANGE_P_1:
		p_o_mm[0] = 1;
		p_o_mm[1] = 1;
		break;
	case RANGE_P_PM0:
		p_o_mm[0] = p_i;
		p_o_mm[1] = p_i;
		break;
	case RANGE_P_PM1:
		p_o_mm[0] = max(p_i-1,op_info->p_ref[0]);
		p_o_mm[1] = min(p_i+1,op_info->p_ref[1]);
		break;
	case RANGE_P_ALL:
		x_mm[0] = op_info->p_ref[0];
		x_mm[1] = op_info->p_ref[1];
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}
}

// Level 2 ********************************************************************************************************** //

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
