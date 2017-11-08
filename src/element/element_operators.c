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
#include "definitions_core.h"
#include "definitions_element_operators.h"
#include "definitions_elements.h"
#include "definitions_bases.h"

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
	(ptrdiff_t*const ind_values,           ///< Defined for \ref set_operator_fptr.
	 const struct Multiarray_Operator* op, ///< Defined for \ref set_operator_fptr.
	 const struct Operator_Info* op_info,  ///< Defined for \ref set_operator_fptr.
	 const struct Simulation* sim          ///< Defined for \ref set_operator_fptr.
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
int convert_to_type
	(const char* name_type ///< Defined for \ref constructor_Operator_Info.
	);

/** \brief Convert the `char*` input to the appropriate definition of OP_R_D_*.
 *  \return See brief. */
int convert_to_range_d
	(const char* name_type ///< Defined for \ref constructor_Operator_Info.
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
	(const char type_range,      ///< The type of range parameter. Options: 'h', 'p'.
	 const char*const name_range ///< Defined for \ref constructor_Operator_Info.
	);

/** \brief Set Op_IO::p_rel for each of the \ref Operator_Info::op_io members.
 *  The relative degree is set based on the values of the p_\*_\*_p members in \ref Simulation corresponding to the
 *  computational element and kind of the input basis of the operator. */
static void set_op_info_p_rel
	(struct Operator_Info* op_info, ///< \ref Operator_Info.
	 const struct Simulation* sim   ///< \ref Simulation.
	);

/// \brief Set up \ref Operator_Info extents_* members.
static void set_up_extents
	(struct Operator_Info* op_info ///< \ref Operator_Info.
	);

/// \brief Set up \ref Operator_Info::values_op.
static void set_up_values_op
	(struct Operator_Info* op_info ///< \ref Operator_Info.
	);

// Interface functions ********************************************************************************************** //

const struct Multiarray_Operator* constructor_operators
	(const char*const name_type, const char*const name_in, const char*const name_out, const char*const name_range,
	 const int p_ref[2], const struct const_Element* element, const struct Simulation* sim)
{
	struct Operator_Info* op_info =
		constructor_Operator_Info(name_type,name_in,name_out,name_range,p_ref,element,sim); // destructed
	assert((op_info->range_d == OP_R_D_0) || (op_info->range_d == OP_R_D_ALL));

	const struct Multiarray_Operator* op = constructor_empty_Multiarray_Operator_V(op_info->extents_op); // returned

	const ptrdiff_t row_max = op_info->values_op->ext_0;
	for (ptrdiff_t row = 0; row < row_max; ) // row is incremented when setting the operators.
		op_info->set_operator(&row,op,op_info,sim);
	destructor_Operator_Info(op_info);

	return op;
}

const struct Multiarray_Operator* constructor_operators2
	(const char*const name_type, const char*const name_in, const char*const name_out, const char*const name_range,
	 const struct const_Element* element, const struct Simulation* sim)
{
	struct Operator_Info* op_info =
		constructor_Operator_Info2(name_type,name_in,name_out,name_range,element,sim); // destructed
	assert((op_info->range_d == OP_R_D_0) || (op_info->range_d == OP_R_D_ALL));

	const struct Multiarray_Operator* op = constructor_empty_Multiarray_Operator_V(op_info->extents_op); // returned

	const ptrdiff_t row_max = op_info->values_op->ext_0;
	for (ptrdiff_t row = 0; row < row_max; ) // row is incremented when setting the operators.
		op_info->set_operator(&row,op,op_info,sim);
	destructor_Operator_Info(op_info);

	return op;
}

const struct const_Multiarray_Vector_i* constructor_operators_nc
	(const int ind_f_elem, const char*const name_in, const char*const name_out, const char*const name_range,
	 const int p_ref[2], const struct const_Element* element, const struct Simulation* sim)
{
	struct Operator_Info* op_info =
		constructor_Operator_Info("UNUSED0",name_in,name_out,name_range,p_ref,element,sim); // destructed

	const struct const_Vector_i* e_o = op_info->extents_op;

	const struct const_Element* f_element = element->face_element[ind_f_elem];
	assert(f_element != NULL);

	// Add an additional index for the possible permutations
	const ptrdiff_t len_e = e_o->ext_0+1;
	struct Vector_i* e_o_p1 = constructor_empty_Vector_i(len_e); // destructed
	e_o_p1->data[0] = get_n_perm_corr(f_element->d,f_element->s_type);
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
	struct Operator_Info* op_info =
		constructor_Operator_Info("UNUSED0",name_in,name_out,name_range,p_ref,element,sim); // destructed

	const struct const_Multiarray_Vector_d* w =
		constructor_empty_const_Multiarray_Vector_d_V(false,op_info->extents_op); // returned

	const ptrdiff_t row_max = op_info->values_op->ext_0;
	for (ptrdiff_t row = 0; row < row_max; ++row)
		set_weights(row,w,op_info,sim);
	destructor_Operator_Info(op_info);

	return w;
}

struct Operator_Info* constructor_Operator_Info
	(const char*const name_type, const char*const name_in, const char*const name_out, const char*const name_range,
	 const int p_ref[2], const struct const_Element* element, const struct Simulation* sim)
{
UNUSED(p_ref);
	struct Operator_Info* op_info = malloc(sizeof *op_info); // returned

	op_info->element = element;

	const int op_type = convert_to_type(name_type);

	const_cast_i(&op_info->op_type,op_type);

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

	const_cast_i1(op_info->p_ref,sim->p_ref,2);
	for (int i = 0; i < 2; ++i)
		const_cast_i(&op_info->op_io[i].p_rel,0);

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

struct Operator_Info* constructor_Operator_Info2
	(const char*const name_type, const char*const name_in, const char*const name_out, const char*const name_range,
	 const struct const_Element* element, const struct Simulation* sim)
{
	struct Operator_Info* op_info = malloc(sizeof *op_info); // returned

	op_info->element = element;

	const int op_type = convert_to_type(name_type);

	const_cast_i(&op_info->op_type,op_type);

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

	const_cast_i1(op_info->p_ref,sim->p_ref,2);
	set_op_info_p_rel(op_info,sim);
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
	case 'g': // geometry
		if (op_sc == 's')
			return 1;

		const_cast_c(&op_io->kind,'s');
		const int p_s = compute_p_basis(op_io,sim);
		const_cast_c(&op_io->kind,'g');

		if (strcmp(sim->geom_rep,"isoparametric") == 0)
			return p_s;
		else if (strcmp(sim->geom_rep,"superparametric") == 0)
			return p_s+1;
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
	case 'p': // plotting
		return p_op;
		break;
	} case 'c': // fallthrough
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

bool check_op_info_loss (const int*const op_values)
{
	if ((op_values[OP_IND_H+1] > op_values[OP_IND_H]) ||
	    (op_values[OP_IND_P+1] > op_values[OP_IND_P]))
		return true;
	return false;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Compute the super type of the nodes based on the kind of operator.
 *  \return See brief. */
static int compute_super_type_op
	(const struct Op_IO* op_io,          ///< \ref Op_IO.
	 const struct const_Element* element ///< \ref const_Element.
	);

/** \brief Constructor for the coefficient to value operator, selecting the operator basis based on `kind`.
 *  \return See brief. */
static const struct const_Matrix_d* constructor_cv
	(const struct const_Matrix_d* rst, ///< The rst coordinates.
	 const struct Op_IO* op_io,        ///< \ref Op_IO.
	 const struct Simulation* sim      ///< \ref Simulation.
	);

/** \brief Constructor for a \ref Multiarray_Matrix_d\* of order 1 with `owns_data = false` to hold locally computed
 *         operators.
 *  \return See brief. */
static const struct const_Multiarray_Matrix_d* constructor_op_MMd
	(const bool owns_data, ///< Defined in \ref Multiarray_Matrix_d.
	 const ptrdiff_t ext_0 ///< The size of the single extent.
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

/** \brief Get the value for the maximum number of h-refinements.
 *  \return See brief. */
static int get_n_ref_max
	(const struct Operator_Info* op_info ///< \ref Operator_Info.
	);

/** \brief Return the number of operators to be set depending on \ref Operator_Info::range_d.
 *  \return See brief. */
static int get_n_op
	(const int range_d, ///< \ref Operator_Info::range_d.
	 const int d        ///< The dimension of the element for which operators are being constructed.
	);

/// \brief Set the values of \ref Op_IO for the current operator.
static void set_current_op_io
	(const struct Operator_Info* op_info, ///< \ref Operator_Info.
	 const int* op_values                 ///< One line of values of \ref Operator_Info::values_op.
	);

/** \brief Return the index of the operator currently being set.
 *  \return See brief. */
static ptrdiff_t get_ind_op
	(const struct Operator_Info* op_info, ///< \ref Operator_Info.
	 const int* op_values,                ///< One line of values of \ref Operator_Info::values_op.
	 const struct Multiarray_Operator* op ///< Defined for \ref set_operator_fptr.
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

	const bool info_loss = check_op_info_loss(op_values);
	const struct const_Element* element = op_info->element;

	const int n_op = get_n_op (op_info->range_d,element->d);

	const struct const_Multiarray_Matrix_d* op_ioN = constructor_op_MMd(false,n_op); // destructed
	if (!info_loss) {
// separate function in future
		set_current_op_io(op_info,op_values);
		const struct Op_IO* op_io = op_info->op_io;

		const int p_i = compute_p_basis(&op_io[OP_IND_I],sim);

		const int s_type = op_io[OP_IND_I].s_type;
		constructor_basis_fptr constructor_basis = get_constructor_basis_by_super_type(s_type,"orthonormal");

		/* Compute `op_cvNr`: [operator, coefficients to values, differentiation order N, reference basis].
		 * op_cvNr == basis(p_i,rst_o). */
		const struct const_Nodes* nodes_o = constructor_const_Nodes_h(OP_IND_O,op_io,element,sim); // destructed

		const struct const_Multiarray_Matrix_d* op_cvNr = NULL;
		if (op_info->range_d == OP_R_D_0) {
			const struct const_Matrix_d* cv0r = constructor_basis(p_i,nodes_o->rst); // moved

			op_cvNr = constructor_op_MMd(true,n_op); // destructed
			const_constructor_move_const_Matrix_d(&op_cvNr->data[0],cv0r); // destructed
		} else if (op_info->range_d == OP_R_D_ALL) {
			constructor_grad_basis_fptr constructor_grad_basis =
				get_constructor_grad_basis_by_super_type(s_type,"orthonormal");
			op_cvNr = constructor_grad_basis(p_i,nodes_o->rst); // destructed
		} else {
			EXIT_ERROR("Unsupported: %d\n",op_info->range_d);
		}

		/* Compute `op_cvN`: [operator, coefficients to values, differentiation order N].
		 * op_cvN == op_cvNr*T_i (Corollary 2.2, \cite Zwanenburg2016). */
		const struct const_Nodes* nodes_i = constructor_const_Nodes_h(OP_IND_I,op_io,element,sim); // destructed

		const struct const_Matrix_d* cv0r_ii = constructor_basis(p_i,nodes_i->rst);   // destructed
		const struct const_Matrix_d* cv0_ii  = constructor_cv(nodes_i->rst,&op_io[OP_IND_I],sim); // destructed
		destructor_const_Nodes(nodes_i);

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
			const struct const_Matrix_d* cv0_oo = constructor_cv(nodes_o->rst,&op_io[OP_IND_O],sim); // destructed

			for (int i = 0; i < n_op; ++i) {
				const struct const_Matrix_d* cv0 = op_cvN->data[i];
			const struct const_Matrix_d* co0 = constructor_sgesv_const_Matrix_d(cv0_oo,cv0); // moved
				const_constructor_move_const_Matrix_d(&op_coN->data[i],co0); // destructed
			}
			destructor_const_Matrix_d(cv0_oo);
		} else {
			const_cast_bool(&op_cvN->owns_data,false);
			for (int i = 0; i < n_op; ++i)
				const_constructor_move_const_Matrix_d(&op_coN->data[i],op_cvN->data[i]); // destructed
		}
		destructor_const_Nodes(nodes_o);
		destructor_const_Multiarray_Matrix_d(op_cvN);

		/* Compute `op_ioN`: [operator, input (coefs/values) to output (coefs/values), differentiation order N].
		 * op_ioN == op_coN*inv_cv0_ii. */
		if (op_type == OP_T_VV || op_type == OP_T_VC) {
			const struct const_Matrix_d* inv_cv0_ii = constructor_inverse_const_Matrix_d(cv0_ii); // destructed
			for (int i = 0; i < n_op; ++i) {
				const struct const_Matrix_d* co0 = op_coN->data[i];
				const struct const_Matrix_d* io0 = constructor_mm_NN1R_const_Matrix_d(co0,inv_cv0_ii); // moved
				const_constructor_move_const_Matrix_d(&op_ioN->data[i],io0); // moved
			}
			destructor_const_Matrix_d(inv_cv0_ii);
		} else {
			const_cast_bool(&op_coN->owns_data,false);
			for (int i = 0; i < n_op; ++i)
				const_constructor_move_const_Matrix_d(&op_ioN->data[i],op_coN->data[i]); // moved
		}
		destructor_const_Matrix_d(cv0_ii);
		destructor_const_Multiarray_Matrix_d(op_coN);
	} else {
		EXIT_ADD_SUPPORT; // L2 projection operators.
	}

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
		scale_Matrix_by_Vector_d('R',1.0,op_std,nodes_o->w,false);
	}
	destructor_const_Nodes(nodes_o);

	if (sim->collocated) {
		const struct const_Nodes* nodes_i = constructor_const_Nodes_h(OP_IND_I,op_io,element,sim); // destructed
		assert(nodes_i->has_weights);

		for (int i = 0; i < n_op; ++i) {
			struct Matrix_d* op_std = (struct Matrix_d*) op->data[ind_op+i]->op_std;
			scale_Matrix_by_Vector_d('L',1.0,op_std,nodes_i->w,true);
		}
		destructor_const_Nodes(nodes_i);
		EXIT_ERROR("Ensure that all is working as expected.\n");
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

int convert_to_type (const char* name_type)
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

int convert_to_range_d (const char* name_type)
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

int convert_to_range_ce (const char ce_i, const char ce_o)
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
		if (strstr(name_range,"P_1P"))
			if (strstr(name_range,"P_1PPM1"))
				return OP_R_P_1PPM1;
			else
				return OP_R_P_1P;
		else if (strstr(name_range,"P_1"))
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

static void set_op_info_p_rel (struct Operator_Info* op_info, const struct Simulation* sim)
{
	const char ce_i   = op_info->op_io[OP_IND_I].ce,
	           kind_i = op_info->op_io[OP_IND_I].kind;

	int p_rel = -1;
	switch (ce_i) {
	case 'v':
		switch (kind_i) {
		case 'g': // fallthrough
		case 'm': // fallthrough
		case 's':
			p_rel = sim->p_s_v_p;
			break;
		case 'r':
			p_rel = sim->p_sg_v_p;
			break;
		default:
			EXIT_ERROR("Unsupported: %c\n",kind_i);
			break;
		}
		break;
	case 'f':
		switch (kind_i) {
		case 'f':
			p_rel = sim->p_s_f_p;
			break;
		case 't':
			p_rel = sim->p_sg_f_p;
			break;
		default:
			EXIT_ERROR("Unsupported: %c\n",kind_i);
			break;
		}
		break;
	default:
		EXIT_ERROR("Unsupported: %c\n",ce_i);
		break;
	}

	for (int i = 0; i < 2; ++i)
		const_cast_i(&op_info->op_io[i].p_rel,p_rel);
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
		} default:
			EXIT_UNSUPPORTED;
			break;
	}

	switch (op_info->range_p) { // p_o, p_i
		case OP_R_P_1:
			push_back_Vector_i(extents_op,2,false,false);
			push_back_Vector_i(extents_op,2,false,false);
			break;
		case OP_R_P_1P:    // fallthrough
		case OP_R_P_1PPM1:
			push_back_Vector_i(extents_op,op_info->p_ref[1]+1,false,false);
			push_back_Vector_i(extents_op,2,false,false);
			break;
		case OP_R_P_PM0: // fallthrough
		case OP_R_P_PM1: // fallthrough
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
		for (int ce_i = ce_i_mm[0]; ce_i < ce_i_mm[1]; ++ce_i) {
		for (int ce_o = ce_o_mm[0]; ce_o < ce_o_mm[1]; ++ce_o) {
		for (int d = d_mm[0]; d < d_mm[1]; ++d) {
			set_row_Matrix_i(row,values,(int[]){d,ce_o,ce_i,h_o,h_i,p_o,p_i});
			++row;
		}}}}}}
	}
	values->ext_0 = row;

	op_info->values_op = (const struct const_Matrix_i*) values;
//print_const_Matrix_i(op_info->values_op);
}

// Level 1 ********************************************************************************************************** //

static int compute_super_type_op (const struct Op_IO* op_io, const struct const_Element* element)
{
	const int sub_e_type = compute_elem_type_sub_ce(element->type,op_io->ce,op_io->h_op);
	return compute_super_from_elem_type(sub_e_type);
}

static const struct const_Matrix_d* constructor_cv
	(const struct const_Matrix_d* rst, const struct Op_IO* op_io, const struct Simulation* sim)
{
	int basis_type = -1;
	switch (op_io->kind) {
	case 'g': // fallthrough
	case 'm':
		basis_type = get_basis_i_from_s(sim->basis_geom);
		break;
	case 's': // fallthrough
	case 'f':
		basis_type = get_basis_i_from_s(sim->basis_sol);
		break;
	default:
		EXIT_ERROR("Unsupported: %c.\n",op_io->kind);
		break;
	}

	const int s_type = op_io->s_type;
	const int p_basis = compute_p_basis(op_io,sim);

	const struct const_Matrix_d* cv = NULL;
	if (basis_type == BASIS_LAGRANGE) {
		cv = constructor_identity_const_Matrix_d('R',rst->ext_0); // returned
	} else {
		constructor_basis_fptr constructor_basis = get_constructor_basis_by_super_type_i(s_type,basis_type);
		cv = constructor_basis(p_basis,rst); // returned
	}

	return cv;
}

static const struct const_Multiarray_Matrix_d* constructor_op_MMd (const bool owns_data, const ptrdiff_t ext_0)
{
	struct Multiarray_Matrix_d* op_MMd = constructor_empty_Multiarray_Matrix_d(false,1,&ext_0); // returned
	op_MMd->owns_data = owns_data;

	return (const struct const_Multiarray_Matrix_d*) op_MMd;
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
		default:
			EXIT_UNSUPPORTED;
			break;
		}
		break;
	case 'p':
		// In general, p_o depends on p_i and requires more information than that provided here.
		assert(var_io == 'i');

		switch (op_info->range_p) {
		case OP_R_P_1: // fallthrough
		case OP_R_P_1P:
		case OP_R_P_1PPM1:
			x_mm[0] = 1;
			x_mm[1] = 1+1;
			break;
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
	case OP_R_P_1:
		p_o_mm[0] = 1;
		p_o_mm[1] = 1+1;
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

static int get_n_ref_max (const struct Operator_Info* op_info)
{
	// There are currently no h-adaptation operators for other `range_ce` values.
	if (op_info->range_ce == OP_R_CE_VV)
		return op_info->element->n_ref_max_v;
	else if (op_info->range_ce == OP_R_CE_VF)
		return op_info->element->n_ref_max_f;
	EXIT_UNSUPPORTED;
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
		const_cast_i(&op_io[i].s_type,compute_super_type_op(&op_io[i],op_info->element));
	}
}

static int get_n_op (const int range_d, const int d)
{
	return ( range_d == OP_R_D_0 ? 1 : d );
}

static ptrdiff_t get_ind_op
	(const struct Operator_Info* op_info, const int* op_values, const struct Multiarray_Operator* op)
{
	const int order_op = op_info->extents_op->ext_0;
	const struct const_Vector_i* indices_op = constructor_indices_Vector_i(order_op,op_values,NULL); // destructed

	ptrdiff_t ind_op = compute_index_sub_container_pi(op->order,0,op->extents,indices_op->data);
	destructor_const_Vector_i(indices_op);

	return ind_op;
}

static ptrdiff_t get_ind_nc
	(const struct Operator_Info* op_info, const int* op_values, const struct const_Multiarray_Vector_i* nc)
{
	const int order_op = op_info->extents_op->ext_0;
	const struct const_Vector_i* indices_op = constructor_indices_Vector_i(order_op,op_values,NULL); // destructed

	ptrdiff_t ind_nc = compute_index_sub_container_pi(nc->order,1,nc->extents,indices_op->data);
	destructor_const_Vector_i(indices_op);

	return ind_nc;
}

static ptrdiff_t get_ind_w
	(const struct Operator_Info* op_info, const int* op_values, const struct const_Multiarray_Vector_d* w)
{
	const int order_op = op_info->extents_op->ext_0;
	const struct const_Vector_i* indices_op = constructor_indices_Vector_i(order_op,op_values,NULL); // destructed

	ptrdiff_t ind_w = compute_index_sub_container_pi(w->order,0,w->extents,indices_op->data);
	destructor_const_Vector_i(indices_op);

	return ind_w;
}
