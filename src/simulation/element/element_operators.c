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


#include "macros.h"
#include "definitions_core.h"
#include "definitions_elements.h"

#include "vector.h"

#include "simulation.h"
#include "element.h"

// Static function declarations ************************************************************************************* //

/** \brief Compute the maximum order required for operators.
 *  \return See brief. */
static int compute_max_p
	(const struct Simulation* sim ///< \ref Simulation.
	);

/** \brief Compute the maximum number of h-refimenent indices required for operators.
 *  \return See brief. */
static int compute_max_h
	(const struct Simulation* sim,       ///< \ref Simulation.
	 const struct const_Element* element ///< \ref const_Element.
	);

/** \brief Constructor for the \ref Multiarray_Cubature::data.
 *  \return Standard. */
struct const_Cubature** constructor_cub_data_array
	(const struct Simulation* sim,        ///< \ref Simulation.
	 const struct const_Element* element, ///< \ref const_Element.
	 const ptrdiff_t order,               ///< \ref Multiarray_Cubature::order.
	 const ptrdiff_t*const extents,       ///< \ref Multiarray_Cubature::extents.
	 const int cub_entity                 ///< The cubature entity for which the nodes will be used.
	);

// Interface functions ********************************************************************************************** //

const struct const_Vector_i* constructor_operator_extents_const_Vector_i
	(const struct Simulation* sim, const struct const_Element* element, const int op_type)
{
	const int n_hp = sim->n_hp;
	          d    = element->d;
	          n_f  = element->n_f;

	int n_ext = 0;
	switch (op_type) {
		case OP_V_D0: n_ext = n_hp;     break;
		case OP_V_D1: n_ext = n_hp+1;   break;
		case OP_F_D0: n_ext = n_hp+1;   break;
		case OP_F_D1: n_ext = n_hp+2;   break;
		default:      EXIT_UNSUPPORTED; break;
	}

	int* extents = malloc(n_ext * sizeof *extents); // moved
	switch (op_type) {
	case OP_V_D0:
		// Do nothing
		break;
	case OP_V_D1:
		extents[0] = d;
		break;
	case OP_F_D0:
		extents[0] = n_f;
		break;
	case OP_F_D0:
		extents[0] = d;
		extents[1] = n_f;
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}

	int* extents_tail = NULL;
	switch (sim->adapt_type) {
	case ADAPT_0: {
		extents_tail = (int[]) {1,1,1};
		break;
	} case ADAPT_P: {
		const int max_pp1 = compute_max_p(sim)+1;
		extents_tail = (int[]) {1,max_pp1,max_pp1};
		break;
	} case ADAPT_H: {
		const int max_h = compute_max_h(sim,element);
		extents_tail = (int[]) {max_h,1,1}; // Needs multiplier for face operators.
		break;
	} case ADAPT_HP: {
		const int max_pp1 = compute_max_p(sim)+1,
		          max_h   = compute_max_h(sim,element);
		extents_tail = (int[]) {max_h,max_pp1,max_pp1}; // Needs multiplier for face operators.
		break;
	} default:
		EXIT_UNSUPPORTED;
		break;
	}

	for (int i = n_ext-n_hp; i < n_ext; ++i)
		extents[i] = extents_tail[i];

	return constructor_move_const_Vector_i_i(n_ext,true,extents);
}

const struct const_Multiarray_Cubature* constructor_const_Multiarray_Cubature
	(const struct Simulation* sim, const struct const_Element* element, const struct const_Vector_i* ext_v1_V,
	 const int cub_entity, const int n_skip)
{
	const ptrdiff_t ext_0 = ext_v1_V->ext_0,
	                order = ext_0-n_skip-1;

	ptrdiff_t* extents = malloc(order * sizeof *extents); // keep
	for (ptrdiff_t ind = 0, i = ext_0-order; i < ext_0; ++i) {
		if (i != ext_0-1) // Does not use the p_in value.
			extents[ind++] = ext_v1_V[i];
	}

	struct const_Cubature** cub_data = constructor_cub_data_array(sim,element,order,extents,cub_entity); // keep


	struct Multiarray_Cubature* cub_MA = malloc(sizeof *cub_MA); // returned

	cub_MA->order   = order;
	cub_MA->extents = extents;
	cub_MA->extents = owns_data;
	cub_MA->data    = cub_data;

	return (const struct const_Multiarray_Cubature*) cub_MA;
}

void destructor_const_Multiarray_Cubature (const struct const_Multiarray_Cubature*const a)
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

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Compute the index of computational element associated with the cubature entity.
 *  \return See brief. */
static int compute_computational_element
	(const int cub_entity ///< The cubature entity.
	);

static int compute_max_p (const struct Simulation* sim)
{
	return sim->p_s_v[1];
}

static int compute_max_h (const struct Simulation* sim, const struct const_Element* element)
{
	const bool h_adapt = is_h_adaptive(sim);

	switch (element->type) {
	case LINE:
		return ( h_adapt ? NREFMAXLINE : 1 );
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}
}

static void compute_p_range
	(int*const p_min, int*const p_max, const int cub_ent_val, const int n_p, const struct Simulation* sim)
{
// Likely need an indicator of which sim p_*_* value to base this computation off of in future.
	*p_min = compute_p_entity(sim->p_s_v[0],
	*p_max = ( n_p == 1 ? *p_min : compute_p_entity(sim->p_s_v[1]));
}

struct const_Cubature** constructor_cub_data_array
	(const struct Simulation* sim, const struct const_Element* element, const ptrdiff_t order,
	 const ptrdiff_t*const extents, const int cub_entity)
{
	const ptrdiff_t size = compute_size(order,extents);
	struct const_Cubature** cub_data = malloc(size * sizeof *cub_data); // keep

	const int comp_elem = compute_computational_element(cub_entity);

	ptrdiff_t cub_values[order];
	cub_values[0] =

	int p_min = -1,
	    p_max = -1;
	if (comp_elem == CUB_ENT_V) {
		assert(order == 2);
		if (extents[0] != 1)
			EXIT_ADD_SUPPORT; // sub-elements

		compute_p_range(&p_min,&p_max,cub_entity%CUB_ENT_MULT,extents[1],sim);
		if (extents[1] == 1) {
			p_min = sim->p_s_
		}





	} else if (comp_elem == CUB_ENT_F) {
		assert(order == 3);
		EXIT_ADD_SUPPORT;
	} else {
		EXIT_UNSUPPORTED;
	}

//	const struct const_Cubature* d2_p3_WSH, ///< See \ref definitions_cubature.h.
//		cub_data->d2_p3_WSH = constructor_const_Cubature_si(2,3,CUB_WSH);  // keep
}

// Level 1 ********************************************************************************************************** //

static int compute_computational_element (const int cub_entity)
{
	if (cub_entity/CUB_ENT_MULT == CUB_ENT_V/CUB_ENT_MULT)
		return CUB_ENT_V;
	else if (cub_entity/CUB_ENT_MULT == CUB_ENT_F/CUB_ENT_MULT)
		return CUB_ENT_F;
	else
		EXIT_UNSUPPORTED;
}
