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
#include "cubature.h"

// Static function declarations ************************************************************************************* //

/** \brief Constructor for the \ref Multiarray_Cubature::data.
 *  \return Standard. */
const struct const_Cubature** constructor_cub_data_array
	(const struct Simulation* sim,        ///< \ref Simulation.
	 const struct const_Element* element, ///< \ref const_Element.
	 const ptrdiff_t order,               ///< \ref Multiarray_Cubature::order.
	 const ptrdiff_t*const extents,       ///< \ref Multiarray_Cubature::extents.
	 const struct Operator_Info* op_info  ///< \ref Operator_Info.
	);

/// \brief Destructor for a \ref Multiarray_Cubature\* container.
static void destructor_Multiarray_Cubature
	(struct Multiarray_Cubature* a ///< Standard.
	);

// Interface functions ********************************************************************************************** //

struct Operator_Info* constructor_Operator_Info
	(const int range_d, const int range_f, const int range_p, const int range_h, const int cub_type,
	 const int p_ref[2])
{
	struct Operator_Info* op_info = malloc(sizeof *op_info); // returned

	const_cast_i(&op_info->range_d,range_d);
	const_cast_i(&op_info->range_f,range_f);
	const_cast_i(&op_info->range_h,range_h);
	const_cast_i(&op_info->range_p,range_p);

	const_cast_i(&op_info->cub_type,cub_type);
	const_cast_i1(op_info->p_ref,p_ref,2);

	return op_info;
}

void destructor_Operator_Info (struct Operator_Info* op_info)
{
	free(op_info);
}

const struct const_Multiarray_Cubature* constructor_const_Multiarray_Cubature
	(const struct Simulation* sim, const struct const_Element* element, const struct Operator_Info* op_info)
{
	// Compute order and extents.
	int order = 0;
	if (op_info->range_f != RANGE_F_0)
		++order;
	++order; // range_h
	++order; // range_p

	ptrdiff_t* extents = malloc(order * sizeof *extents); // keep
	int ind = 0;
	switch (op_info->range_f) {
		case RANGE_F_0:   /* Do nothing */               break;
		case RANGE_F_ALL: extents[ind++] = element->n_f; break;
		default:          EXIT_UNSUPPORTED;              break;
	}

	switch (op_info->range_h) {
		case RANGE_H_1:   extents[ind++] = 1;                  break;
		case RANGE_H_ALL: extents[ind++] = element->n_ref_max; break;
		default:          EXIT_UNSUPPORTED;                    break;
	}

	switch (op_info->range_p) {
		case RANGE_P_1:   extents[ind++] = 1;                   break;
		case RANGE_P_PM0: extents[ind++] = 1;                   break;
		case RANGE_P_PM1: /* fallthrough */
		case RANGE_P_ALL: extents[ind++] = op_info->p_ref[1]+1; break;
		default:          EXIT_UNSUPPORTED;                     break;
	}

	// Compute \ref Cubature data.
	const struct const_Cubature** cub_data = constructor_cub_data_array(sim,element,order,extents,op_info); // keep

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
	(const struct Simulation* sim, const struct const_Element* element, const ptrdiff_t order,
	 const ptrdiff_t*const extents, const struct Operator_Info* op_info)
{
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
