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

#include "cubature_operators.h"
#include "element_operators.h"

#include <assert.h>
#include <string.h>
#include "gsl/gsl_math.h"

#include "macros.h"
#include "definitions_cubature.h"
#include "definitions_element_operators.h"

#include "matrix.h"
#include "vector.h"

#include "simulation.h"
#include "element.h"
#include "bases.h"
#include "cubature.h"
#include "element_operators.h"
#include "const_cast.h"

// Static function declarations ************************************************************************************* //

/** \brief Compute the node type of the cubature based on the kind of operator.
 *  \return See brief. */
static int compute_node_type_cub
	(const struct Op_IO* op_io,           ///< \ref Op_IO.
	 const struct const_Element* element, ///< \ref const_Element.
	 const struct Simulation* sim         ///< \ref Simulation.
	);

/** \brief Compute the dimension of the cubature based off of the computational element and the simulation dimension.
 *  \return See brief. */
static int compute_d_cub
	(const char ce,   ///< \ref Op_IO::ce.
	 const int elem_d ///< \ref Element::d.
	);

/** \brief Compute the order of the cubature nodes to be computed based on the reference order and the kind of operator.
 *  \return See brief. */
static int compute_p_cub
	(const struct Op_IO* op_io,   ///< \ref Op_IO.
	 const int node_type,         ///< \ref Cubature::node_type.
	 const struct Simulation* sim ///< \ref Simulation.
	);

/** \brief Constructor for the rst coordinates associated with vertices of the (potentially h-refined) reference
 *         element.
 *  \return See brief. */
const struct const_Matrix_d* constructor_rst_ve
	(const int s_type,            ///< \ref Element::s_type.
	 const int d_i,               ///< The dimension of the input basis coordinates.
	 const int d_io,              ///< The dimension of the input/output rst coordinates.
	 const int ind_h,             ///< The h-refinement index.
	 const char ce,               ///< \ref Op_IO::ce.
	 const struct Simulation* sim ///< \ref Simulation.
	);

/** \brief Constructor for rst coordinates as a projection of
 *  \return See brief. */
struct Matrix_d* constructor_rst_proj
	(const int e_type,                      ///< \ref Element::type.
	 const int d_i,                         ///< The dimension of the input basis coordinates.
	 const int d_io,                        ///< The dimension of the input/output rst coordinates.
	 const int ind_h,                       ///< The h-refinement index.
	 const char ce,                         ///< \ref Op_IO::ce.
	 const struct const_Matrix_d* b_coords, ///< The barycentric coordinates of the element from which to project.
	 const struct Simulation* sim           ///< \ref Simulation.
	);

// Constructor functions ******************************************************************************************** //

// Possibly return a Cubature_Multiarray here.
const struct const_Cubature* constructor_const_Cubature_h
	(const int ind_io, const struct Op_IO op_io[2], const struct const_Element* element,
	 const struct Simulation* sim)
{
	const int s_type_io    = op_io[ind_io].s_type;
	const int node_type_io = compute_node_type_cub(&op_io[ind_io],element,sim);

	// Always use the input to establish the dimension of the cubature nodes as it is in bases of this dimension in
	// which the nodes will be used.
	const int d_i          = compute_d_cub(op_io[OP_IND_I].ce,element->d);
	const int d_io         = compute_d_cub(op_io[ind_io].ce,element->d);
	const int p_io         = compute_p_cub(&op_io[ind_io],node_type_io,sim);
	const int h_io         = op_io[ind_io].h_op;
	const char ce_io       = op_io[ind_io].ce;

	cubature_fptr constructor_Cubature = get_cubature_by_super_type(s_type_io);
	basis_fptr    constructor_basis    = get_basis_by_super_type(s_type_io,"ortho");

	const struct const_Cubature* cub_io = constructor_Cubature(d_io,p_io,node_type_io);

	const struct const_Matrix_d* rst_ve = constructor_rst_ve(s_type_io,d_i,d_io,0,ce_io,sim); // tbd

	// vXX: (v)olume XX (Arbitrary, would be replaced with [op_io.kind,op_io.sc] if specified)
	const struct const_Matrix_d* phi_vvs_vXX     = constructor_basis(1,cub_io->rst); // tbd
	const struct const_Matrix_d* phi_vvs_vvs     = constructor_basis(1,rst_ve); // tbd
	const struct const_Matrix_d* phi_vvs_vvs_inv = constructor_inverse_const_Matrix_d(phi_vvs_vvs); // tbd
	const struct const_Matrix_d* cv0_vvs_vXX =
		constructor_mm_const_Matrix_d('N','N',1.0,0.0,phi_vvs_vXX,phi_vvs_vvs_inv,'R'); // tbd

	struct Cubature* cubature = malloc(sizeof* cubature); // returned

	cubature->has_weights = false;
	cubature->w           = NULL;

	cubature->p           = p_io;
	cubature->node_type   = node_type_io;

	const char ce_i = op_io[OP_IND_I].ce,
	           ce_o = op_io[OP_IND_O].ce;
	if ((ind_io == OP_IND_I) || (ce_i == 'v')) { // ((vv || ff || ee) || (vf || ve))
		const struct const_Matrix_d* rst_ve_io =
			constructor_rst_ve(s_type_io,d_i,d_io,h_io,ce_io,sim); // destructed
		cubature->rst = constructor_mm_Matrix_d('N','N',1.0,0.0,cv0_vvs_vXX,rst_ve_io,'C'); // keep
		destructor_const_Matrix_d(rst_ve_io);
	} else if (ce_o == 'v') { // (fv || ev)
		cubature->rst = constructor_rst_proj(element->type,d_i,d_io,h_io,ce_i,cv0_vvs_vXX,sim); // keep
	} else { // (fe || ef)
		EXIT_ADD_SUPPORT;
	}

printf("ccCh: rst\n");
print_Matrix_d(cubature->rst,1e-15);

	return (const struct const_Cubature*) cubature;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Compute the node type of the cubature based on the kind of operator for standard (non-collocated) nodes.
 *  \return See brief. */
static int compute_node_type_cub_std
	(const struct Op_IO* op_io,           ///< \ref Op_IO.
	 const struct const_Element* element, ///< \ref const_Element.
	 const struct Simulation* sim         ///< \ref Simulation.
	);

/** \brief Compute the node type of the cubature based on the kind of operator for collocated nodes.
 *  \return See brief. */
static int compute_node_type_cub_collocated
	(const struct Op_IO* op_io,           ///< \ref Op_IO.
	 const struct const_Element* element, ///< \ref const_Element.
	 const struct Simulation* sim         ///< \ref Simulation.
	);

/** \brief Constructor for the barycentric coordinates of the p2 reference element of the input type.
 *  \return See brief.
 *
 *  \todo Modify the b_coords such that they correspond to the documentation from python.
 *  \todo Move this to be with the h-refinement related connectivity information?
 */
const struct const_Matrix_d* constructor_b_coords
	(const int e_type ///< \ref Element::type.
	);

/** \brief Constructor for the indices of the barycentric coordinates of the p2 reference element of the input type to
 *         be used.
 *
 *  \todo Add images of the reference elements with associated coordinates and volume/face/edge h_ref indices, the nodes
 *  correspond to those used in the python documentation file.
 *  \todo Move this to be with the h-refinement related connectivity information?
 *
 *  \return See brief. */
const struct const_Vector_i* constructor_ind_h_b_coords
	(const int e_type,            ///< \ref Element::type.
	 const char ce,               ///< \ref Op_IO::ce.
	 const int ind_h,             ///< The index of the h-refinement.
	 const struct Simulation* sim ///< \ref Simulation.
	);

static int compute_node_type_cub
	(const struct Op_IO* op_io, const struct const_Element* element, const struct Simulation* sim)
{
	switch (op_io->kind) {
	case 's': // fallthrough
	case 'g': // fallthrough
	case 'm': // fallthrough
		return compute_node_type_cub_std(op_io,element,sim);
		break;
	case 'c':
		if (!sim->collocated)
			return compute_node_type_cub_std(op_io,element,sim);
		else
			return compute_node_type_cub_collocated(op_io,element,sim);
		break;
	default:
		EXIT_ERROR("Unsupported: %c\n",op_io->kind);
		break;
	}
	EXIT_UNSUPPORTED;
	return -1;
}

static int compute_d_cub (const char ce, const int elem_d)
{
	assert(elem_d > 0);

	switch (ce) {
		case 'v': return elem_d;              break;
		case 'f': return elem_d-1;            break;
		case 'e': return GSL_MAX(elem_d-2,0); break;
		default:
			EXIT_ERROR("Unsupported: %c\n",ce);
			break;
	}
	EXIT_UNSUPPORTED;
	return -1;
}

static int compute_p_cub (const struct Op_IO* op_io, const int node_type, const struct Simulation* sim)
{
	const int cub_kind = op_io->kind;
// Add Simulation::p_X_p for each kind, X, for variable orders for a given reference order in future.
	switch (cub_kind) {
	case 's': // fallthrough
	case 'g':
	case 'm':
		return compute_p_basis(op_io,sim);
		break;
	case 'c':
UNUSED(node_type);
		// Note that cubature order depends on the node_type (May be p_op or p_c_x*p_op+p_c_p).
		// 1) Make function to judge which is the case.
		// 2) Set based on the node_type
		EXIT_ADD_SUPPORT;
		break;
	default:
		EXIT_ERROR("Unsupported: %c\n",cub_kind);
		break;
	}
	EXIT_UNSUPPORTED;
	return -1;
}

const struct const_Matrix_d* constructor_rst_ve
	(const int s_type, const int d_i, const int d_io, const int ind_h, const char ce, const struct Simulation* sim)
{
	const int e_type = compute_elem_from_super_type(s_type,d_io);
	if (e_type == PYR && ce != 'v')
		EXIT_ADD_SUPPORT;

	const struct const_Matrix_d* b_coords       = constructor_b_coords(e_type); // destructed
	const struct const_Vector_i* ind_h_b_coords = constructor_ind_h_b_coords(e_type,ce,ind_h,sim); // destructed

	const struct const_Matrix_d* vv0_vvs_vvs =
		constructor_subset_const_Matrix_d(b_coords,ind_h_b_coords); // destructed

	cubature_fptr constructor_Cubature = get_cubature_by_super_type(s_type);
// May need to use d_io below for fv, ev
	const struct const_Cubature* cub_ve = constructor_Cubature(d_io,1,CUB_VERTEX); // destructed

	const struct const_Matrix_d* rst_ve  =
		constructor_mm_const_Matrix_d('N','N',1.0,0.0,vv0_vvs_vvs,cub_ve->rst,'C'); // returned
	const_cast_ptrdiff(&rst_ve->ext_1,d_i);

	destructor_const_Matrix_d(b_coords);
	destructor_const_Vector_i(ind_h_b_coords);
	destructor_const_Matrix_d(vv0_vvs_vvs);
	destructor_const_Cubature(cub_ve);

	return rst_ve;
}

struct Matrix_d* constructor_rst_proj
	(const int e_type, const int d_i, const int d_io, const int ind_h, const char ce,
	 const struct const_Matrix_d* b_coords, const struct Simulation* sim)
{
EXIT_ERROR("Ensure that all is working as expected.\n");
	struct Matrix_d* rst_proj = NULL;

	const int s_type = compute_super_from_elem_type(e_type);
	if ((strcmp(sim->geom_blending[s_type],"gordon_hall") == 0) ||
	    (strcmp(sim->geom_blending[s_type],"szabo_babuska_gen") == 0))
	{
		const struct const_Matrix_d* rst_ve_proj = constructor_rst_ve(s_type,d_i,d_io,ind_h,ce,sim); // destructed
		rst_proj = constructor_mm_Matrix_d('N','N',1.0,0.0,b_coords,rst_ve_proj,'C'); // returned
		destructor_const_Matrix_d(rst_ve_proj);
	} else if (strcmp(sim->geom_blending[s_type],"scott") == 0) {
		EXIT_ADD_SUPPORT;
	} else if ((strcmp(sim->geom_blending[s_type],"lenoir") == 0) ||
	           (strcmp(sim->geom_blending[s_type],"nielson") == 0)) {
		EXIT_ADD_SUPPORT;
	} else {
		EXIT_ERROR("Unsupported: (%d, %s).\n",s_type,sim->geom_blending[s_type]);
	}

	return rst_proj;
}

// Level 1 ********************************************************************************************************** //

static int compute_node_type_cub_std
	(const struct Op_IO* op_io, const struct const_Element* element, const struct Simulation* sim)
{
	const char cub_kind = op_io->kind,
	           cub_ce   = op_io->ce;
	const int s_type = element->s_type;
	switch (cub_kind) {
	case 's':
		assert(cub_ce == 'v'); // Can be updated to include 'f' in future.
		if (strcmp(sim->nodes_interp[s_type],"GL") == 0)
			return CUB_GL;
		else if (strcmp(sim->nodes_interp[s_type],"GLL") == 0)
			return CUB_GLL;
		else if (strcmp(sim->nodes_interp[s_type],"AO") == 0)
			return CUB_AO;
		else if (strcmp(sim->nodes_interp[s_type],"WSH") == 0)
			return CUB_WSH;
		else if (strcmp(sim->nodes_interp[s_type],"EQ") == 0)
			return CUB_EQ;
		else
			EXIT_ERROR("Unsupported: %s\n",sim->nodes_interp[s_type]);
	case 'g':
		assert(cub_ce == 'v');
		switch (s_type) {
			case ST_TP:  return CUB_GLL; break;
			case ST_SI:  return CUB_AO;  break;
			case ST_PYR: return CUB_GLL; break;
			default:     EXIT_ERROR("Unsupported: %d\n",s_type); break;
		}
		break;
	case 'm':
		assert(cub_ce == 'v');
		EXIT_ADD_SUPPORT;
		break;
	case 'c':
		switch (s_type) {
		case ST_TP:
			return CUB_GL;
			break;
		case ST_SI:
			switch (sim->d) {
			case 2:
				switch (cub_ce) {
					case 'v': return CUB_WV; break;
					case 'f': return CUB_GL; break;
					case 'e': // fallthrough
					default:  EXIT_ERROR("Unsupported: %c\n",cub_ce); break;
				}
				break;
			case 3:
				switch (cub_ce) {
					case 'v': return CUB_WV; break;
					case 'f': return CUB_WV; break;
					case 'e': return CUB_GL; break;
					default:  EXIT_ERROR("Unsupported: %c\n",cub_ce); break;
				}
				break;
			default:
				EXIT_ERROR("Unsupported: %d\n",sim->d); break;
				break;
			}
			break;
		case ST_PYR: {
			const int out_e_type = compute_elem_type_sub_ce(element->type,cub_ce,op_io->h_op);

			switch (out_e_type) {
				case TRI:  return CUB_WV;  break;
				case QUAD: return CUB_GL;  break;
				case TET:  return CUB_WV;  break;
				case PYR:  return CUB_GJW; break;
				default:   EXIT_ERROR("Unsupported: %d\n",out_e_type);
			}
			break;
		} default:
			EXIT_ERROR("Unsupported: %d\n",s_type); break;
			break;
		}
		break;
	default:
		EXIT_ERROR("Unsupported: %c\n",cub_kind); break;
		break;
	}
	EXIT_UNSUPPORTED;
	return -1;
}

static int compute_node_type_cub_collocated
	(const struct Op_IO* op_io, const struct const_Element* element, const struct Simulation* sim)
{
	const char cub_kind = op_io->kind,
	           cub_ce   = op_io->ce;
	const int s_type = element->s_type;
	switch (cub_kind) {
	case 'c':
		switch (s_type) {
		case ST_TP:
			if (strcmp(sim->nodes_interp[s_type],"GL") == 0)
				return CUB_GL;
			else if (strcmp(sim->nodes_interp[s_type],"GLL") == 0)
				return CUB_GLL;
			else
				EXIT_ERROR("Unsupported: %s\n",sim->nodes_interp[s_type]);
			break;
		case ST_SI:
			if (strcmp(sim->nodes_interp[s_type],"WSH") != 0)
				EXIT_ERROR("Unsupported: %s\n",sim->nodes_interp[s_type]);

			switch (sim->d) {
			case 2:
				switch (cub_ce) {
					case 'v': return CUB_WSH; break;
					case 'f': return CUB_GL;  break;
					case 'e': // fallthrough
					default:  EXIT_ERROR("Unsupported: %c\n",cub_ce); break;
				}
				break;
			case 3:
				// If wedges are present, WSH nodes must be used for the face cubature for consistency with
				// tet faces. If not present, the higher strength WV nodes are used.
				if (wedges_present(sim->elements))
					return CUB_WSH;
				else
					return CUB_WV;
				break;
			default:
				break;
			}
			break;
		case ST_PYR: // fallthrough
		default:
			EXIT_ERROR("Unsupported: %d\n",s_type); break;
			break;
		}
		break;
	case 's': // fallthrough
	case 'g': // fallthrough
	case 'm': // fallthrough
	default:
		EXIT_ERROR("Unsupported: %c\n",cub_kind); break;
		break;
	}
	EXIT_UNSUPPORTED;
	return -1;
}

const struct const_Matrix_d* constructor_b_coords (const int e_type)
{
	/* `ext_0` is not necessarily equal to the number of p2 rst coordinate nodes for the `e_type`. Additional nodes
	 * may be included for certain refinements (such as splitting a TET into 12 TETs). */
	ptrdiff_t ext_0 = -1,
	          ext_1 = -1; // Number of vertices of the `e_type`.
	const double* b_coords = NULL;

	switch (e_type) {
	case LINE:
		ext_0 = 3;
		ext_1 = 2;
		b_coords = (double[])
			{ 1.0 , 0.0 ,
			  0.5 , 0.5 ,
			  0.0 , 1.0 ,};
		break;
	case TRI:
		ext_0 = 6;
		ext_1 = 3;
		b_coords = (double[])
			{ 1.0 , 0.0 , 0.0 ,
			  0.5 , 0.5 , 0.0 ,
			  0.0 , 1.0 , 0.0 ,
			  0.5 , 0.0 , 0.5 ,
			  0.0 , 0.5 , 0.5 ,
			  0.0 , 0.0 , 1.0 ,};
		break;
	case TET:
		ext_0 = 11;
		ext_1 = 4;
		b_coords = (double[])
			{ 1.0 , 0.0 , 0.0 , 0.0 ,
			  0.5 , 0.5 , 0.0 , 0.0 ,
			  0.0 , 1.0 , 0.0 , 0.0 ,
			  0.5 , 0.0 , 0.5 , 0.0 ,
			  0.0 , 0.5 , 0.5 , 0.0 ,
			  0.0 , 0.0 , 1.0 , 0.0 ,
			  0.5 , 0.0 , 0.0 , 0.5 ,
			  0.0 , 0.5 , 0.0 , 0.5 ,
			  0.0 , 0.0 , 0.5 , 0.5 ,
			  0.0 , 0.0 , 0.0 , 1.0 ,
			  0.25, 0.25, 0.25, 0.25,};
		break;
	case PYR:
		ext_0 = 14;
		ext_1 = 5;
		b_coords = (double[])
			{ 1.0 , 0.0 , 0.0 , 0.0 , 0.0 ,
			  0.5 , 0.5 , 0.0 , 0.0 , 0.0 ,
			  0.0 , 1.0 , 0.0 , 0.0 , 0.0 ,
			  0.5 , 0.0 , 0.5 , 0.0 , 0.0 ,
			  0.25, 0.25, 0.25, 0.25, 0.0 ,
			  0.0 , 0.5 , 0.0 , 0.5 , 0.0 ,
			  0.0 , 0.0 , 1.0 , 0.0 , 0.0 ,
			  0.0 , 0.0 , 0.5 , 0.5 , 0.0 ,
			  0.0 , 0.0 , 0.0 , 1.0 , 0.0 ,
			  0.5 , 0.0 , 0.0 , 0.0 , 0.5 ,
			  0.0 , 0.5 , 0.0 , 0.0 , 0.5 ,
			  0.0 , 0.0 , 0.5 , 0.0 , 0.5 ,
			  0.0 , 0.0 , 0.0 , 0.5 , 0.5 ,
			  0.0 , 0.0 , 0.0 , 0.0 , 1.0 ,};
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n.",e_type);
		break;
	}

	return constructor_copy_const_Matrix_d_d('R',ext_0,ext_1,b_coords);
}

const struct const_Vector_i* constructor_ind_h_b_coords
	(const int e_type, const char ce, const int ind_h, const struct Simulation* sim)
{
	const struct const_Element* element = get_element_by_type(sim->elements,e_type);

	ptrdiff_t ext_0 = -1;
	const int* ind_h_b_coords = NULL;

	switch (e_type) {
	case LINE:
		switch (ce) {
		case 'v': {
			enum { n_ve_ce = 2, n_ref_max = 3, };
			assert(n_ve_ce == element->n_ve);
			assert(n_ref_max == element->n_ref_max_v);
			assert(ind_h < n_ref_max);

			static const int ind_h_b_coords_all[n_ref_max][n_ve_ce] =
				{{ 0, 2,},

				 { 0, 1,},
				 { 1, 2,},
				};

			ext_0 = n_ve_ce;
			ind_h_b_coords = ind_h_b_coords_all[ind_h];
			break;
		} case 'f': { // fallthrough
		} case 'e': {
			enum { n_ve_ce = 1, n_ref_max = 2*1, };
			assert(n_ref_max == element->n_ref_max_f);
			assert(ind_h < n_ref_max);

			static const int ind_h_b_coords_all[n_ref_max][n_ve_ce] =
				{{ 0, },
				 { 2, },
				};

			ext_0 = n_ve_ce;
			ind_h_b_coords = ind_h_b_coords_all[ind_h];
			break;
		} default:
			EXIT_ERROR("Unsupported: %c.\n",ce);
			break;
		}
		break;
	case TRI:
		switch (ce) {
		case 'v': {
			enum { n_ve_ce = 3, n_ref_max = 5, };
			assert(n_ve_ce == element->n_ve);
			assert(n_ref_max == element->n_ref_max_v);
			assert(ind_h < n_ref_max);

			static const int ind_h_b_coords_all[n_ref_max][n_ve_ce] =
				{{ 0, 2, 5,},

				 { 0, 1, 3,},
				 { 1, 2, 4,},
				 { 3, 4, 5,},
				 { 4, 3, 1,},
				};

			ext_0 = n_ve_ce;
			ind_h_b_coords = ind_h_b_coords_all[ind_h];
			break;
		} case 'f': { // fallthrough
		} case 'e': {
			enum { n_ve_ce = 2, n_ref_max = 3*3, };
			assert(n_ref_max == element->n_ref_max_f);
			assert(ind_h < n_ref_max);

			static const int ind_h_b_coords_all[n_ref_max][n_ve_ce] =
				{{ 2, 5,},
				 { 0, 5,},
				 { 0, 2,},

				 { 2, 4,}, { 4, 5,},
				 { 0, 3,}, { 3, 5,},
				 { 0, 1,}, { 1, 2,},
				};

			ext_0 = n_ve_ce;
			ind_h_b_coords = ind_h_b_coords_all[ind_h];
			break;
		} default:
			EXIT_ERROR("Unsupported: %c.\n",ce);
			break;
		}
		break;
	case TET:
		switch (ce) {
		case 'v': {
			enum { n_ve_ce = 4, n_ref_max = 9, };
			assert(n_ve_ce == element->n_ve);
			assert(n_ref_max == element->n_ref_max_v);
			assert(ind_h < n_ref_max);

			static const int ind_h_b_coords_all[n_ref_max][n_ve_ce] =
				{{  0,  2,  5,  9,},

				 {  0,  1,  3,  6,},
				 {  1,  2,  4,  7,},
				 {  3,  4,  5,  8,},
				 {  6,  7,  8,  9,},
				 {  1,  8,  7,  4,},
				 {  8,  1,  6,  3,},
				 {  7,  6,  8,  1,},
				 {  4,  3,  1,  8,},
				};

			ext_0 = n_ve_ce;
			ind_h_b_coords = ind_h_b_coords_all[ind_h];
			break;
		} case 'f': {
			enum { n_ve_ce = 3, n_ref_max = 4*5, };
			assert(n_ref_max == element->n_ref_max_f);
			assert(ind_h < n_ref_max);

			static const int ind_h_b_coords_all[n_ref_max][n_ve_ce] =
				{{  2,  5,  9,},
				 {  0,  5,  9,},
				 {  0,  2,  9,},
				 {  0,  2,  5,},

				 {  2,  4,  7,}, {  4,  5,  8,}, {  7,  8,  9,}, {  8,  7,  4,},
				 {  0,  3,  6,}, {  3,  5,  8,}, {  6,  8,  9,}, {  8,  6,  3,},
				 {  0,  1,  6,}, {  1,  2,  7,}, {  6,  7,  9,}, {  7,  6,  1,},
				 {  0,  1,  3,}, {  1,  2,  4,}, {  3,  4,  5,}, {  4,  3,  1,},
				};

			ext_0 = n_ve_ce;
			ind_h_b_coords = ind_h_b_coords_all[ind_h];
		} case 'e': {
			enum { n_ve_ce = 2, n_ref_max = 6*3, };
			assert(n_ref_max == element->n_ref_max_e);
			assert(ind_h < n_ref_max);

			static const int ind_h_b_coords_all[n_ref_max][n_ve_ce] =
				{{  2,  5,},
				 {  0,  5,},
				 {  0,  2,},
				 {  0,  9,},
				 {  2,  9,},
				 {  5,  9,},

				 {  2,  4,}, {  4,  5,},
				 {  0,  3,}, {  3,  5,},
				 {  0,  1,}, {  1,  2,},
				 {  0,  6,}, {  6,  9,},
				 {  2,  7,}, {  7,  9,},
				 {  5,  8,}, {  8,  9,},
				};

			ext_0 = n_ve_ce;
			ind_h_b_coords = ind_h_b_coords_all[ind_h];
			break;
		} default:
			EXIT_ERROR("Unsupported: %c.\n",ce);
			break;
		}
	case PYR:
		switch (ce) {
		case 'v': {
			enum { n_ve_ce_max = 5, n_ref_max = 11, };
			assert(n_ve_ce_max == element->n_ve);
			assert(n_ref_max == element->n_ref_max_v);
			assert(ind_h < n_ref_max);

			static const int ind_h_b_coords_all[n_ref_max][n_ve_ce_max] =
				{{  0,  2,  6,  8, 13,},

				 {  0,  1,  3,  4,  9,},
				 {  1,  2,  4,  5, 10,},
				 {  3,  4,  6,  7, 11,},
				 {  4,  5,  7,  8, 12,},
				 {  3,  4, 11,  9, -1,},
				 {  4,  5, 12, 10, -1,},
				 { 10,  9,  4,  1, -1,},
				 { 12, 11,  7,  4, -1,},
				 { 10,  9, 12, 11,  4,},
				 {  9, 10, 11, 12, 13,},
				};
			static const int n_ve_ce_all[n_ref_max] = {5, 5,5,5,5,4,4,4,4,5,5, };

			ext_0 = n_ve_ce_all[ind_h];
			ind_h_b_coords = ind_h_b_coords_all[ind_h];
			break;
		} case 'f': {
			enum { n_ve_ce_max = 4, n_ref_max = 5*5, };
			assert(n_ref_max == element->n_ref_max_f);
			assert(ind_h < n_ref_max);

			static const int ind_h_b_coords_all[n_ref_max][n_ve_ce_max] =
				{{  0,  6, 13, -1,},
				 {  2,  8, 13, -1,},
				 {  0,  2, 13, -1,},
				 {  6,  8, 13, -1,},
				 {  0,  2,  6,  8, },

				 {  0,  3,  9, -1,}, {  3,  6, 11, -1,}, {  9, 11, 13, -1,}, { 11,  9,  3, -1,},
				 {  2,  5, 10, -1,}, {  5,  8, 12, -1,}, { 10, 12, 13, -1,}, { 12, 10,  5, -1,},
				 {  0,  1,  9, -1,}, {  1,  2, 10, -1,}, {  9, 10, 13, -1,}, { 10,  9,  1, -1,},
				 {  6,  7, 11, -1,}, {  7,  8, 12, -1,}, { 11, 12, 13, -1,}, { 12, 11,  7, -1,},
				 {  0,  1,  3,  4,}, {  1,  2,  4,  5,}, {  3,  4,  6,  7,}, {  4,  5,  7,  8,},
				};
			static const int n_ve_ce_all[n_ref_max] = {3,3,3,3,4, 3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,};

			ext_0 = n_ve_ce_all[ind_h];
			ind_h_b_coords = ind_h_b_coords_all[ind_h];
		} case 'e': {
			enum { n_ve_ce = 2, n_ref_max = 8*3, };
			assert(n_ref_max == element->n_ref_max_e);
			assert(ind_h < n_ref_max);

			static const int ind_h_b_coords_all[n_ref_max][n_ve_ce] =
				{{  0,  6,},
				 {  2,  8,},
				 {  0,  2,},
				 {  6,  8,},
				 {  0, 13,},
				 {  2, 13,},
				 {  6, 13,},
				 {  8, 13,},

				 {  0,  3,}, {  3,  6,},
				 {  2,  5,}, {  5,  8,},
				 {  0,  1,}, {  1,  2,},
				 {  6,  7,}, {  7,  8,},
				 {  0,  9,}, {  9, 13,},
				 {  0, 10,}, { 10, 13,},
				 {  0, 11,}, { 11, 13,},
				 {  0, 12,}, { 12, 13,},
				};


			ext_0 = n_ve_ce;
			ind_h_b_coords = ind_h_b_coords_all[ind_h];
			break;
		} default:
			EXIT_ERROR("Unsupported: %c.\n",ce);
			break;
		}
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n.",e_type);
		break;
	}

	return constructor_copy_const_Vector_i_i(ext_0,ind_h_b_coords);
}
