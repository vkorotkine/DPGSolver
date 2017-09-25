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

#include "matrix.h"

#include "simulation.h"
#include "element.h"
#include "bases.h"
#include "cubature.h"

// Static function declarations ************************************************************************************* //

/** \brief Compute the super type of the cubature based on the kind of operator.
 *  \return See brief. */
static int compute_super_type_cub
	(const struct Op_IO* op_io,          ///< \ref Op_IO.
	 const struct const_Element* element ///< \ref const_Element.
	);

/** \brief Compute the node type of the cubature based on the kind of operator.
 *  \return See brief. */
static int compute_node_type_cub
	(const struct Op_IO* op_io,           ///< \ref Op_IO.
	 const struct const_Element* element, ///< \ref const_Element.
	 const struct Simulation *sim         ///< \ref Simulation.
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
	 const int s_type,            ///< \ref Element::s_type.
	 const int node_type,         ///< \ref Cubature::node_type.
	 const struct Simulation *sim ///< \ref Simulation.
	);

/** \brief Constructor for the rst coordinates associated with vertices of the (potentially h-refined) reference
 *         element.
 *  \return See brief. */
const struct const_Matrix_d* constructor_rst_ve
	(const int s_type, ///< \ref Element::s_type.
	 const int d,      ///< The dimension of the rst coordinates.
	 const int ind_h,  ///< The h-refinement index.
	 const char ce     ///< \ref Op_IO::ce.
	);

// Constructor functions ******************************************************************************************** //

// Possibly return a Cubature_Multiarray here.
const struct const_Cubature* constructor_const_Cubature_h
	(const int ind_io, const struct Op_IO op_io[2], const struct const_Element* element,
	 const struct Simulation* sim)
{
	const int s_type_io    = compute_super_type_cub(&op_io[ind_io],element);
	const int node_type_io = compute_node_type_cub(&op_io[ind_io],element,sim);
	const int d_io         = compute_d_cub(op_io[ind_io].ce,element->d);
	const int p_io         = compute_p_cub(&op_io[ind_io],s_type_io,node_type_io,sim);
	const int h_io         = op_io[ind_io].h_op;

	cubature_fptr constructor_Cubature = get_cubature_by_super_type(s_type_io);
	basis_fptr    constructor_basis    = get_basis_by_super_type(s_type_io,"ortho");

	const struct const_Cubature* cub_io = constructor_Cubature(d_io,p_io,node_type_io);

	const struct const_Matrix_d* rst_ve = constructor_rst_ve(s_type_io,d_io,0,op_io[ind_io].ce); // tbd

	// vXX: (v)olume XX (Arbitrary, would be replated with [op_io.kind,op_io.sc] if specified)
	const struct const_Matrix_d* phi_vvs_vXX     = constructor_basis(1,cub_io->rst); // tbd
	const struct const_Matrix_d* phi_vvs_vvs     = constructor_basis(1,rst_ve); // tbd
	const struct const_Matrix_d* phi_vvs_vvs_inv = constructor_inverse_const_Matrix_d(phi_vvs_vvs); // tbd
	const struct const_Matrix_d* cv0_vvs_vXX =
		constructor_mm_const_Matrix_d('N','N',1.0,0.0,phi_vvs_vXX,phi_vvs_vvs_inv,'R'); // tbd

	const struct const_Matrix_d* rst_ve_io = constructor_rst_ve(s_type_io,d_io,h_io,op_io[ind_io].ce); // tbd

	struct Cubature* cubature = malloc(sizeof* cubature); // returned

	cubature->p           = p_io;
	cubature->node_type   = node_type_io;
	cubature->rst         = constructor_mm_Matrix_d('N','N',1.0,0.0,cv0_vvs_vXX,rst_ve_io,'C'); // moved
	cubature->has_weights = false;
	cubature->w           = NULL;

	return (const struct const_Cubature*) cubature;


//	const struct const_Cubature* cub_ve = constructor_Cubature(d_i,1,CUB_VERTEX);
// Need function to return appropriate subset of cub_ve.



}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Compute the node type of the cubature based on the kind of operator for standard (non-collocated) nodes.
 *  \return See brief. */
static int compute_node_type_cub_std
	(const struct Op_IO* op_io,           ///< \ref Op_IO.
	 const struct const_Element* element, ///< \ref const_Element.
	 const struct Simulation *sim         ///< \ref Simulation.
	);

/** \brief Compute the node type of the cubature based on the kind of operator for collocated nodes.
 *  \return See brief. */
static int compute_node_type_cub_collocated
	(const struct Op_IO* op_io,           ///< \ref Op_IO.
	 const struct const_Element* element, ///< \ref const_Element.
	 const struct Simulation *sim         ///< \ref Simulation.
	);

static int compute_super_type_cub (const struct Op_IO* op_io, const struct const_Element* element)
{
	const int sub_e_type = compute_elem_type_sub_ce(element->type,op_io->ce,op_io->h_op);
	return compute_super_from_elem_type(sub_e_type);
}

static int compute_node_type_cub
	(const struct Op_IO* op_io, const struct const_Element* element, const struct Simulation *sim)
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

static int compute_p_cub
	(const struct Op_IO* op_io, const int s_type, const int node_type, const struct Simulation *sim)
{
	const int p_ref  = op_io->p_op,
	          cub_kind = op_io->kind,
	          cub_sc   = op_io->sc;
// Add Simulation::p_X_p for each kind, X, for variable orders for a given reference order in future.
	switch (cub_kind) {
	case 's': // solution
		return p_ref;
		break;
	case 'g':
		if (cub_sc == 's')
			return 1;

		if (strcmp(sim->geom_rep,"isoparametric") == 0)
			return p_ref;
		else if (strcmp(sim->geom_rep,"superparametric") == 0)
			return p_ref+1;
		else if (strstr(sim->geom_rep,"fixed"))
			EXIT_ADD_SUPPORT; // Find number in geom_rep (use something similar to 'convert_to_range_d').
		else
			EXIT_ERROR("Unsupported: %s\n",sim->geom_rep);
		break;
	case 'm': {
		EXIT_ADD_SUPPORT;
UNUSED(s_type);

		break;
	} case 'c':
UNUSED(node_type);
		// Note that cubature order depends on the node_type (May be p_ref or p_c_x*p_ref+p_c_p).
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

const struct const_Matrix_d* constructor_rst_ve (const int s_type, const int d, const int ind_h, const char ce)
{
	if (ce != 'v')
		EXIT_ADD_SUPPORT;

	cubature_fptr constructor_Cubature = get_cubature_by_super_type(s_type);
	const struct const_Cubature* cub_ve = constructor_Cubature(d,1,CUB_VERTEX); // destructed

	const ptrdiff_t ext_1 = cub_ve->rst->ext_0;
	ptrdiff_t ext_0 = -1;
	double* sub_rst = NULL;

	const int e_type = compute_elem_from_super_type(s_type,d);
/// \todo Move this to be with the h-refinement related connectivity information?
	switch (e_type) {
	case LINE: {
		static const double barycentric_coords[3][2] = {{ 1.0 , 0.0 },
		                                                { 0.0 , 1.0 },
		                                                { 0.5 , 0.5 },};
		static const int ind_h_coords_all[3][2] = {{ 0 , 1 },
		                                           { 0 , 2 },
		                                           { 2 , 1 },};

		const int*const ind_h_coords = ind_h_coords_all[ind_h];
		ext_0 = ext_1;

		sub_rst = malloc(ext_1*ext_0 * sizeof *sub_rst); // moved
		double* sub_rst_ptr = sub_rst;

		for (ptrdiff_t i = 0; i < ext_0; ++i) {
		for (ptrdiff_t j = 0; j < ext_1; ++j) {
			*sub_rst_ptr++ = barycentric_coords[ind_h_coords[i]][j];
		}}
		break;
	} case TRI: {
		static const double barycentric_coords[6][3] = {{ 1.0 , 0.0 , 0.0 },
		                                                { 0.0 , 1.0 , 0.0 },
		                                                { 0.0 , 0.0 , 1.0 },
		                                                { 0.5 , 0.5 , 0.0 },
		                                                { 0.5 , 0.0 , 0.5 },
		                                                { 0.0 , 0.5 , 0.5 },};
		static const int ind_h_coords_all[5][3] = {{ 0 , 1 , 2 },
		                                           { 0 , 3 , 4 },
		                                           { 3 , 1 , 5 },
		                                           { 4 , 5 , 2 },
		                                           { 5 , 4 , 3 },};

		const int*const ind_h_coords = ind_h_coords_all[ind_h];
		ext_0 = ext_1;

		sub_rst = malloc(ext_1*ext_0 * sizeof *sub_rst); // moved
		double* sub_rst_ptr = sub_rst;

		for (ptrdiff_t i = 0; i < ext_0; ++i) {
		for (ptrdiff_t j = 0; j < ext_1; ++j) {
			*sub_rst_ptr++ = barycentric_coords[ind_h_coords[i]][j];
		}}
		break;
	} case QUAD: { // For PYR faces
		static const double barycentric_coords[9][4] = {{ 1.0 , 0.0 , 0.0 , 0.0 },
		                                                { 0.0 , 1.0 , 0.0 , 0.0 },
		                                                { 0.0 , 0.0 , 1.0 , 0.0 },
		                                                { 0.0 , 0.0 , 0.0 , 1.0 },
		                                                { 0.5 , 0.5 , 0.0 , 0.0 },
		                                                { 0.5 , 0.0 , 0.5 , 0.0 },
		                                                { 0.0 , 0.5 , 0.0 , 0.5 },
		                                                { 0.0 , 0.0 , 0.5 , 0.5 },
		                                                { 0.25, 0.25, 0.25, 0.25},};
		static const int ind_h_coords_all[9][4] = {{ 0 , 1 , 2 , 3 },
		                                           { 0 , 4 , 5 , 8 },
		                                           { 4 , 1 , 8 , 6 },
		                                           { 5 , 8 , 2 , 7 },
		                                           { 8 , 6 , 7 , 3 },
		                                           { 0 , 4 , 2 , 7 },
		                                           { 4 , 1 , 7 , 3 },
		                                           { 0 , 1 , 5 , 6 },
		                                           { 5 , 6 , 2 , 3 },};

		const int*const ind_h_coords = ind_h_coords_all[ind_h];
		ext_0 = ext_1;

		sub_rst = malloc(ext_1*ext_0 * sizeof *sub_rst); // moved
		double* sub_rst_ptr = sub_rst;

		for (ptrdiff_t i = 0; i < ext_0; ++i) {
		for (ptrdiff_t j = 0; j < ext_1; ++j) {
			*sub_rst_ptr++ = barycentric_coords[ind_h_coords[i]][j];
		}}
		break;
	} case TET: {
		static const double barycentric_coords[11][4] = {{ 1.0 , 0.0 , 0.0 , 0.0 },
		                                                 { 0.0 , 1.0 , 0.0 , 0.0 },
		                                                 { 0.0 , 0.0 , 1.0 , 0.0 },
		                                                 { 0.0 , 0.0 , 0.0 , 1.0 },
		                                                 { 0.5 , 0.5 , 0.0 , 0.0 },
		                                                 { 0.5 , 0.0 , 0.5 , 0.0 },
		                                                 { 0.0 , 0.5 , 0.5 , 0.0 },
		                                                 { 0.5 , 0.0 , 0.0 , 0.5 },
		                                                 { 0.0 , 0.5 , 0.0 , 0.5 },
		                                                 { 0.0 , 0.0 , 0.5 , 0.5 },
		                                                 { 0.25, 0.25, 0.25, 0.25},};
		static const int ind_h_coords_all[9][4] = {{ 0 , 1 , 2 , 3 },
		                                           { 0 , 4 , 5 , 7 },
		                                           { 4 , 1 , 6 , 8 },
		                                           { 5 , 6 , 2 , 9 },
		                                           { 7 , 8 , 9 , 3 },
		                                           { 4 , 9 , 8 , 6 },
		                                           { 9 , 4 , 7 , 5 },
		                                           { 8 , 7 , 9 , 4 },
		                                           { 6 , 5 , 4 , 9 },};

		const int*const ind_h_coords = ind_h_coords_all[ind_h];
		ext_0 = ext_1;

		sub_rst = malloc(ext_1*ext_0 * sizeof *sub_rst); // moved
		double* sub_rst_ptr = sub_rst;

		for (ptrdiff_t i = 0; i < ext_0; ++i) {
		for (ptrdiff_t j = 0; j < ext_1; ++j) {
			*sub_rst_ptr++ = barycentric_coords[ind_h_coords[i]][j];
		}}
		break;
	} case PYR: {
		static const double barycentric_coords[14][5] = {{ 1.0 , 0.0 , 0.0 , 0.0 , 0.0 },
		                                                 { 0.0 , 1.0 , 0.0 , 0.0 , 0.0 },
		                                                 { 0.0 , 0.0 , 1.0 , 0.0 , 0.0 },
		                                                 { 0.0 , 0.0 , 0.0 , 1.0 , 0.0 },
		                                                 { 0.0 , 0.0 , 0.0 , 0.0 , 1.0 },
		                                                 { 0.5 , 0.5 , 0.0 , 0.0 , 0.0 },
		                                                 { 0.5 , 0.0 , 0.5 , 0.0 , 0.0 },
		                                                 { 0.25, 0.25, 0.25, 0.25, 0.0 },
		                                                 { 0.0 , 0.5 , 0.0 , 0.5 , 0.0 },
		                                                 { 0.0 , 0.0 , 0.5 , 0.5 , 0.0 },
		                                                 { 0.5 , 0.0 , 0.0 , 0.0 , 0.5 },
		                                                 { 0.0 , 0.5 , 0.0 , 0.0 , 0.5 },
		                                                 { 0.0 , 0.0 , 0.5 , 0.0 , 0.5 },
		                                                 { 0.0 , 0.0 , 0.0 , 0.5 , 0.5 },};
		static const int ind_h_coords_all[11][5] = {{ 0 , 1 , 2 , 3 , 4 },
		                                            { 0 , 5 , 6 , 7 , 10},
		                                            { 5 , 1 , 7 , 8 , 11},
		                                            { 6 , 7 , 2 , 9 , 12},
		                                            { 7 , 8 , 9 , 3 , 13},
		                                            { 6 , 7 , 12, 10, -1},
		                                            { 7 , 8 , 13, 11, -1},
		                                            { 11, 10, 7 , 5 , -1},
		                                            { 13, 12, 9 , 7 , -1},
		                                            { 11, 10, 13, 12, 7 },
		                                            { 10, 11, 12, 13, 4 },};
		static const ptrdiff_t n_ve_o_all[11] = {5,5,5,5,5,4,4,4,4,5,5,};

		const int*const ind_h_coords = ind_h_coords_all[ind_h];
		ext_0 = n_ve_o_all[ind_h];

		sub_rst = malloc(ext_1*ext_0 * sizeof *sub_rst); // moved
		double* sub_rst_ptr = sub_rst;

		for (ptrdiff_t i = 0; i < ext_0; ++i) {
		for (ptrdiff_t j = 0; j < ext_1; ++j) {
			*sub_rst_ptr++ = barycentric_coords[ind_h_coords[i]][j];
		}}
		break;
	} default:
		EXIT_ERROR("Unsupported: %d\n",e_type);
		break;
	}

	const struct const_Matrix_d* sub_rst_M =
		constructor_move_const_Matrix_d_d('R',ext_0,ext_1,true,sub_rst); // destructed
	const struct const_Matrix_d* rst_ve  =
		constructor_mm_const_Matrix_d('N','N',1.0,0.0,sub_rst_M,cub_ve->rst,'C'); // returned

	destructor_const_Matrix_d(sub_rst_M);
	destructor_const_Cubature(cub_ve);

	return rst_ve;
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
