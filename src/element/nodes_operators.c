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

#include "nodes_operators.h"
#include "element_operators.h"

#include <assert.h>
#include <string.h>
#include "gsl/gsl_math.h"

#include "macros.h"
#include "definitions_core.h"
#include "definitions_nodes.h"
#include "definitions_element_operators.h"

#include "matrix.h"
#include "multiarray.h"
#include "vector.h"

#include "bases.h"
#include "const_cast.h"
#include "element.h"
#include "element_operators.h"
#include "nodes.h"
#include "nodes_correspondence.h"
#include "simulation.h"

// Static function declarations ************************************************************************************* //

/** \brief Compute the node type based on the kind of operator.
 *  \return See brief. */
static int compute_node_type
	(const struct Op_IO* op_io,           ///< \ref Op_IO.
	 const struct const_Element* element, ///< \ref const_Element.
	 const struct Simulation* sim         ///< \ref Simulation.
	);

/** \brief Compute the dimension of the nodes based off of the computational element and the simulation dimension.
 *  \return See brief. */
static int compute_d_nodes
	(const char ce,   ///< \ref Op_IO::ce.
	 const int elem_d ///< \ref Element::d.
	);

/** \brief Compute the order of the nodes to be computed based on the reference order and the kind of operator.
 *  \return See brief. */
static int compute_p_nodes
	(const struct Op_IO* op_io,   ///< \ref Op_IO.
	 const int node_type,         ///< \ref Nodes::node_type.
	 const struct Simulation* sim ///< \ref Simulation.
	);

/** \brief Constructor for the rst coordinates associated with vertices of the (potentially h-refined) reference
 *         element.
 *  \return See brief. */
static const struct const_Matrix_d* constructor_rst_ve
	(const int s_type,            ///< \ref Element::s_type.
	 const int d_i,               ///< The dimension of the input basis coordinates.
	 const int d_io,              ///< The dimension of the input/output rst coordinates.
	 const int ind_h,             ///< The h-refinement index.
	 const int ind_ce,            ///< The computational element index.
	 const char ce,               ///< \ref Op_IO::ce.
	 const struct Simulation* sim ///< \ref Simulation.
	);

/** \brief Constructor for rst coordinates as a projection of
 *  \return See brief. */
static struct Matrix_d* constructor_rst_proj
	(const int e_type,                      ///< \ref Element::type.
	 const int d_i,                         ///< The dimension of the input basis coordinates.
	 const int d_io,                        ///< The dimension of the input/output rst coordinates.
	 const int ind_h,                       ///< The h-refinement index.
	 const int ind_ce,                      ///< The computational element index.
	 const char ce,                         ///< \ref Op_IO::ce.
	 const struct const_Matrix_d* b_coords, ///< The barycentric coordinates of the element from which to project.
	 const struct Simulation* sim           ///< \ref Simulation.
	);

// Constructor functions ******************************************************************************************** //

const struct const_Nodes* constructor_const_Nodes_h
	(const int ind_io, const struct Op_IO op_io[2], const struct const_Element* element,
	 const struct Simulation* sim)
{
	const int s_type_i     = op_io[OP_IND_I].s_type,
	          s_type_io    = op_io[ind_io].s_type,
	          node_type_io = compute_node_type(&op_io[ind_io],element,sim);

	// Always use the input to establish the dimension of the nodes as it is in bases of this dimension in which the
	// nodes will be used.
	const int d_i       = compute_d_nodes(op_io[OP_IND_I].ce,element->d),
	          d_io      = compute_d_nodes(op_io[ind_io].ce,element->d),
	          p_io      = compute_p_nodes(&op_io[ind_io],node_type_io,sim),
	          ind_h_io  = op_io[ind_io].h_op,
	          ind_ce_io = op_io[ind_io].ce_op;

	constructor_Nodes_fptr constructor_Nodes = get_constructor_Nodes_by_super_type(s_type_io);
	constructor_basis_fptr constructor_basis = get_constructor_basis_by_super_type(s_type_io,"orthonormal");

	const struct const_Nodes* nodes_io = constructor_Nodes(d_io,p_io,node_type_io); // destructed

	const struct const_Matrix_d* rst_ve = constructor_rst_ve(s_type_io,d_io,d_io,0,0,'v',sim); // destructed

	// vXX: (v)olume XX (Arbitrary, would be replaced with [op_io.kind,op_io.sc] if specified)
	const struct const_Matrix_d* cv0r_vvs_vvs     = constructor_basis(1,rst_ve);                      // destructed
	const struct const_Matrix_d* inv_cv0r_vvs_vvs = constructor_inverse_const_Matrix_d(cv0r_vvs_vvs); // destructed
	const struct const_Matrix_d* cv0r_vvs_vXX     = constructor_basis(1,nodes_io->rst);               // destructed
	const struct const_Matrix_d* cv0_vvs_vXX =
		constructor_mm_const_Matrix_d('N','N',1.0,cv0r_vvs_vXX,inv_cv0r_vvs_vvs,'R'); // destructed
	destructor_const_Matrix_d(rst_ve);
	destructor_const_Matrix_d(cv0r_vvs_vvs);
	destructor_const_Matrix_d(inv_cv0r_vvs_vvs);
	destructor_const_Matrix_d(cv0r_vvs_vXX);

	struct Nodes* nodes = malloc(sizeof* nodes); // returned

	nodes->has_weights = nodes_io->has_weights;
	nodes->w           = (struct Vector_d*) nodes_io->w;

	nodes->has_symms = nodes_io->has_symms;
	nodes->s         = (struct Vector_i*) nodes_io->s;

	const_cast_b(&nodes_io->has_weights,false);
	const_cast_b(&nodes_io->has_symms,false);
	destructor_const_Nodes(nodes_io);

	nodes->p           = p_io;
	nodes->node_type   = node_type_io;

	const char ce_i  = op_io[OP_IND_I].ce,
	           ce_o  = op_io[OP_IND_O].ce,
	           ce_io = op_io[ind_io].ce;
/// \todo Attempt to remove one of the (potentially redundant ce_x == ce_y) tests in the if condition(s) below.
	if ((ce_i == ce_o) || (ce_i == 'v') || (ce_i == ce_io)) { // ((vv || ff || ee) || (vf || ve))
		// Compute the output rst coordinates by multiplying the barycentric coordinates of the nodes with the
		// appropriate (sub)set of reference element vertices.
		const char ce_n = ((ce_i == ce_o) || (ce_i == ce_io) ? 'v' : op_io[ind_io].ce);
		const struct const_Matrix_d* rst_ve_io =
			constructor_rst_ve(s_type_i,d_i,d_io,ind_h_io,ind_ce_io,ce_n,sim); // destructed
		nodes->rst = constructor_mm_Matrix_d('N','N',1.0,cv0_vvs_vXX,rst_ve_io,'C'); // keep
		destructor_const_Matrix_d(rst_ve_io);
	} else if (ce_o == 'v') { // (fv || ev)
	      const int ind_ce_i = op_io[OP_IND_I].ce_op;
		nodes->rst = constructor_rst_proj(element->type,d_i,d_io,ind_h_io,ind_ce_i,ce_i,cv0_vvs_vXX,sim); // keep
	} else { // (fe || ef)
		EXIT_ERROR("Add support: %c %c\n",ce_i,ce_o);
	}
	destructor_const_Matrix_d(cv0_vvs_vXX);

	return (const struct const_Nodes*) nodes;
}

const struct const_Multiarray_Vector_i* constructor_nodes_face_corr_op
	(const struct Op_IO* op_io, const struct const_Element* element, const struct Simulation* sim)
{
	const int node_type = compute_node_type(op_io,element,sim),
	          s_type    = op_io->s_type,
	          d         = compute_d_nodes(op_io->ce,element->d),
	          p         = compute_p_nodes(op_io,node_type,sim);
	return constructor_nodes_face_corr(d,p,node_type,s_type);
}

const struct const_Vector_d* constructor_weights
	(const struct Op_IO* op_io, const struct const_Element* element, const struct Simulation* sim)
{
	const int s_type    = op_io->s_type,
	          node_type = compute_node_type(op_io,element,sim),
	          d         = compute_d_nodes(op_io->ce,element->d),
	          p         = compute_p_nodes(op_io,node_type,sim);

	constructor_Nodes_fptr constructor_Nodes = get_constructor_Nodes_by_super_type(s_type);

	const struct const_Nodes* nodes = constructor_Nodes(d,p,node_type); // destructed
	assert(nodes->has_weights);

	const struct const_Vector_d* w = nodes->w;

	const_cast_b(&nodes->has_weights,false);
	destructor_const_Nodes(nodes);

	return w;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Compute the node type based on the kind of operator for standard (non-collocated) nodes.
 *  \return See brief. */
static int compute_node_type_std
	(const struct Op_IO* op_io,           ///< \ref Op_IO.
	 const struct const_Element* element, ///< \ref const_Element.
	 const struct Simulation* sim         ///< \ref Simulation.
	);

/** \brief Compute the node type based on the kind of operator for collocated nodes.
 *  \return See brief. */
static int compute_node_type_collocated
	(const struct Op_IO* op_io,           ///< \ref Op_IO.
	 const struct const_Element* element, ///< \ref const_Element.
	 const struct Simulation* sim         ///< \ref Simulation.
	);

/** \brief Compute the index of the type of processing required to get from the reference order to the cubature
 *         (integration) order.
 *  \return See brief. */
static int compute_cub_c_type
	(const int node_type,         ///< The type of cubature nodes.
	 const struct Simulation* sim ///< \ref Simulation.
	);

/** \brief Constructor for the barycentric coordinates of the p2 reference element of the input type.
 *  \return See brief.
 *
 *  \todo Move this to be with the h-refinement related connectivity information?
 */
static const struct const_Matrix_d* constructor_b_coords
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
static const struct const_Vector_i* constructor_ind_h_b_coords
	(const int e_type,            ///< \ref Element::type.
	 const char ce,               ///< \ref Op_IO::ce.
	 const int ind_h,             ///< The index of the h-refinement.
	 const int ind_ce,            ///< The index of the computational element.
	 const struct Simulation* sim ///< \ref Simulation.
	);

/** \brief Compute scaling factors to be used for the computation of the projected rst coordinates for tensor-product
 *         elements.
 *  \return A statically allocated \ref const_Matrix_T containing the required data. */
static struct const_Matrix_d get_scaling_factors
	(const int e_type, ///< The element type.
	 const char ce,    ///< The computational element type.
	 const int ind_ce  ///< The computational element index.
	);

static int compute_node_type
	(const struct Op_IO* op_io, const struct const_Element* element, const struct Simulation* sim)
{
	switch (op_io->kind) {
	case 's': // fallthrough
	case 'f': // fallthrough
	case 'r': // fallthrough
	case 't': // fallthrough
	case 'g': // fallthrough
	case 'm': // fallthrough
	case 'v':
		return compute_node_type_std(op_io,element,sim);
		break;
	case 'p':
		return NODES_PLOT;
	case 'c':
		if (!sim->collocated)
			return compute_node_type_std(op_io,element,sim);
		else
			return compute_node_type_collocated(op_io,element,sim);
		break;
	default:
		EXIT_ERROR("Unsupported: %c\n",op_io->kind);
		break;
	}
	EXIT_UNSUPPORTED;
	return -1;
}

static int compute_d_nodes (const char ce, const int elem_d)
{
	switch (ce) {
		case 'v': return elem_d;              break;
		case 'f': return elem_d-1;            break;
		case 'e':
			switch (elem_d) {
			case 3:
				return elem_d-2;
				break;
			case 2: // fallthrough
			case 1:
				return compute_d_nodes('f',elem_d);
				break;
			default:
				EXIT_ERROR("Unsupported: %d.\n",elem_d);
				break;
			}
		default:
			EXIT_ERROR("Unsupported: %c\n",ce);
			break;
	}
	EXIT_UNSUPPORTED;
	return -1;
}

static int compute_p_nodes (const struct Op_IO* op_io, const int node_type, const struct Simulation* sim)
{
	const int node_kind = op_io->kind;
	switch (node_kind) {
	case 's': // fallthrough
	case 'f': // fallthrough
	case 'r': // fallthrough
	case 't': // fallthrough
	case 'g': // fallthrough
	case 'm': // fallthrough
	case 'p': // fallthrough
	case 'v':
		return compute_p_basis(op_io,sim);
		break;
	case 'c': {
		// Note that cubature order depends on the node_type (May be p_op or p_c_x*p_op+p_c_p).
		const int cub_c_type = compute_cub_c_type(node_type,sim),
		          p_op       = op_io->p_op;
		if (cub_c_type == CUB_C_COL) {
			return p_op;
		} else {
			const int ind_p_c = ( op_io->sc == 's' ? 0 : 1 );
			const int p_c_x = sim->p_c_x[ind_p_c],
			          p_c_p = sim->p_c_p[ind_p_c];

			const int cub_order = p_c_x*p_op+p_c_p;

			if (cub_c_type == CUB_C_STD)
				return cub_order;
			else if (cub_c_type == CUB_C_DIV2)
				return cub_order/2;
			else
				EXIT_UNSUPPORTED;
		}
		break;
	} default:
		EXIT_ERROR("Unsupported: %c\n",node_kind);
		break;
	}
	EXIT_UNSUPPORTED;
	return -1;
}

static const struct const_Matrix_d* constructor_rst_ve
	(const int s_type, const int d_i, const int d_io, const int ind_h, const int ind_ce, const char ce,
	 const struct Simulation* sim)
{
/// \todo Clean up this function.
UNUSED(d_io);
	const int e_type = compute_elem_from_super_type(s_type,d_i);
	if (e_type == PYR && ce != 'v')
		EXIT_ADD_SUPPORT;

	const struct const_Matrix_d* b_coords       = constructor_b_coords(e_type); // destructed
	const struct const_Vector_i* ind_h_b_coords = constructor_ind_h_b_coords(e_type,ce,ind_h,ind_ce,sim); // destructed

	const struct const_Matrix_d* vv0_vvs_vvs =
		constructor_subset_const_Matrix_d(b_coords,ind_h_b_coords); // destructed

	constructor_Nodes_fptr constructor_Nodes = get_constructor_Nodes_by_super_type(s_type);
// May need to use d_io below for fv, ev
	const struct const_Nodes* nodes_ve = constructor_Nodes(d_i,1,NODES_VERTEX); // destructed
//	const struct const_Nodes* nodes_ve = constructor_Nodes(d_io,1,NODES_VERTEX); // destructed

	const struct const_Matrix_d* rst_ve  = constructor_mm_NN1C_const_Matrix_d(vv0_vvs_vvs,nodes_ve->rst); // returned

	destructor_const_Matrix_d(b_coords);
	destructor_const_Vector_i(ind_h_b_coords);
	destructor_const_Matrix_d(vv0_vvs_vvs);
	destructor_const_Nodes(nodes_ve);

	return rst_ve;
}

static struct Matrix_d* constructor_rst_proj
	(const int e_type, const int d_i, const int d_io, const int ind_h, const int ind_ce, const char ce,
	 const struct const_Matrix_d* b_coords, const struct Simulation* sim)
{
	assert(d_io > d_i);
	assert(ce == 'f' || ce == 'e');

	const ptrdiff_t n_n = b_coords->ext_0;
	struct Matrix_d* rst_proj = constructor_zero_Matrix_d('C',n_n,d_i); // returned

	const struct const_Element*const element = get_element_by_type(sim->elements,e_type);

	const struct const_Vector_i*const b_ve = ( ce == 'f' ? element->f_ve->data[ind_ce] : element->e_ve->data[ind_ce] );

	const int s_type    = compute_super_from_elem_type(e_type),
	          e_type_ce = compute_elem_type_sub_ce(e_type,ce,ind_ce),
	          s_type_ce = compute_super_from_elem_type(e_type_ce);

	assert(ind_h == 0); // Ensure that all is working as expected otherwise.
	if (strcmp(sim->geom_blending[s_type],"gordon_hall") == 0) {
		const struct const_Matrix_d scaling_factors = get_scaling_factors(e_type,ce,ind_ce);
		mm_d('N','N',1.0,0.0,b_coords,&scaling_factors,rst_proj);
	} else if (strcmp(sim->geom_blending[s_type],"szabo_babuska_gen") == 0) {
		const struct const_Matrix_d* rst_ve_d_ce = constructor_rst_ve(s_type_ce,d_i,d_i,ind_h,0,'v',sim); // d.
		for (int n = 0; n < n_n; ++n) {
			const double*const data_b_coords = get_row_const_Matrix_d(n,b_coords);
			for (int d = 0; d < d_i; ++d) {
				const double*const data_rst_ve = get_col_const_Matrix_d(d,rst_ve_d_ce);
				double*const data_rst_proj     = get_col_Matrix_d(d,rst_proj);
				for (int ve = 0; ve < b_ve->ext_0; ++ve)
					data_rst_proj[n] += data_b_coords[b_ve->data[ve]]*data_rst_ve[ve];
			}
		}
		destructor_const_Matrix_d(rst_ve_d_ce);
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

static int compute_node_type_std
	(const struct Op_IO* op_io, const struct const_Element* element, const struct Simulation* sim)
{
	const char node_kind = op_io->kind,
	           node_ce   = op_io->ce;
	const int s_type = compute_super_type_op(node_ce,op_io->h_op,element);
	switch (node_kind) {
	case 's': // fallthrough
	case 'f': // fallthrough
	case 'r': // fallthrough
	case 't':
		if (strcmp(sim->nodes_interp[s_type],"GL") == 0) {
			return NODES_GL;
		} else if (strcmp(sim->nodes_interp[s_type],"GLL") == 0) {
			if (op_io->p_op > 0)
				return NODES_GLL;
			else
				return NODES_GL;
		} else if (strcmp(sim->nodes_interp[s_type],"AO") == 0) {
			return NODES_AO;
		} else if (strcmp(sim->nodes_interp[s_type],"WSH") == 0) {
			return NODES_WSH;
		} else if (strcmp(sim->nodes_interp[s_type],"EQ") == 0) {
			return NODES_EQ;
		} else {
			EXIT_ERROR("Unsupported: %s\n",sim->nodes_interp[s_type]);
		}
	case 'v':
		assert(op_io->p_op >= 0 && op_io->p_op <= 2);
		return NODES_VERTEX;
		break;
	case 'm': // fallthrough
	case 'g':
		switch (s_type) {
			case ST_TP:  return NODES_GLL; break;
			case ST_SI:  return NODES_AO;  break;
			case ST_PYR: return NODES_GLL; break;
			default:     EXIT_ERROR("Unsupported: %d\n",s_type); break;
		}
		break;
	case 'c':
		switch (s_type) {
		case ST_TP:
			return NODES_GL;
			break;
		case ST_SI:
			switch (DIM) {
			case 2:
				switch (node_ce) {
					case 'v': return NODES_WV; break;
					case 'f': return NODES_GL; break;
					case 'e': // fallthrough
					default:  EXIT_ERROR("Unsupported: %c\n",node_ce); break;
				}
				break;
			case 3:
				switch (node_ce) {
					case 'v': return NODES_WV; break;
					case 'f': return NODES_WV; break;
					case 'e': return NODES_GL; break;
					default:  EXIT_ERROR("Unsupported: %c\n",node_ce); break;
				}
				break;
			default:
				EXIT_ERROR("Unsupported: %d\n",DIM); break;
				break;
			}
			break;
		case ST_PYR: {
			const int out_e_type = compute_elem_type_sub_ce(element->type,node_ce,op_io->h_op);

			switch (out_e_type) {
				case TRI:  return NODES_WV;  break;
				case QUAD: return NODES_GL;  break;
				case TET:  return NODES_WV;  break;
				case PYR:  return NODES_GJW; break;
				default:   EXIT_ERROR("Unsupported: %d\n",out_e_type);
			}
			break;
		} default:
			EXIT_ERROR("Unsupported: %d\n",s_type); break;
			break;
		}
		break;
	default:
		EXIT_ERROR("Unsupported: %c\n",node_kind); break;
		break;
	}
	EXIT_UNSUPPORTED;
	return -1;
}

static int compute_node_type_collocated
	(const struct Op_IO* op_io, const struct const_Element* element, const struct Simulation* sim)
{
	const char node_kind = op_io->kind,
	           node_ce   = op_io->ce;
	const int s_type = element->s_type;
	switch (node_kind) {
	case 'c':
		switch (s_type) {
		case ST_TP:
			if (strcmp(sim->nodes_interp[s_type],"GL") == 0)
				return NODES_GL;
			else if (strcmp(sim->nodes_interp[s_type],"GLL") == 0)
				return NODES_GLL;
			else
				EXIT_ERROR("Unsupported: %s\n",sim->nodes_interp[s_type]);
			break;
		case ST_SI:
			if (strcmp(sim->nodes_interp[s_type],"WSH") != 0)
				EXIT_ERROR("Unsupported: %s\n",sim->nodes_interp[s_type]);

			switch (DIM) {
			case 2:
				switch (node_ce) {
					case 'v': return NODES_WSH; break;
					case 'f': return NODES_GL;  break;
					case 'e': // fallthrough
					default:  EXIT_ERROR("Unsupported: %c\n",node_ce); break;
				}
				break;
			case 3:
				// If wedges are present, WSH nodes must be used for the face cubature for consistency with
				// tet faces. If not present, the higher strength WV nodes are used.
				if (wedges_present(sim->elements))
					return NODES_WSH;
				else
					return NODES_WV;
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
		EXIT_ERROR("Unsupported: %c\n",node_kind); break;
		break;
	}
	EXIT_UNSUPPORTED;
	return -1;
}

static int compute_cub_c_type (const int node_type, const struct Simulation* sim)
{
	switch (node_type) {
	case NODES_GL: // fallthrough
	case NODES_GLL:
		if (!sim->collocated)
			return CUB_C_DIV2;
		else
			return CUB_C_COL;
		break;
	case NODES_WSH:
		assert(sim->collocated);
		return CUB_C_COL;
		break;
	case NODES_WV:
		return CUB_C_STD;
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",node_type);
		break;
	}
}

static const struct const_Matrix_d* constructor_b_coords (const int e_type)
{
	/* `ext_0` is not necessarily equal to the number of p2 rst coordinate nodes for the `e_type`. Additional nodes
	 * may be included for certain refinements (such as splitting a TET into 12 TETs). */
	ptrdiff_t ext_0 = -1,
	          ext_1 = -1; // Number of vertices of the `e_type`.
	const double* b_coords = NULL;

	switch (e_type) {
	case POINT:
		ext_0 = 1;
		ext_1 = 1;
		b_coords = (double[])
			{ 1.0, };
		break;
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
	case QUAD:
		ext_0 = 9;
		ext_1 = 4;
		b_coords = (double[])
			{ 1.0 , 0.0 , 0.0 , 0.0 ,
			  0.5 , 0.5 , 0.0 , 0.0 ,
			  0.0 , 1.0 , 0.0 , 0.0 ,
			  0.5 , 0.0 , 0.5 , 0.0 ,
			  0.25, 0.25, 0.25, 0.25,
			  0.0 , 0.5 , 0.0 , 0.5 ,
			  0.0 , 0.0 , 1.0 , 0.0 ,
			  0.0 , 0.0 , 0.5 , 0.5 ,
			  0.0 , 0.0 , 0.0 , 1.0 ,};
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

static const struct const_Vector_i* constructor_ind_h_b_coords
	(const int e_type, const char ce, const int ind_h, const int ind_ce, const struct Simulation* sim)
{
	const struct const_Element* element = get_element_by_type(sim->elements,e_type);

	ptrdiff_t ext_0 = -1;
	const int* ind_h_b_coords = NULL;

	switch (e_type) {
	case POINT:
		switch (ce) {
		case 'v': // fallthrough
		case 'f': // fallthrough
		case 'e': {
			enum { n_ve_ce = 1, n_ref_max = 1, };
			assert(n_ve_ce == element->n_ve);
			assert(n_ref_max == element->n_ref_max_v);
			assert(ind_h < n_ref_max);

			static const int ind_h_b_coords_all[n_ref_max][n_ve_ce] =
				{{ 0, },
				};

			ext_0 = n_ve_ce;
			ind_h_b_coords = ind_h_b_coords_all[ind_h];
			break;
		} default:
			EXIT_ERROR("Unsupported: %c.\n",ce);
			break;
		}
		break;
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
			enum { n_ve_ce = 1, n_ref_max = 1, n_ce_max = 2, };
			assert(n_ref_max == element->n_ref_max_f);
			assert(n_ce_max == element->n_f);
			assert(ind_h < n_ref_max);
			assert(ind_ce < n_ce_max);

			static const int ind_h_b_coords_all[n_ce_max][n_ref_max][n_ve_ce] =
				{{{0,} // f0
				 },
				 {{2,} // f1
				 },
				};

			ext_0 = n_ve_ce;
			ind_h_b_coords = ind_h_b_coords_all[ind_ce][ind_h];
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
			enum { n_ve_ce = 2, n_ref_max = 3, n_ce_max = 3, };
			assert(n_ref_max == element->n_ref_max_f);
			assert(n_ce_max == element->n_f);
			assert(ind_h < n_ref_max);
			assert(ind_ce < n_ce_max);

			static const int ind_h_b_coords_all[n_ce_max][n_ref_max][n_ve_ce] =
				{{{2,5,}, // f0
				  {2,4,}, {4,5,},
				 },
				 {{0,5,}, // f1
				  {0,3,}, {3,5,},
				 },
				 {{0,2,}, // f2
				  {0,1,}, {1,2,},
				 },
				};

			ext_0 = n_ve_ce;
			ind_h_b_coords = ind_h_b_coords_all[ind_ce][ind_h];
			break;
		} default:
			EXIT_ERROR("Unsupported: %c.\n",ce);
			break;
		}
		break;
	case QUAD:
		switch (ce) {
		case 'v': {
			// Should only be used for operators from lower to higher dimensional computational elements and thus
			// not be required for h-refinement.
			enum { n_ve_ce = 4, n_ref_max = 1, };
			assert(n_ve_ce == element->n_ve);
			assert(ind_h == 0);

			static const int ind_h_b_coords_all[n_ref_max][n_ve_ce] =
				{{ 0, 2, 6, 8,},
				};

			ext_0 = n_ve_ce;
			ind_h_b_coords = ind_h_b_coords_all[ind_h];
			break;
		} case 'f': { // fallthrough
		} case 'e': {
			enum { n_ve_ce = 2, n_ref_max = 3, n_ce_max = 4, };
			assert(n_ref_max == element->n_ref_max_f);
			assert(n_ce_max == element->n_f);
			assert(ind_h < n_ref_max);
			assert(ind_ce < n_ce_max);

			static const int ind_h_b_coords_all[n_ce_max][n_ref_max][n_ve_ce] =
				{{{0,6,}, // f0
				  {0,3,}, {3,6,},
				 },
				 {{2,8,}, // f1
				  {2,5,}, {5,8,},
				 },
				 {{0,2,}, // f2
				  {0,1,}, {1,2,},
				 },
				 {{6,8,}, // f3
				  {6,7,}, {7,8,},
				 },
				};

			ext_0 = n_ve_ce;
			ind_h_b_coords = ind_h_b_coords_all[ind_ce][ind_h];
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
			enum { n_ve_ce = 3, n_ref_max = 5, n_ce_max = 4, };
			assert(n_ref_max == element->n_ref_max_f);
			assert(n_ce_max == element->n_f);
			assert(ind_h < n_ref_max);
			assert(ind_ce < n_ce_max);

			static const int ind_h_b_coords_all[n_ce_max][n_ref_max][n_ve_ce] =
				{{{ 2, 5, 9,}, // f0
				  { 2, 4, 7,}, { 4, 5, 8,}, { 7, 8, 9,}, { 8, 7, 4,},
				 },
				 {{ 0, 5, 9,}, // f1
				  { 0, 3, 6,}, { 3, 5, 8,}, { 6, 8, 9,}, { 8, 6, 3,},
				 },
				 {{ 0, 2, 9,}, // f2
				  { 0, 1, 6,}, { 1, 2, 7,}, { 6, 7, 9,}, { 7, 6, 1,},
				 },
				 {{ 0, 2, 5,}, // f3
				  { 0, 1, 3,}, { 1, 2, 4,}, { 3, 4, 5,}, { 4, 3, 1,},
				 },
				};

			ext_0 = n_ve_ce;
			ind_h_b_coords = ind_h_b_coords_all[ind_ce][ind_h];
		} case 'e': {
			enum { n_ve_ce = 2, n_ref_max = 3, n_ce_max = 6, };
			assert(n_ref_max == element->n_ref_max_f);
			assert(n_ce_max == element->n_f);
			assert(ind_h < n_ref_max);
			assert(ind_ce < n_ce_max);

			static const int ind_h_b_coords_all[n_ce_max][n_ref_max][n_ve_ce] =
				{{{ 2, 5,}, // e0
				  { 2, 4,}, { 4, 5,},
				 },
				 {{ 0, 5,}, // e1
				  { 0, 3,}, { 3, 5,},
				 },
				 {{ 0, 2,}, // e2
				  { 0, 1,}, { 1, 2,},
				 },
				 {{ 0, 9,}, // e3
				  { 0, 6,}, { 6, 9,},
				 },
				 {{ 2, 9,}, // e4
				  { 2, 7,}, { 7, 9,},
				 },
				 {{ 5, 9,}, // e5
				  { 5, 8,}, { 8, 9,},
				 },
				};

			ext_0 = n_ve_ce;
			ind_h_b_coords = ind_h_b_coords_all[ind_ce][ind_h];
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
			enum { n_ve_ce = 4, n_ref_max = 5, n_ce_max = 5, };
			assert(n_ref_max == element->n_ref_max_f);
			assert(n_ce_max == element->n_f);
			assert(ind_h < n_ref_max);
			assert(ind_ce < n_ce_max);

			static const int ind_h_b_coords_all[n_ce_max][n_ref_max][n_ve_ce] =
				{{{ 0, 6,13,-1,}, // f0
				  { 0, 3, 9,-1,}, { 3, 6,11,-1,}, { 9,11,13,-1,}, {11, 9, 3,-1,},
				 },
				 {{ 2, 8,13,-1,}, // f1
				  { 2, 5,10,-1,}, { 5, 8,12,-1,}, {10,12,13,-1,}, {12,10, 5,-1,},
				 },
				 {{ 0, 2,13,-1,}, // f2
				  { 0, 1, 9,-1,}, { 1, 2,10,-1,}, { 9,10,13,-1,}, {10, 9, 1,-1,},
				 },
				 {{ 6, 8,13,-1,}, // f3
				  { 6, 7,11,-1,}, { 7, 8,12,-1,}, {11,12,13,-1,}, {12,11, 7,-1,},
				 },
				 {{ 0, 2, 6, 8,}, // f4
				  { 0, 1, 3, 4,}, { 1, 2, 4, 5,}, { 3, 4, 6, 7,}, { 4, 5, 7, 8,},
				 },
				};
			static const int n_ve_ce_all[n_ce_max] = {3,3,3,3,4,};

			ext_0 = n_ve_ce_all[ind_ce];
			ind_h_b_coords = ind_h_b_coords_all[ind_ce][ind_h];
		} case 'e': {
			enum { n_ve_ce = 2, n_ref_max = 3, n_ce_max = 8, };
			assert(n_ref_max == element->n_ref_max_f);
			assert(n_ce_max == element->n_f);
			assert(ind_h < n_ref_max);
			assert(ind_ce < n_ce_max);

			static const int ind_h_b_coords_all[n_ce_max][n_ref_max][n_ve_ce] =
				{{{ 0, 6,}, // e0
				  { 0, 3,}, { 3, 6,},
				 },
				 {{ 2, 8,}, // e1
				  { 2, 5,}, { 5, 8,},
				 },
				 {{ 0, 2,}, // e2
				  { 0, 1,}, { 1, 2,},
				 },
				 {{ 6, 8,}, // e3
				  { 6, 7,}, { 7, 8,},
				 },
				 {{ 0,13,}, // e4
				  { 0, 9,}, { 9,13,},
				 },
				 {{ 2,13,}, // e5
				  { 0,10,}, {10,13,},
				 },
				 {{ 6,13,}, // e6
				  { 0,11,}, {11,13,},
				 },
				 {{ 8,13,}, // e7
				  { 0,12,}, {12,13,},
				 },
				};

			ext_0 = n_ve_ce;
			ind_h_b_coords = ind_h_b_coords_all[ind_ce][ind_h];
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

static struct const_Matrix_d get_scaling_factors (const int e_type, const char ce, const int ind_ce)
{
	assert(ce == 'f' || ce == 'e');

	static struct Matrix_d factors = { .layout = 'C', .ext_0 = 0, .ext_1 = 0, .owns_data = false, .data = NULL, };
	switch (e_type) {
	case QUAD:
		factors.ext_0 = 4;
		factors.ext_1 = 1;
		switch (ind_ce) {
		case 0: // fallthrough
		case 1: {
			static double data[] = (double[]) { -1.0, -1.0,  1.0,  1.0, };
			factors.data = data;
			break;
		} case 2: // fallthrough
		  case 3: {
			static double data[] = (double[]) { -1.0,  1.0, -1.0,  1.0, };
			factors.data = data;
			break;
		} default:
			EXIT_ERROR("Unsupported: %d\n",ind_ce);
			break;
		}
		break;
	case HEX:
		EXIT_ADD_SUPPORT;
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",e_type);
		break;
	}
	return *(struct const_Matrix_d*)&factors;
}
