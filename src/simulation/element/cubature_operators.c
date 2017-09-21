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

#include <assert.h>
#include <string.h>
#include "gsl/gsl_math.h"

#include "macros.h"

#include "element_operators.h"

// Static function declarations ************************************************************************************* //

/** \brief Compute the node type of the cubature based on the kind of operator.
 *  \return See brief. */
static int compute_node_type_cub
	(const char cub_kind,         ///< \ref Op_IO::kind.
	 const int s_type,            ///< \ref Element::s_type.
	 const struct Simulation *sim ///< \ref Simulation.
	);

/** \brief Compute the dimension of the cubature based off of the computational element and the simulation dimension.
 *  \return See brief. */
static int compute_d_cub
	(const char ce,  ///< \ref Op_IO::ce.
	 const int sim_d ///< \ref Simulation::d.
	);

/** \brief Compute the order of the cubature nodes to be computed based on the reference order and the kind of operator.
 *  \return See brief. */
static int compute_p_cub
	(const int p_ref,             ///< The reference polynomial order.
	 const char cub_kind,         ///< \ref Op_IO::kind.
	 const char cub_sc,           ///< \ref Op_IO:sc.
	 const int s_type,            ///< \ref Element::s_type.
	 const int node_type,         ///< \ref Cubature::node_type.
	 const struct Simulation *sim ///< \ref Simulation.
	);

// Constructor functions ******************************************************************************************** //

const struct const_Cubature* constructor_const_Cubature_h
	()
{
	const int node_type = compute_node_type_cub(op_io[ind_io].kind,s_type,sim);
	const int d         = compute_d_cub(op_io[ind_io].ce,sim->d);
	const int p         = compute_p_cub(p_op,op_io[ind_io].kind,op_io[ind_io].sc,s_type,node_type,sim);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Compute the node type of the cubature based on the kind of operator for standard (non-collocated) nodes.
 *  \return See brief. */
static int compute_node_type_cub_std
	(const char cub_kind,         ///< \ref Op_IO::kind.
	 const int s_type,            ///< \ref Element::s_type.
	 const struct Simulation *sim ///< \ref Simulation.
	);

/** \brief Compute the node type of the cubature based on the kind of operator for collocated nodes.
 *  \return See brief. */
static int compute_node_type_cub_collocated
	(const char cub_kind,         ///< \ref Op_IO::kind.
	 const int s_type,            ///< \ref Element::s_type.
	 const struct Simulation *sim ///< \ref Simulation.
	);

static int compute_node_type_cub (const char cub_kind, const int s_type, const struct Simulation *sim)
{
	switch (cub_kind) {
	case 'g':
		return compute_node_type_cub_std(cub_kind,s_type,sim);
		break;
	case 'c':
		if (!sim->collocated)
			return compute_node_type_cub_std(cub_kind,s_type,sim);
		else
			return compute_node_type_cub_collocated(cub_kind,s_type,sim);
		break;
	default:
		EXIT_ERROR("Unsupported: %c\n",cub_kind);
		break;
	}
	EXIT_UNSUPPORTED;
	return -1;
}

static int compute_d_cub (const char ce, const int sim_d)
{
	assert(sim_d > 0);

	switch (ce) {
		case 'v': return sim_d;              break;
		case 'f': return sim_d-1;            break;
		case 'e': return GSL_MAX(sim_d-2,0); break;
		default:
			EXIT_ERROR("Unsupported: %c\n",ce);
			break;
	}
	EXIT_UNSUPPORTED;
	return -1;
}

static int compute_p_cub
	(const int p_ref, const char cub_kind, const char cub_sc, const int s_type, const int node_type,
	 const struct Simulation *sim)
{
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

// Level 1 ********************************************************************************************************** //

static int compute_node_type_cub_std (const char cub_kind, const int s_type, const struct Simulation *sim)
{
	switch (cub_kind) {
	case 's':
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
		switch (s_type) {
			case ST_TP:  return CUB_GLL;   break;
			case ST_SI:  return CUB_AO;    break;
			case ST_PYR: return CUB_GLL;   break;
			default:     EXIT_ERROR("Unsupported: %d\n",s_type); break;
		}
		break;
	case 'm':
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
					case 'e': /* fallthrough */
/// \todo Replace all possible EXIT_UNSUPPORTED with EXIT_ERROR with name.
					default:  EXIT_UNSUPPORTED; break;
				}
				break;
			case 3:
				switch (cub_ce) {
					case 'v': return CUB_WV; break;
					case 'f': return CUB_WV; break;
					case 'e': return CUB_GL; break;
					default:  EXIT_UNSUPPORTED; break;
				}
				break;
			default:
				EXIT_UNSUPPORTED;
			}
			break;
		case ST_PYR:
			return CUB_GJW;
			break;
		default:
			EXIT_UNSUPPORTED;
			break;
		}
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}
	EXIT_UNSUPPORTED;
	return -1;
}

static int compute_node_type_cub_collocated (const char cub_kind, const int s_type, const struct Simulation *sim)
{
	switch (cub_kind) {
	case 'c':
		switch (s_type) {
		case ST_TP:
			if (strcmp(sim->nodes_interp[s_type],"GL") == 0)
				return CUB_GL;
			else if (strcmp(sim->nodes_interp[s_type],"GLL") == 0)
				return CUB_GLL;
			else
				EXIT_UNSUPPORTED;
			break;
		case ST_SI:
			if (strcmp(sim->nodes_interp[s_type],"WSH") != 0)
				EXIT_UNSUPPORTED;

			switch (sim->d) {
			case 2:
				switch (cub_ce) {
					case 'v': return CUB_WSH; break;
					case 'f': return CUB_GL;  break;
					case 'e': /* fallthrough */
					default:  EXIT_UNSUPPORTED; break;
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
		case ST_PYR: /* fallthrough */
		default:
			EXIT_UNSUPPORTED;
			break;
		}
		break;
	case 's': /* fallthrough */
	case 'g': /* fallthrough */
	case 'm': /* fallthrough */
	default:
		EXIT_UNSUPPORTED;
		break;
	}
	EXIT_UNSUPPORTED;
	return -1;
}
