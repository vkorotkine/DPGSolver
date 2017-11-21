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

#ifndef DPG__definitions_nodes_h__INCLUDED
#define DPG__definitions_nodes_h__INCLUDED
/** \file
 *  \brief Provides the definitions related to the reference element nodes.
 *
 *  Definitions are provided for:
 *  - the node types;
 *  - the cubature order from reference order indicators.
 */

///\{ \name Nodes available for several element super types.
#define NODES_VERTEX 101 ///< Vertex coordinates.
#define NODES_PLOT   102 ///< Plotting nodes.

#define NODES_EQ 111 ///< Equally spaced.
#define NODES_WV 112 ///< Witherden-Vincent.
///\}

///\{ \name The available tensor-product nodes.
#define NODES_GL  201 ///< Gauss-Legendre.
#define NODES_GLL 202 ///< Gauss-Lobatto-Legendre.
///\}

///\{ \name The available simplex nodes.
#define NODES_AO  301 ///< Alpha-optimized.
#define NODES_WSH 302 ///< Williams-Shunn-Ham.
///\}

///\{ \name The available pyramid nodes.
#define NODES_GLW    401 ///< Gauss-Legendre with Weights.
#define NODES_GLLW   402 ///< Gauss-Lobatto-Legendre with Weights.
#define NODES_GJW    403 ///< Gauss-Jacobi with Weights.
#define NODES_WVHToP 404 ///< Witherden-Vincent, HEX to PYR.
///\}


///\{ \name The available options for the required processing to get from p_ref to the cubature order.
#define CUB_C_STD  1001 ///< Standard:     p_cub_std  = p_c_x*p_ref+p_c_p.
#define CUB_C_DIV2 1002 ///< Divided by 2: p_cub_div2 = p_cub_std/2.
#define CUB_C_COL  1003 ///< Collocated:   p_cub_col  = p_ref.
///\}

#endif // DPG__definitions_nodes_h__INCLUDED
