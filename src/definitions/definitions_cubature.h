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

#ifndef DPG__definitions_cubature_h__INCLUDED
#define DPG__definitions_cubature_h__INCLUDED
/** \file
 *  \brief Provides the definitions related to the cubature.
 */

///\{ \name Nodes available for several element super types.
#define CUB_VERTEX 101 ///< Vertex coordinates.

#define CUB_EQ 111 ///< Equally spaced.
#define CUB_WV 112 ///< Witherden-Vincent.
///\}

///\{ \name The available tensor-product nodes.
#define CUB_GL  201 ///< Gauss-Legendre.
#define CUB_GLL 202 ///< Gauss-Lobatto-Legendre.
///\}

///\{ \name The available simplex nodes.
#define CUB_AO  301 ///< Alpha-optimized.
#define CUB_WSH 302 ///< Williams-Shunn-Ham.
///\}

///\{ \name The available pyramid nodes.
#define CUB_GLW    401 ///< Gauss-Legendre with Weights.
#define CUB_GLLW   402 ///< Gauss-Lobatto-Legendre with Weights.
#define CUB_GJW    403 ///< Gauss-Jacobi with Weights.
#define CUB_WVHToP 404 ///< Witherden-Vincent, HEX to PYR.
///\}

#endif // DPG__definitions_cubature_h__INCLUDED
