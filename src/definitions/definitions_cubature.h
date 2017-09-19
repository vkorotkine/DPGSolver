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
#define CUB_EQ 101 ///< Equally spaced.
#define CUB_WV 102 ///< Witherden-Vincent.
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



///\{ \name The base/multiplier constants for the various indices associated with the cubature node classification.
#define CUB_ENT_BASE 100   ///< Entities.
#define CUB_CE_MULT  1000  ///< Computational elements.
#define CUB_SC_MULT  10000 ///< Straight/Curved.
///\}

///\{ \name Entities for which nodes are used.
#define CUB_ENT_S CUB_ENT_BASE+1 ///< Solution.
#define CUB_ENT_C CUB_ENT_BASE+2 ///< Cubature.
#define CUB_ENT_G CUB_ENT_BASE+3 ///< Geometry.
#define CUB_ENT_M CUB_ENT_BASE+4 ///< Metric.
#define CUB_ENT_P CUB_ENT_BASE+5 ///< Plotting.
///\}

///\{ \name The computational element for which the nodes may be set.
#define CUB_CE_V 1*CUB_CE_MULT ///< Volume.
#define CUB_CE_F 2*CUB_CE_MULT ///< Face.
///\}

///\{ \name The straight/curved indicator.
#define CUB_SC_S 1*CUB_SC_MULT ///< Straight.
#define CUB_SC_C 2*CUB_SC_MULT ///< Curved.
///\}


///\{ \name Compound entities for which nodes are used.
#define CUB_CMP_VGS CUB_CE_V+CUB_ENT_G+CUB_SC_S ///< Volume Geometry which is Straight.
///\}

#endif // DPG__definitions_cubature_h__INCLUDED
