// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

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

#endif // DPG__definitions_cubature_h__INCLUDED
