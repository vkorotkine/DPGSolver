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

#ifndef DPG__definitions_element_operators_h__INCLUDED
#define DPG__definitions_element_operators_h__INCLUDED
/** \file
 *  \brief Provides the definitions relating to the element operators.
 *
 *  \todo Move these comments.
 *  The definitions provided here are related to parameters needed to specify the element operators:
 *  - type     (OP_T_*);
 *  - basis    (OP_B_*);
 *  - cubature (OP_C_*);
 *  - range    (OP_R_*).
 *
 *  Operator names take the general form: [type][0]_[1][2](3)_[4][5](6) where entries in [square brackets] are required
 *  and those in (round brackets) are optional:
 *  - type:
 *  	- cv: coefficients to values
 *  	- cc: coefficients to coefficients
 *  	- vv: values       to values
 *  	- vc: values       to coefficients
 *  - [0]:   value for the order of differentiation (0 for no differentiation).
 *  - [1/4]: character denoting the type of computational element:
 *  	- v: volume
 *  	- f: face
 *  	- e: edge
 *  - [2/5]: character denoting the type of basis/cubature to be used:
 *  	- g: geometry
 *  	- m: metric
 *  	- s: solution
 *  	- p: plotting (Not available as a basis [2])
 *  	- c: cubature (Not available as a basis [2])
 *  - (3/6): character denoting whether the basis/cubature is meant to be used straight or curved elements:
 *  	- s: straight
 *  	- c: curved
 *
 *  The optional straight/curved parameters should be replaced with 'A'll if not present when passed to the
 *  constructor function.
 *
 *  Each operator also has an associated range with a maximum order of \ref OP_ORDER_MAX with the following parameters
 *  (d)(f)[h_o][h_i][p_o][p_i], where entries in square and round brackets are once again required and optional,
 *  respectively.
 */

///\{ \name Operator related parameters
#define OP_ORDER_MAX 6 // d, f, h_o, h_i, p_o, p_i
#define OP_IND_H     2
#define OP_IND_P     4
///\}

// Operator types *************************************************************************************************** //

///\{ \name Operator types.
#define OP_T_CV 100 ///< Coefficients to values.
#define OP_T_CC 101 ///< Coefficients to coefficients.
#define OP_T_VV 102 ///< Values       to values.
#define OP_T_VC 103 ///< Values       to coefficients.

#define OP_T_DG_WEAK_VV 110 ///< (D)iscontinuous-(G)alerkin (Weak) (V)olume to (V)olume.
///\}


///\{ \name Dimension range options.
#define RANGE_D_0   0  ///< No dimensions (NULL range)
#define RANGE_D_ALL 10 ///< All dimensions.
///\}

///\{ \name Face range options.
#define RANGE_F_0   0  ///< No faces (NULL range)
#define RANGE_F_ALL 10 ///< All faces.
///\}

///\{ \name Refinement (h) range options.
#define RANGE_H_1  10 ///< Standard operator (no h-refinement) only.
#define RANGE_H_CF 11 ///< Coarse to fine + Standard.
#define RANGE_H_FC 12 ///< Fine to coarse.
///\}

///\{ \name Order (p) range options.
#define RANGE_P_1   10 ///< Order = 1 only.
#define RANGE_P_PM0 11 ///< Order = p_reference +/- 0
#define RANGE_P_PM1 12 ///< Order = p_reference +/- 1
#define RANGE_P_ALL 13 ///< Order = 0:p_max
///\}

///\{ \name Supported compound range options.
#define RNG_CMP_ALL_0_1_1 ///< D_ALL, F_0, H_1, P_1.
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





#endif // DPG__definitions_element_operators_h__INCLUDED
