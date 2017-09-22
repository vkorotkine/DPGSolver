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
 *  The definitions provided here are related to parameters needed to specify the element operators:
 *  - type  (OP_T_*);
 *  - range (OP_R_*).
 *
 *  \todo Remove unused definitions.
 */

///\{ \name Operator related parameters indices.
#define OP_ORDER_MAX 7 ///< Maximum operator order (d, ce_o, ce_i, h_o, h_i, p_o, p_i).

#define OP_IND_I     1 ///< Index of the input  parameter relative to OP_IND_*.
#define OP_IND_O     0 ///< Index of the output parameter relative to OP_IND_*.

#define OP_IND_CE    1 ///< Index of the first computational element parameter.
#define OP_IND_H     3 ///< Index of the first h-refinement parameter.
#define OP_IND_P     5 ///< Index of the first p-refinement parameter.
///\}


// Operator Type Options ******************************************************************************************** //

///\{ \name Operator types.
#define OP_T_CV 100 ///< Coefficients to values.
#define OP_T_CC 101 ///< Coefficients to coefficients.
#define OP_T_VV 102 ///< Values       to values.
#define OP_T_VC 103 ///< Values       to coefficients.

//#define OP_T_DG_WEAK_VV 110 ///< (D)iscontinuous-(G)alerkin (Weak) (V)olume to (V)olume.
///\}


// Operator Range Options ******************************************************************************************* //

///\{ \name Dimension range options.
#define OP_R_D_0   0  ///< No dimensions (NULL range)
#define OP_R_D_ALL 10 ///< All dimensions.
///\}

///\{ \name Computational element range options.
#define OP_R_CE_VV 10 ///< Volume to volume.
#define OP_R_CE_VF 11 ///< Volume to face.
#define OP_R_CE_VE 12 ///< Volume to edge.
#define OP_R_CE_FV 20 ///< Face to volume.
#define OP_R_CE_FF 21 ///< Face to face.
#define OP_R_CE_FE 22 ///< Face to edge.
#define OP_R_CE_EV 30 ///< Edge to volume.
#define OP_R_CE_EF 31 ///< Edge to face.
#define OP_R_CE_EE 32 ///< Edge to edge.
///\}

///\{ \name Refinement (h) range options.
#define OP_R_H_1  10 ///< Standard operator (no h-refinement) only.
#define OP_R_H_CF 11 ///< Coarse to fine + Standard.
#define OP_R_H_FC 12 ///< Fine to coarse.
///\}

///\{ \name Order (p) range options.
#define OP_R_P_1   10 ///< Order = 1 only.
#define OP_R_P_PM0 11 ///< Order = p_reference +/- 0
#define OP_R_P_PM1 12 ///< Order = p_reference +/- 1
#define OP_R_P_ALL 13 ///< Order = 0:p_max
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
