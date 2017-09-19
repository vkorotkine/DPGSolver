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
 *  It is necessary that the values of *_0 parameters are 0.
 */

///\{ \name Dimension range options.
#define RANGE_D_0   0  ///< No dimensions (NULL range)
#define RANGE_D_ALL 10 ///< All dimensions.
///\}

///\{ \name Face range options.
#define RANGE_F_0   0  ///< No faces (NULL range)
#define RANGE_F_ALL 10 ///< All faces.
///\}

///\{ \name Order (p) range options.
#define RANGE_P_1   10 ///< Order = 1 only.
#define RANGE_P_PM0 11 ///< Order = p_reference +/- 0
#define RANGE_P_PM1 12 ///< Order = p_reference +/- 1
#define RANGE_P_ALL 13 ///< Order = 0:p_max
///\}

///\{ \name Refinement (h) range options.
#define RANGE_H_1   10 ///< Standard operator (no h-refinement) only.
#define RANGE_H_SF  11 ///< Sum-factorized refinement.
#define RANGE_H_ALL 12 ///< All refinements.
///\}


///\{ \name Operator types.
#define OP_T_CC 100 ///< Coefficients to coefficients.
#define OP_T_CV 101 ///< Coefficients to values.
#define OP_T_VC 102 ///< Values       to coefficients.
#define OP_T_VV 103 ///< Values       to values.

#define OP_T_DG_WEAK_VV 110 ///< (D)iscontinuous-(G)alerkin (Weak) (V)olume to (V)olume.
///\}

#endif // DPG__definitions_element_operators_h__INCLUDED
