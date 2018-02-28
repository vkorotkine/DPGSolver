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

#ifndef DPG__definitions_adaptation_h__INCLUDED
#define DPG__definitions_adaptation_h__INCLUDED
/** \file
 *  \brief Provides the definitions relating to the adaptation functionality.
 */

///\{ \name Flag for computational elements which should not be adapted.
#define ADAPT_NONE 0
///\}

///\{ \name Supported adaptation options.
#define ADAPT_0  0
#define ADAPT_P  1
#define ADAPT_H  2
#define ADAPT_HP 3
///\}

///\{ \name Supported adaptation types for computational elements.
#define ADAPT_P_REFINE 101
#define ADAPT_P_COARSE 102
#define ADAPT_H_REFINE 103
#define ADAPT_H_COARSE 104

#define ADAPT_HP_CC 105 ///< Coarsen in h, coarsen in p.
#define ADAPT_HP_CR 106 ///< Coarsen in h, refine  in p.
#define ADAPT_HP_RC 107 ///< Refine  in h, coarsen in p.
#define ADAPT_HP_RR 108 ///< Refine  in h, refine  in p.

#define ADAPT_H_CREATE  111
#define ADAPT_H_DESTROY 112
///\}

///\{ \name Supported adaptation strategies for the entire domain.
#define ADAPT_S_P_REFINE 1011 ///< Uniform p refinement for the entire domain.
#define ADAPT_S_P_COARSE 1012 ///< Uniform p coarsening for the entire domain.
#define ADAPT_S_H_REFINE 1013 ///< Uniform h refinement for the entire domain.
#define ADAPT_S_H_COARSE 1014 ///< Uniform h coarsening for the entire domain.

#define ADAPT_S_XYZ_VE 1021 ///< Adapt around the specified vertices until the required mesh level is reached.
///\}

///\{ \name Support h-refinement types for tetrahedral elements.
#define TET8  8  ///< Middle octohedron split into 4 TETs.
#define TET12 12 ///< Middle octohedron split into 8 TETs.
#define TET6  6  ///< Middle octohedron split into 2 PYRs.
///\}

#endif // DPG__definitions_adaptation_h__INCLUDED
