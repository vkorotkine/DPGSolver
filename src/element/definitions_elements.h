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

#ifndef DPG__definitions_elements_h__INCLUDED
#define DPG__definitions_elements_h__INCLUDED
/**	\file
 *	\brief Provides the definitions relating to the \ref Element container.
 *
 *	The more extensive list of the element types supported by Gmsh can be found under the 'elm-type' header of the
 *	[File formats][gmsh_ff] section of the Gmsh manual. The curved element types are not be supported here however as
 *	the curving treatment is performed within the code.
 *
 *	<!-- References: -->
 *	[gmsh_ff]: http://gmsh.info/doc/texinfo/gmsh.html#File-formats
 */

#include "definitions_core.h"

///\{ \name Gmsh element types
#define POINT 15
#define LINE  1
#define TRI   2
#define QUAD  3
#define TET   4
#define HEX   5
#define WEDGE 6
#define PYR   7
///\}

/**\{ \name The number of element super types which have at least one member which is not formed from a tensor-product
 *          of other super types. */
#define N_ST_STD 3
///\}

///\{ \name Element super types.
#define ST_TP    0 // The values for ST_TP, ST_SI, and ST_PYR are used as indices and should not be changed.
#define ST_SI    1
#define ST_PYR   2
#define ST_WEDGE 10
///\}

///\{ \name Element related values

#if DIM == 1
	#define NVEMAX          2  ///< (MAX)imum (N)umber of (VE)rtices for an element. "LINE"

	#define NEMAX           2  ///< (MAX)imum (N)umber of (E)dges.                                 "LINE"
	#define NFMAX           2  ///< (MAX)imum (N)umber of (F)aces.                                 "LINE"
	#define NFVEMAX         1  ///< (MAX)imum (N)unber of (F)ace (VE)rtices.                       "POINT"
	#define NFREFMAX        1  ///< (MAX)imum (N)umber of (F)ACE (REF)inements.                    "POINT"
	#define NSUBFMAX        1  ///< (MAX)imum (N)umber of h-adaptive (SUB)-(F)aces (on each FACE). "POINT"
#elif DIM == 2
	#define NVEMAX          4  ///< "QUAD"

	#define NEMAX           4  ///< "QUAD"
	#define NFMAX           4  ///< "QUAD"
	#define NFVEMAX         2  ///< "LINE"
	#define NFREFMAX        3  ///< "LINE"
	#define NSUBFMAX        2  ///< "LINE (Isotropic)"
#elif DIM == 3
	#define NVEMAX          8  ///< "HEX"

	#define NEMAX           12 ///< "HEX"
	#define NFMAX           6  ///< "HEX"
	#define NFVEMAX         4  ///< "QUAD"
	#define NFREFMAX        5  ///< "QUAD"
	#define NSUBFMAX        4  ///< "QUAD/TRI (Isotropic)"
#endif
///\}

///\{ \name h-refinement related values
#define NREFMAXPOINT 1
#define NREFMAXLINE  3
#define NREFMAXTRI   5
#define NREFMAXQUAD  5
#define NREFMAXTET   13
#define NREFMAXHEX   9
#define NREFMAXWEDGE 9
#define NREFMAXPYR   11
///\}

///\{ \name Supported operator types. /// \todo Delete if compiling when commented.
//#define OP_V_D0 101
//#define OP_V_D1 102
//#define OP_F_D0 103
//#define OP_F_D1 104
///\}

#endif // DPG__definitions_elements_h__INCLUDED
