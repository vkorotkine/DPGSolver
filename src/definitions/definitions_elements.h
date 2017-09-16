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

///\{ \name Element super types.
#define ST_TP    1
#define ST_SI    2
#define ST_PYR   3
#define ST_WEDGE 4
///\}

///\{ \name Element related values

// Vertex
#define NVEMAX          8  ///< (MAX)imum (N)umber of (VE)rtices for an element. "HEX"

// Face
#define NFMAX           6  ///< (MAX)imum (N)umber of (F)aces.                                 "HEX"
#define NFVEMAX         4  ///< (MAX)imum (N)unber of (F)ace (VE)rtices.                       "QUAD"
#define NFREFMAX        9  ///< (MAX)imum (N)umber of (F)ACE (REF)inements.                    "QUAD"
#define NSUBFMAX        4  ///< (MAX)imum (N)umber of h-adaptive (SUB)-(F)aces (on each FACE). "QUAD/TRI (Isotropic)"
///\}

#endif // DPG__definitions_elements_h__INCLUDED
