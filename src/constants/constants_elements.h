// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__constants_elements_h__INCLUDED
#define DPG__constants_elements_h__INCLUDED
/**	\file
 *	\brief Provides the definition of constants associated with \ref Element containers.
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

///\{ \name Element related values

// Vertex
#define NVEMAX          8  ///< (MAX)imum (N)umber of (VE)rtices for an element. "HEX"

// Face
#define NFMAX           6  ///< (MAX)imum (N)umber of (F)aces.                                 "HEX"
#define NFVEMAX         4  ///< (MAX)imum (N)unber of (F)ace (VE)rtices.                       "QUAD"
#define NFREFMAX        9  ///< (MAX)imum (N)umber of (F)ACE (REF)inements.                    "QUAD"
#define NSUBFMAX        4  ///< (MAX)imum (N)umber of h-adaptive (SUB)-(F)aces (on each FACE). "QUAD/TRI (Isotropic)"
///\}

#endif // DPG__constants_elements_h__INCLUDED
