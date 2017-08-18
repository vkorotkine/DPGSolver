// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__constants_elements_h__INCLUDED
#define DPG__constants_elements_h__INCLUDED
/**	\file
 *	Define constants associated with elements.
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
#define NFMAX           6  ///< (MAX)imum (N)umber of (F)aces.           "HEX"
#define NFVEMAX         4  ///< (MAX)imum (N)unber of (F)ace (VE)rtices. "QUAD"
///\}

#endif // DPG__constants_elements_h__INCLUDED
