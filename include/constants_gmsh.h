// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__constants_gmsh_h__INCLUDED
#define DPG__constants_gmsh_h__INCLUDED
/**	\file
 *	Define constants associated with gmsh.
 *
 *	The more extensive list of the element types supported by Gmsh can be found under the 'elm-type' header of the
 *	[File formats][gmsh_ff] section of the Gmsh manual.
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

#endif // DPG__constants_gmsh_h__INCLUDED
