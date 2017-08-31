// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__constants_mesh_h__INCLUDED
#define DPG__constants_mesh_h__INCLUDED
/**	\file
 *	Define constants associated with the mesh.
 */

///\{ \name The supported domain types.
#define DOM_STRAIGHT   1
#define DOM_CURVED     2
#define DOM_PARAMETRIC 3
///\}

///\{ \name The node coordinate tolerance of the mesh vertices
#define NODETOL_MESH 1.0e-5
///\}

#endif // DPG__constants_mesh_h__INCLUDED
