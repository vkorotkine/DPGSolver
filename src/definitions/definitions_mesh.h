// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__definitions_mesh_h__INCLUDED
#define DPG__definitions_mesh_h__INCLUDED
/**	\file
 *	\brief Provides the definitions relating to the \ref Mesh container.
 */

///\{ \name The supported domain types.
#define DOM_STRAIGHT   1
#define DOM_CURVED     2
#define DOM_PARAMETRIC 3
///\}

///\{ \name The node coordinate tolerance of the mesh vertices
#define NODETOL_MESH 1.0e-5
///\}

#endif // DPG__definitions_mesh_h__INCLUDED
