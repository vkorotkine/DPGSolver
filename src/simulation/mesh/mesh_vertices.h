// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__mesh_vertices_h__INCLUDED
#define DPG__mesh_vertices_h__INCLUDED
/**	\file
 *	\brief Provides the interface to mesh vertex containers and functions.
 *
 *	\section s1_mesh_vert Discussion Concerning the use of Exact Geometry
 *
 *	For curved domains which are not mapped, the mesh vertices may optionally be corrected such that they are located on
 *	the domain boundary to a much smaller tolerance than that provided by the mesh generator based on the `unrealistic`
 *	flag. This flag is so named because the tolerance of the geometry based on the CAD or eventual manufacturing would
 *	almost certainly be much larger than machine precision. One is thus running an inherently unrealistic simulation
 *	when this correction is made.
 *
 *	The functionality is used for assessing convergence with extremely small error levels on curved domains for research
 *	cases.
 *
 *	\warning In light of the above discussion, very special consideration should be made before implementing a
 *	         correction of vertex coordinates to the exact boundary.
 *
 *	\section s2_mesh_vert Future Extensions
 *
 *	If functionality is eventually provided for the CAD representation of domain geometry to be read directly, each
 *	vertex should store additional information:
 *	- an index specifying a curved surface on which the vertex is located;
 *	- an array of parametric coordinates relating to the vertex position on the surface.
 */

struct const_Intrusive_List;
struct Mesh;
struct Mesh_Input;

/// \brief Container for additional information relating to the mesh vertices.
struct Mesh_Vertices {
	const struct const_Vector_i*const ve_curved;        ///< Flags for which vertices are located on curved boundaries.
	const struct const_Vector_i*const ve_boundary;      ///< Flags for which vertices are located on boundaries.
	const struct const_Multiarray_Vector_i*const ve_bc; /**< Holds the sorted values of all unique boundary conditions
	                                                     *   associated with the vertices. */
};

/** \brief Constructor for the \ref Mesh_Vertices with correction for exact geometry under special circumstances.
 *	\return Standard.
 */
struct Mesh_Vertices* constructor_Mesh_Vertices
	(const struct Mesh*const mesh,                ///< \ref Mesh.
	 const struct const_Intrusive_List* elements, ///< The base \ref Element list.
	 const struct Mesh_Input*const mesh_input     ///< \ref Mesh_Input.
	);

/// \brief Destructor for \ref Mesh_Vertices.
void destructor_Mesh_Vertices
	(struct Mesh_Vertices* mesh_vert ///< Standard.
	);

#endif // DPG__mesh_vertices_h__INCLUDED
