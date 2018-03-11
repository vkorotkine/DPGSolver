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
/** \file
 *  \brief Provides the interface to templated functions used for blended geometry processing.
 */

#include <stdbool.h>

struct Solver_Volume_T;
struct Simulation;

struct Blended_Parametric_Data_T;

/// \brief Container for boundary computational element related data.
struct Boundary_Comp_Elem_Data_T {
	const struct const_Multiarray_R* xyz_ve; ///< \ref Volume::xyz_ve.

	int s_type; ///< \ref Element::s_type.

	const struct const_Vector_i* bc_boundaries;   ///< \ref Volume::bc_edges or \ref Volume::bc_faces.
	const struct const_Multiarray_Vector_i* b_ve; ///< \ref Element::e_ve or \ref Element::f_ve.

	/** See notation in \ref element_operators.h. The unknown parameter (X) may be set as:
	 *  - 'v'ertex (p2 vertices);
	 *  - 'g'eometry (standard blending).
	 */
	const struct Operator* vv0_vv_vX;

	const struct Multiarray_Operator* vv0_bv_vX; ///< See \ref Boundary_Comp_Elem_Data_T::vv0_vv_vX.

	const struct Operator* vv0_vX_bX;   ///< See \ref Boundary_Comp_Elem_Data_T::vv0_vv_vX.
	const struct Operator* vv0_bX_vX;   ///< See \ref Boundary_Comp_Elem_Data_T::vv0_vv_vX.
	const struct Operator* vv0_vv_bX;   ///< See \ref Boundary_Comp_Elem_Data_T::vv0_vv_vX.
	const struct Operator* vv0_vv_bv;   ///< See \ref Boundary_Comp_Elem_Data_T::vv0_vv_vX.
	const struct Operator* vv0_vv_fcc;  ///< See notation in \ref element_operators.h.
	const struct Operator* vv0_bv_bX;   ///< See \ref Boundary_Comp_Elem_Data_T::vv0_vv_vX.
	const struct Operator* cv0_vgc_bgc; ///< See notation in \ref element_operators.h.
};

/** \brief Container for data required to compute curved surface geometry values for various methods of surface
 *         parametrization.
 */
struct Blended_Parametric_Data_T {
	const int geom_parametrization; ///< \ref Test_Case_T::geom_parametrization.

	const struct const_Multiarray_R* xyz_ve; ///< \ref Volume::xyz_ve.

	const struct Operator* vv0_vv_bX;  ///< See \ref Boundary_Comp_Elem_Data_T.
	const struct Operator* vv0_vv_bv;  ///< See \ref Boundary_Comp_Elem_Data_T.
	const struct Operator* vv0_vv_fcc; ///< See \ref Boundary_Comp_Elem_Data_T.
	const struct Operator* vv0_bv_bX;  ///< See \ref Boundary_Comp_Elem_Data_T.

	const struct const_Vector_R* normal; ///< Outward pointing unit normal vector.

	const char n_type; ///< The 'n'ode "type" of the surface nodes. Options: 'g'eometry, 'v'ertex, 'c'ubature.
	const int domain_type; ///< \ref Simulation::domain_type.

	constructor_xyz_fptr_T constructor_xyz; ///< \ref Test_Case_T::constructor_xyz.
};


/** \brief Version of \ref constructor_xyz_fptr_T for the blended curved volume geometry.
 *  \return See brief. */
const struct const_Multiarray_R* constructor_xyz_blended_T
	(const char n_type,                      ///< See brief.
	 const struct const_Multiarray_R* xyz_i, ///< See brief.
	 const struct Solver_Volume_T* s_vol,    ///< See brief.
	 const struct Simulation* sim            ///< See brief.
	);

/** \brief Constructor for a statically allocated \ref Boundary_Comp_Elem_Data_T containers.
 *  \return See brief. */
struct Boundary_Comp_Elem_Data_T constructor_static_Boundary_Comp_Elem_Data_T
	(const char ce_type,                      ///< Computational element type. Options: 'e'dge, 'f'ace.
	 const char n_type,                       ///< Defined for \ref constructor_xyz_fptr_T.
	 const int p_geom,                        ///< The order of the geometry node basis.
	 const struct Solver_Volume_T*const s_vol ///< The current volume.
	);

/// \brief Destructor for a statically allocated \ref Boundary_Comp_Elem_Data_T containers.
void destructor_static_Boundary_Comp_Elem_Data_T
	(struct Boundary_Comp_Elem_Data_T*const b_ce_d ///< Standard.
	);

/** \brief Set the operator members of \ref Boundary_Comp_Elem_Data_T which depend on the the current order, `p`, or
 *         boundary computational element index, `b`. */
void set_Boundary_Comp_Elem_operators_T
	(struct Boundary_Comp_Elem_Data_T*const b_ce_d, ///< \ref Boundary_Comp_Elem_Data_T.
	 const struct Solver_Volume_T*const s_vol,      ///< The current volume.
	 const char ce_type,                            ///< The type of computational element.
	 const char n_type,                             ///< Defined for \ref constructor_xyz_fptr_T.
	 const int p,                                   ///< The current geometry correction order.
	 const int b                                    ///< The index of the boundary comp. element under consideration.
	);

/** \brief Constructor for a \ref const_Matrix_T\* storing the difference between the surface xyz values for the current
 *         and exact boundary.
 *  \return See brief. */
const struct const_Matrix_R* constructor_xyz_surf_diff_T
	(const struct Boundary_Comp_Elem_Data_T*const b_ce_d, ///< \ref Boundary_Comp_Elem_Data_T.
	 const struct const_Matrix_R*const xyz_i,             ///< Input xyz coordinates.
	 const struct Solver_Volume_T*const s_vol,            ///< The current \ref Solver_Volume_T.
	 const char n_type,                                   ///< \ref Blended_Parametric_Data_T::n_type.
	 const bool use_existing_as_surf,                     /**< Flag for whether the existing high-order geometry
	                                                       *   representation should be used for surface
							       *   representation. */
	 const struct Simulation*const sim                    ///< \ref Simulation.
	);

/** \brief Correct the internal coefficients for the volume geometry representation using the existing face geometry
 *         representation as the surface values. */
void correct_internal_xyz_blended_T
	(struct Solver_Volume_T*const s_vol, ///< \ref Solver_Volume_T.
	 const struct Simulation*const sim   ///< \ref Simulation.
	);
