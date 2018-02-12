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
 *  \brief Provides the interface to functions used for geometry processing.
 */

struct Simulation;
struct Intrusive_List;
struct Solver_Volume_T;
struct Solver_Face_T;
struct Multiarray_R;
struct const_Multiarray_R;

/** \brief Pointer to functions constructing the physical xyz coordinates from parametric/straight element coordinates.
 *
 *  \param n_type The geometry node type. Options: volume curved 'g'eometry, p2 'v'ertices.
 *  \param xyz_i  The input xyz coordinates.
 *  \param s_vol  \ref Solver_Volume_T.
 *  \param sim    \ref Simulation.
 */
typedef const struct const_Multiarray_R* (*constructor_xyz_fptr_T)
	(const char n_type,
	 const struct const_Multiarray_R* xyz_i,
	 const struct Solver_Volume_T* s_vol,
	 const struct Simulation* sim
	);

/// \brief Set up the solver geometry for all volumes and faces.
void set_up_solver_geometry_T
	(struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Compute the face unit normal vectors at the nodes corresponding to the given face metrics.
void compute_unit_normals_T
	(const int ind_lf,                             ///< Defined for \ref compute_unit_normals_and_det_T.
	 const struct const_Multiarray_R* normals_ref, ///< Defined for \ref compute_unit_normals_and_det_T.
	 const struct const_Multiarray_R* metrics_f,   ///< Defined for \ref compute_unit_normals_and_det_T.
	 struct Multiarray_R* normals_f                ///< Defined for \ref compute_unit_normals_and_det_T.
	);

/** \brief Compute the geometry of the \ref Solver_Volume_T.
 *
 *  The following members are set:
 *  - Solver_Volume_T::metrics_vm;
 *  - Solver_Volume_T::metrics_vc;
 *  - Solver_Volume_T::jacobian_det_vc.
 *
 *  Following the analysis of Kopriva \cite Kopriva2006, the metric terms are computed using the curl-form such that
 *  the free-stream preservation property may be recovered. The consistent symmetric-conservative (CSC) metric of Abe
 *  and Haga (section 5.3) is used for the implementation \cite Abe2015. The steps are repeated below to clarify the
 *  equivalence of the prodecure adopted here with their procedure:
 *  - (step 0-1) As the metric contributions computed in step 0 are computed in a basis of sufficient order to
 *    represent them exactly and are subsequently interpolated to the consistent grid points (CGPs), the metric
 *    contributions are here computed directly at the CGPs.
 *  	- Our terminology for the GPs is R_vg ((R)eference coordinates of the (v)olume (g)eometry nodes).
 *  	- Our terminology for the CGPs is R_vm ((R)eference coordinates of the (v)olume (m)etric nodes).
 *  	- We allow for flexibility in the order of the R_vm nodes such that superparametric geometry can be used on
 *  	  curved domain boundaries; Abe and Haga use an isoparametric partial metric representation **before** the
 *  	  differentiation is applied, resulting in a subparametric metric representation (see eq. (43) \cite Abe2015).
 *  - (step 2) The computed metric terms are interpolated to the solution points (SPs).
 *  	- As the flux reconstruction scheme is collocated (solution interpolation and cubature nodes are coincident),
 *  	  the interpolation to the SPs is equivalent to the interpolation to the cubature nodes. Thus, interpolation
 *  	  to the R_vc ((R)eference coordinates of the (v)olume (c)ubature) is then performed in the implementation
 *  	  here.
 *
 *  \todo Investigate requirement of superparametric geometry on curved surfaces and add comments. Potentially ok by
 *        using over-integration in curved elements.
 *
 *  Given the 3D geometry Jacobian ordering of
 *
 *  \f{eqnarray*}{
 *  	J  = \{ &\{x_r,x_s,x_t\}, &\\
 *  	        &\{y_r,y_s,y_t\}, &\\
 *  	        &\{z_r,z_s,z_t\}  &\},
 *  \f}
 *
 *  using the nonconservative metric (NC) for clarity of exposition (section 5.1 \cite Abe2015), the ordering of the
 *  metric terms is:
 *
 *  \f{eqnarray*}{
 *  	m  = \{ &\{ +(y_s z_t - y_t z_s), -(y_r z_t - y_t z_r), +(y_r z_s - y_s z_r) \}, &\\
 *  	        &\{ -(x_s z_t - x_t z_s), +(x_r z_t - x_t z_r), -(x_r z_s - x_s z_r) \}, &\\
 *  	        &\{ +(x_s y_t - x_t y_s), -(x_r y_t - x_t y_r), +(x_r y_s - x_s y_r) \}, &\}.
 *  \f}
 */
void compute_geometry_volume_T
	(struct Solver_Volume_T* s_vol, ///< \ref Solver_Volume_T.
	 const struct Simulation* sim   ///< \ref Simulation.
	);

/** \brief Compute the geometry of the \ref Solver_Face_T.
 *
 *  The following members are set:
 *  - Solver_Face_T::xyz_fc;
 *  - Solver_Face_T::n_fc;
 *  - Solver_Face_T::jacobian_det_fc.
 */
void compute_geometry_face_T
	(struct Solver_Face_T* s_face, ///< \ref Solver_Face_T.
	 struct Simulation* sim        ///< \ref Simulation.
	);
