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
 *  \brief Provides the interface for the templated \ref Solver_Volume_T container and associated functions.
 */

#include <stddef.h>

/** \brief Container for data relating to the solver volumes.
 *  \note Certain members are declared `const` despite requiring modification for adaptive simulations. Only members
 *        changing with every solver iteration are `mutable`.
 */
struct Solver_Volume_T {
	struct Volume volume; ///< The base \ref Volume.

	/// The index of the first degree of freedom (dof) of the volume in relation to the global dof.
	const ptrdiff_t ind_dof;

	/** The index of the first degree of freedom (dof) of the constraint-related volume members in relation to the
	 *  global dof. */
	const ptrdiff_t ind_dof_constraint;

	/// The index of the first degree of freedom (dof) of test functions for the volume.
	const ptrdiff_t ind_dof_test;

	/// The reference order of the volume. Need not be equal to the order of the solution in the volume.
	const int p_ref;

	const int ml; ///< The mesh level of the volume.

	/** The geometry coefficients of the volume in the \ref Simulation::basis_geom. For each of the supported
	 *  \ref Simulation::domain_type options, geom_coef represents:
	 *  - DOM_STRAIGHT: the projection of xyz_ve into the geometry basis of order 1.
	 *  - DOM_BLENDED:
	 *  	- straight volumes: [See DOM_STRAIGHT];
	 *  	- boundary volumes: the coefficients of the *blended* face geometry of order k_g.
	 *  - DOM_PARAMETRIC: the coefficients of the mapped geometry of order k_g.
	 */
	const struct const_Multiarray_T*const geom_coef;

	/** The geometry coefficients corresponding to \ref Volume::xyz_ve (degree: p = 1). For each of the
	 *  supported \ref Simulation::domain_type options, geom_coef represents:
	 *  - DOM_STRAIGHT: the interpolation of xyz_ve into the geometry basis of order p_g.
	 *  - DOM_PARAMETRIC: the high-order geometry __before__ being mapped through the parametrization.
	 */
	const struct const_Multiarray_T*const geom_coef_p1;

	/// Pointer to function used to construct the parametrized surface geometry values (Set to NULL if not requried).
	constructor_xyz_surface_fptr_T constructor_xyz_surface;

	/// The coefficients of the solution in the \ref Simulation::basis_sol.
	struct Multiarray_T* sol_coef;

	/// The coefficients of the solution gradient in the \ref Simulation::basis_sol.
	struct Multiarray_T* grad_coef;

	/** The metric terms (cofactors of the geometry Jacobian) used for transformation of integrals between physical
	 *  and computational space stored at the (v)olume (m)etric nodes. */
	const struct const_Multiarray_T*const metrics_vm;

	/** The metric terms (cofactors of the geometry Jacobian) used for transformation of integrals between physical
	 *  and computational space stored at the (v)olume (c)ubature nodes. */
	const struct const_Multiarray_T*const metrics_vc;

	/// The determinant of the geometry mapping Jacobian evaluated at the volume cubature nodes.
	const struct const_Multiarray_T*const jacobian_det_vc;

	/// Same as \ref Solver_Volume_T::metrics_vm but for the p1 geometry coefficients.
	const struct const_Multiarray_T*const metrics_vm_p1;

	struct Vector_T* flux_imbalance; ///< The values of the flux imbalances for each equation.
	struct Multiarray_T* l_mult; ///< The values of the Lagrange multipliers used to enforce conservation.

	struct Multiarray_T* rhs;   ///< The rhs terms.
	struct Multiarray_T* rhs_0; ///< The rhs terms computed from the initial solution.

	/** Optimal test function coefficients for the solution. These are included here and not in derived Solver Volumes
	 *  so that they may be used outside of the solver functions (visualization, comparison, etc). */
	struct Multiarray_T* test_s_coef;
};

/// \brief Constructor for a derived \ref Solver_Volume_T.
void constructor_derived_Solver_Volume_T
	(struct Volume* volume_ptr,   ///< Pointer to the volume.
	 const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Destructor for a derived \ref Solver_Volume_T.
void destructor_derived_Solver_Volume_T
	(struct Volume* volume_ptr ///< Pointer to the volume.
	);

/** \brief Get the appropriate sub-range of the \ref Solver_Element::w_vc operators.
 *  \return See brief. */
const struct const_Vector_d* get_operator__w_vc__s_e_T
	(const struct Solver_Volume_T* s_vol ///< The current volume.
	);

/** \brief Constructor for the mass matrix of the input volume.
 *  \return See brief. */
const struct const_Matrix_R* constructor_mass_T
	(const struct Solver_Volume_T*const s_vol ///< Standard.
	);

/** \brief Constructor for the inverse mass matrix of the input volume.
 *  \return See brief. */
const struct const_Matrix_R* constructor_inverse_mass_T
	(const struct Solver_Volume_T*const s_vol, ///< Standard.
	 const struct const_Matrix_R*const mass    ///< Mass matrix. Input if available, otherwise pass `NULL`.
	);
