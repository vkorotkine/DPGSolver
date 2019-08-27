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

#ifndef DPG__Simulation_h__INCLUDED
#define DPG__Simulation_h__INCLUDED
/** \file
 *  \brief Provides the interface for the \ref Simulation container and associated functions.
 *
 *  \note Some of the orders specified in the ctrl file will be disregarded if \ref Simulation::collocated is true.
 */

#include <stdbool.h>
#include <stddef.h>

#include "definitions_alloc.h"
#include "definitions_elements.h"

/// \brief Struct holding data related to the simulation.
struct Simulation {
	const int mpi_size, ///< The number of mpi processes.
	          mpi_rank; ///< The mpi rank of the current processor.

	const char*const ctrl_name;            ///< Name of the control file.
	const char* ctrl_name_full;            ///< Name of the control file (including full path and file extension).
	const char mesh_name_full[STRLEN_MAX]; ///< Name of the mesh    file (including full path and file extension).

	const char pde_name[STRLEN_MIN];  ///< Name of the Partial Differential Equation (PDE).
	const char pde_spec[STRLEN_MAX];  ///< Additional specifications for the PDE.
	const char geom_name[STRLEN_MAX]; ///< Name of the base geometry to be used for the domain.
	const char geom_spec[STRLEN_MAX]; ///< Additional specifications for the geometry.

	/** Extension on the file name from which to read the test case data.
	 *  If this variable is not provided in the input control file, the standard test case data file is used:
	 *  "test_case.data". If it is provided, the file name used is "test_case_{test_case_extension}.data".
	 */
	const char test_case_extension[STRLEN_MAX];

	/** Extension on the file name from which to read the solution data.
	 *  If this variable is not provided in the input control file, the standard solution data file is used:
	 *  "solution.data". If it is provided, the file name used is "solution_{solution_extension}.data".
	 */
	const char solution_extension[STRLEN_MAX];

	/** The type of domain.
	 *  Options:
	 *  	- Straight: All volumes have a geometry order of 1 (affine for simplicies).
	 *  	- Blended:  Volumes along the boundary of the domain may have a higher geometry order.
	 *  	- Mapped:   All volumes are mapped from a simple reference domain to the final simulation domain.
	 *
	 *  For `domain_type = Blended`, it is **necessary** to set the values of the boundary conditions on curved
	 *  surfaces to different multiples of BC_STEP_SC+bc_value for the `boundary` and `curved` flags for
	 *  \ref Volume and \ref Face elements to be properly set.
	 */
	const int domain_type;

	/** The minimal and maximal mesh levels to be used. h-adaptation is enabled if the levels differ. The input mesh
	 *  is chosen based on the value of ml[0]. */
	const int ml[2];

	/** Flag for whether the mesh vertices should be unrealistically corrected to lie on the input domain boundary
	 *  to within a very small tolerance. See \ref mesh_vertices.h additional discussion of this issue. */
	const bool mesh_unrealistic;

	ptrdiff_t n_v, ///< The number of \ref Volume finite elements.
	          n_f; ///< The number of \ref Face   finite elements.


	/** The type of nodes used for polynomial interpolation. Options:
	 *  - interp_tp:  GL, GLL, EQ;
	 *  - interp_si:  AO, WSH, EQ;
	 *  - interp_pyr: GL, GLL.
	 *
	 *  See \ref definitions_nodes.h for the expansion of the acronyms.
	 *
	 *  \warning EQ nodes should only be used for comparison with other codes where GL/GLL nodes are not available as
	 *           they provide an increasingly poor interpolation as the order is increased due to growth of their
	 *           Lebesgue constant.
	 *
	 *  In the case of a collocated scheme, these nodes are also used for the cubature.
	 */
	const char nodes_interp[N_ST_STD][STRLEN_MIN];

	/** The type of basis used for the geometry representation. Options:
	 *	- lagrange;
	 *	- bezier;
	 *	- nurbs.
	 *  see \ref check_necessary_simulation_parameters for supported options.
	 */
	const char basis_geom[STRLEN_MIN];

	/** The type of basis used for the solution representation. Options:
	 *	- orthonormal;
	 *	- lagrange;
	 *	- bezier.
	 */
	const char basis_sol[STRLEN_MIN];

	/** The type of geometry representation used. Options:
	 *	- isoparametric;
	 *	- superparametric        (geometry order one greater than the solution order);
	 *	- superparametric_p_le_1 (geometry order one greater than the solution order for p <= 1 (iso otherwise));
	 *	- fixed# (used for geometry which is of fixed order #).
	 */
	const char geom_rep[STRLEN_MIN];

	/** The type of geometry blending to be used. Options:
	 *	- geom_blending_tp:
	 *		- gordon_hall.
	 *	- geom_blending_si:
	 *		- szabo_babuska_gen (Szabo-Babuska generalized);
	 *		- scott;
	 *		- lenoir (eq. 22, \cite Lenoir1986);
	 *		- lenoir_simple (Excludes intermediate corrections (remark following eq. 22, \cite Lenoir1986);
	 *		- nielson;
	 *	- geom_blending_pyr: currently unsupported; requires investigation.
	 */
	const char geom_blending[N_ST_STD][STRLEN_MIN];

	/** Whether or not NURBS enhanced metrics should be used for the case. 
	 * 	If true, the geometry file provided must contain NURBS information for the
	 *	mesh. The NURBS patch will then be used to compute all the metric terms.
	 * 	By default, if no parameter is set in the input file, this is false.
	 */
	const bool nurbs_enhanced_metrics;

	/** Whether or not use the multipatch mapping. 
	 * If true, the geometry file provided must contain information for 
	 * multiple NURBS patches, as generated by 
	 * input/meshes/NURBS_Patch_Generator/NURBS_Parametric_Domain_Generator_Multiple.py
	 */
	const bool nurbs_multipatch;

	/**
	 * The number of patches to use in the NURBS multipatch. Goes together with 
	 * the nurbs_multipatch bool. 
	 */
	const int nurbs_n_patches; 

	/** The minimal and maximal reference orders. p-adaptation is enabled if the orders differ.
	 *
	 *  \note While this order corresponds to the order of the solution in the domain volumes in the case of the
	 *        standard continuous/discontinuous Galerkin scheme, this is not necessarily the case for other methods.
	 */
	const int p_ref[2];

	/** The minimal and maximal isoparametric geometry orders.
	 *
	 *  These orders differ from \ref Simulation::p_ref in the case where it is desired to use straight internal faces
	 *  for a parametric mesh in the case of the solution order not including p = 1.
	 */
	const int p_ig[2];

	const int p_s_v_p,  ///< The additive ((p)lus: +) constant of the volume solution order relative to p_ref.
	          p_s_f_p,  ///< The additive ((p)lus: +) constant of the face   solution order relative to p_ref.
	          p_sg_v_p, ///< The additive ((p)lus: +) constant of the volume solution gradient order relative to p_ref.
	          p_sg_f_p; ///< The additive ((p)lus: +) constant of the face   solution gradient order relative to p_ref.

	/** The multiplicative ((times): *) constant of the cubature order in relation to the solution order.
	 *  p_c = p_c_x*p_s + p_c_p. */
	const int p_c_x[2];

	/** The additive ((p)lus: +) constant of the cubature order in relation to the solution order.
	 *  p_c = p_c_x*p_s + p_c_p. */
	const int p_c_p[2];

	/** The additive ((p)lus: +) constant of the test function order in relation to the solution order.
	 *  p_t = p_s + p_t_p. */
	const int p_t_p[2];

	/// The minimal and maximal solution orders for volumes. p-adaptation is enabled if the orders differ.
	const int p_s_v[2];

	/// The minimal and maximal solution/flux orders for faces. p-adaptation is enabled if the orders differ.
	const int p_s_f[2];

	/** The minimal and maximal gradient solution orders for volumes. p-adaptation is enabled if the orders
	 *  differ. */
	const int p_sg_v[2];

	/** The minimal and maximal gradient/trace solution orders for faces. p-adaptation is enabled if the orders
	 *  differ. */
	const int p_sg_f[2];

	const int ml_p_curr[2]; ///< The current reference 'm'esh 'l'evel and 'p'olynomial order of the simulation.

	/// Finite element method to be used. Options: 1 (DG), 2 (HDG), 3 (HDPG), 4 (DPG).
	const int method;

	const bool collocated; /**< Whether a collocated interpolation and integration node set is being used.
	                        *   Significant performance increase may be observed when this is `true`. */

	/** The number of required extents for hp operator stored in the various \ref Element\*s. This is currently a
	 *  fixed value determined from the most general hp adaptive case.
	 *  - h-adaptation: Does not add any extents, but increases the range of the first extent such that operators
	 *                  acting between sub-elements can be stored. Unlike the case for the p-adaptive operators, it
	 *                  is quite rare to require an operator from a fine to coarse element and the additional index
	 *                  was thus not added for this case.
	 *  - p-adaptation: Adds two extents such that operators acting between bases of different orders can be stored.
	 *                  The additional extent ordering is [p_in][p_out].
	 *  - hp-adaptation: Add both p and h adaptation indices:
	 *                  The additional extent ordering is [h_range][p_in][p_out].
	 */
	const int n_hp;

	const int adapt_type; ///< The type of adaptation to be used. Set based on the ctrl file paramters.

	struct Test_Case_rc* test_case_rc; ///< Pointer to the \ref Test_Case_rc.

	const struct const_Intrusive_List* elements; ///< Pointer to the head of the Element list.
	struct Intrusive_List* volumes;              ///< Pointer to the head of the Volume  list.
	struct Intrusive_List* faces;                ///< Pointer to the head of the Face    list.

	struct Geo_Data* geometric_data; 			 ///< Pointer to the geometric data struct.

	double mesh_volume_initial; ///< The volume of the initial mesh. \todo Look to place this somewhere with the functionals

};

/** \brief Constructor for \ref Simulation omitting the construction of members dependent upon an input mesh.
 *  \return Standard. */
struct Simulation* constructor_Simulation__no_mesh
	(const char*const ctrl_name ///< The partial name of the control file.
	);

/** \brief Constructor for \ref Simulation.
 *  \return Standard. */
struct Simulation* constructor_Simulation
	(const char*const ctrl_name ///< The partial name of the control file.
	);

/** \brief Constructor for a \ref Simulation with the minimal required information from a restart file.
 *  \return See brief. */
struct Simulation* constructor_Simulation_restart
	(const struct Simulation*const sim_main ///< The main \ref Simulation.
	);

/// \brief Destructor for \ref Simulation.
void destructor_Simulation
	(struct Simulation* sim ///< Standard.
	);

/** \brief Set full control file name (including path and file extension).
 *  \return See brief.
 *
 *  If "TEST" is not included as part of the name (default option):
 *  - it is assumed that `ctrl_name` includes the full name and path from CMAKE_PROJECT_DIR/control_files;
 *  - `ctrl_name_full` -> ../control_files/`ctrl_name`.
 *  otherwise:
 *  - it is assumed that `ctrl_name` includes the full name and path from CMAKE_PROJECT_DIR/testing/control_files;
 *  - `ctrl_name_full` -> ../testing/control_files/`ctrl_name`.
 */
const char* set_ctrl_name_full
	(const char*const ctrl_name ///< Defined in \ref set_simulation_core.
	);

/**	\brief Set the \ref Mesh_Input based on the parameters in the \ref Simulation.
 *	\return See brief. */
struct Mesh_Input set_Mesh_Input
	(const struct Simulation*const sim ///< \ref Simulation.
	);

/// \brief Set \ref Simulation::elements.
void set_Simulation_elements
	(struct Simulation*const sim,          ///< Standard.
	 struct const_Intrusive_List* elements ///< See \ref Simulation.
	);

/** \brief Computes the value of \ref Simulation::adapt_type based on the input order/mesh level parameters.
 *  \return See brief. */
int compute_adapt_type
	(const int p_ref[2], ///< The array of minimal and maximal orders.
	 const int ml[2]     ///< The array of minimal and maximal mesh levels.
	);

/// \brief Set \ref Simulation::ml_p_curr.
void set_ml_p_curr
	(const int ml,          ///< The mesh level.
	 const int p,           ///< The polynomial order.
	 struct Simulation* sim ///< \ref Simulation.
	);

/** \brief Return the value of the statically allocated operator format to be used.
 *  \return See brief.
 *
 *  Passing a non-zero value for `new_format` sets the statically allocated operator format to that value.
 */
char get_set_op_format
	(const char new_format ///< New format. Options: 'd'efault, 's'tandard, 't'ensor-product, 'c'ompressed sparse row.
	);

/** \brief Return a statically allocated `bool` flag indicating whether collocated interpolation and cubature nodes are
 *         being used.
 *  \return See brief.
 *
 *  Passing a non-NULL value for `new_val` sets the statically allocated value to that pointed to by the input.
 */
bool get_set_collocated
	(const bool*const new_val ///< Pointer to new value.
	);

/** \brief Return a statically allocated `int` referring to the solver method being used.
 *  \return See brief.
 *
 *  Passing a non-NULL input for `new_val` sets the statically value to that of the input.
 */
int get_set_method
	(const int*const new_val ///< The new value if non-NULL.
		);

/** \brief Return a statically allocated `int` referring to the domain type being used.
 *  \return See brief.
 *
 *  Passing a non-NULL input for `new_val` sets the statically value to that of the input.
 */
int get_set_domain_type
	(const int*const new_val ///< The new value if non-NULL.
		);

/** \brief Return a statically allocated `int*` for the polynomial degree of the specified basis.
 *  \return See brief; the two indices of the returned array are for straight ([0]) and curved ([1]) values.
 *
 *  Passing a non-NULL input for `new_val` sets the statically value to that of the input.
 */
const int* get_set_degree_poly
	(const int*const new_vals, ///< The new values if non-NULL.
	 const char*const key      /**< The key used to obtain specify the desired variable. See the function definition
				    *   for options; the notation is similar to that of \ref element_operators.h */
	 );

#endif // DPG__Simulation_h__INCLUDED
