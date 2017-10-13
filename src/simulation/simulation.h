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

///\{ \name Definitions for the available solver methods.
#define METHOD_DG  1
#define METHOD_HDG 2
///\}

///\{ \name Definitions for the available PDEs.
#define PDE_ADVECTION    1
#define PDE_POISSON      2
#define PDE_EULER        3
#define PDE_NAVIERSTOKES 4
///\}

/// \brief Struct holding data related to the simulation.
struct Simulation {
	const int mpi_size, ///< The number of mpi processes.
	          mpi_rank; ///< The mpi rank of the current processor.

	const char* ctrl_name_full; ///< Name of the control file (including full path and file extension).
	const char mesh_name_full[STRLEN_MAX]; ///< Name of the mesh    file (including full path and file extension).
	const char input_path[STRLEN_MAX];     ///< The path to the directory containing relevant input files.

	const char pde_name[STRLEN_MIN];  ///< Name of the Partial Differential Equation (PDE).
	const char pde_spec[STRLEN_MAX];  ///< Additional specifications for the PDE.
	const char geom_name[STRLEN_MAX]; ///< Name of the base geometry to be used for the domain.
	const char geom_spec[STRLEN_MAX]; ///< Additional specifications for the geometry.

	const int d; ///< Dimension.

	/** The type of domain.
	 *  Options:
	 *  	- Straight: All volumes have a geometry order of 1 (affine for simplicies).
	 *  	- Curved:   Volumes along the boundary of the domain may have a higher geometry order.
	 *  	- Mapped:   All volumes are mapped from a simple reference domain to the final simulation domain.
	 *
	 *  For `domain_type = Curved`, it is **necessary** to set the values of the boundary conditions on curved
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
	 *	- superparametric1 (geometry order one greater than the solution order);
	 *	- fixed# (used for geometry which is of fixed order #).
	 */
	const char geom_rep[STRLEN_MIN];

	/** The type of geometry blending to be used. Options:
	 *	- geom_blending_tp:
	 *		- gordon_hall.
	 *	- geom_blending_si:
	 *		- szabo_babuska_gen (Szabo-Babuska generalized);
	 *		- scott;
	 *		- lenoir;
	 *		- nielson;
	 *	- geom_blending_pyr: currently unsupported; requires investigation.
	 */
	const char geom_blending[N_ST_STD][STRLEN_MIN];

	/** The minimal and maximal reference orders. p-adaptation is enabled if the orders differ.
	 *
	 *  \note While this order corresponds to the order of the solution in the domain volumes in the case of the
	 *        standard continuous/discontinuous Galerkin scheme, this is not necessarily the case for other methods.
	 */
	const int p_ref[2];

	/// The minimal and maximal solution orders for volumes. p-adaptation is enabled if the orders differ.
	const int p_s_v[2];

	/// The minimal and maximal solution orders for faces. p-adaptation is enabled if the orders differ.
	const int p_s_f[2];

	/** The minimal and maximal gradient solution orders for volumes. p-adaptation is enabled if the orders
	 *  differ. */
	const int p_sg_v[2];

	/** The minimal and maximal gradient solution orders for faces. p-adaptation is enabled if the orders
	 *  differ. */
	const int p_sg_f[2];

	/** The multiplicative ((times): *) constant of the cubature order in relation to the solution order.
	 *  p_c = p_c_x*p_s + p_c_p. */
	const int p_c_x;

	/** The additive ((p)lus: +) constant of the cubature order in relation to the solution order.
	 *  p_c = p_c_x*p_s + p_c_p. */
	const int p_c_p;

	/** The additive ((p)lus: +) constant of the test function order in relation to the solution order.
	 *  p_t = p_s + p_t_p. */
	const int p_t_p;

	/** The number of required extents for hp operator stored in the various \ref Element\*s. This is currently a
	 *  fixed value determined from the most general hp adaptive case.
	 *  - h-adaptation: Does not add any extents, but increases the range of the first extent such that operators
	 *                  acting between sub-elements can be stored. Unlike the case for the p-adaptive operators, it is
	 *                  quite rare to require an operator from a fine to coarse element and the additional index was
	 *                  thus not added for this case.
	 *  - p-adaptation: Adds two extents such that operators acting between bases of different orders can be stored.
	 *                  The additional extent ordering is [p_in][p_out].
	 *  - hp-adaptation: Add both p and h adaptation indices:
	 *                  The additional extent ordering is [h_range][p_in][p_out].
	 */
	const int n_hp;

	const int adapt_type; ///< The type of adaptation to be used. Set based on the ctrl file paramters.

/// \todo maybe move to a testcase struct.
	const int pde_index; ///< Index corresponding to \ref pde_name.

	const int n_var, ///< Number of variables in the PDE under consideration.
	          n_eq;  ///< Number of equations in the PDE under consideration.

	const struct const_Intrusive_List* elements; ///< Pointer to the head of the Element list.
	struct Intrusive_List* volumes;              ///< Pointer to the head of the Volume  list.
	struct Intrusive_List* faces;                ///< Pointer to the head of the Face    list.

// ---------------------------- //
	const bool collocated; /**< Whether a collocated interpolation and integration node set is being used. Significant
	                        *   performance increase may be observed when this is `true`. */

	const int method;     /**< Solver method to be used.
	                       *   	Options: 1 (DG), 2 (HDG), 3 (HDPG), 4 (DPG). */
};

/** \brief Constructor for \ref Simulation.
 *	\return Standard.
 *
 *	As the struct contains many data members, only the memory allocation is performed as part of the constructor. The
 *	appropriate setter functions must be called to define the members.
 */
struct Simulation* constructor_Simulation
	(const char*const ctrl_name ///< The partial name of the control file.
	);

/// \brief Destructor for \ref Simulation.
void destructor_Simulation
	(struct Simulation* sim ///< Standard.
	);

/** \brief Set full control file name (including path and file extension).
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

#endif // DPG__Simulation_h__INCLUDED
