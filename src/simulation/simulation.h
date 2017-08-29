// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__Simulation_h__INCLUDED
#define DPG__Simulation_h__INCLUDED
/**	\file
 *	\brief Provides the struct for top-level simulation related information and its associated functions.
 */

#include <stdbool.h>
#include <stddef.h>

#include "constants_alloc.h"

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
	const char ctrl_name_full[STRLEN_MAX]; ///< Name of the control file (including the full path and file extension).
	const char mesh_name_full[STRLEN_MAX]; ///< Name of the mesh    file (including the full path and file extension).
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
	 *	For `domain_type = Curved`, it is **necessary** to set the values of the boundary conditions on curved surfaces
	 *	to different multiples of BC_STEP_SC+bc_value for the `boundary` and `curved` flags for \ref Volume and \ref
	 *	Face elements to be properly set.
	 */
	const int domain_type;

	const bool unrealistic; /**< Flag for whether the mesh vertices should be unrealistically corrected to lie on the
	                          *  input domain boundary to within a very small tolerance. See \ref mesh_vertices.h
	                          *  additional discussion of this issue. */

	ptrdiff_t n_v, ///< The number of \ref Volume finite elements.
	          n_f; ///< The number of \ref Face   finite elements.



	const char node_type[STRLEN_MIN];  /**< Type of nodes to be used for interpolation.
	                                    *   The node_type input should be of the form (1)_(2) where (1) and (2) denote
	                                    *   the node type to be used for tensor-product and simplex elements,
	                                    *   respectively.
	                                    *   	- tensor-product: GL (Gauss-Legendre), GLL (Gauss-Lobatto-Legendre), EQ
	                                    *   	  (Equally spaced)
	                                    *   	- simplex: AO (Alpha-optimized), WHS (Williams-Ham-Shunn), EQ (Equally
	                                    *   	  spaced)
	                                    */
	const char basis_type[STRLEN_MIN]; /**< Type of basis to be used for solution representation.
	                                    *   	Options: Nodal (Lagrange), Modal (Orthonormal). */

	const bool vectorized, /**< Whether vectorization is being used. When this is enabled, memory for certain variables
	                        *   is allocated contiguously such that fewer blas3 calls on larger arrays are made. This
	                        *   can result in significant performance increase for fixed order, fixed mesh level runs,
	                        *   but advantages may be lost for adapative runs. For this reason, this functionality is
	                        *   currently not supported.
	                        *
	                        *	\todo Investigate further and modify comments above.
	                        */
	           collocated; /**< Whether a collocated interpolation and integration node set is being used. Significant
	                        *   performance increase may be observed when this is `true`. */

	const int method,     /**< Solver method to be used.
	                       *   	Options: 1 (DG), 2 (HDG), 3 (HDPG), 4 (DPG). */
	          adapt_type, /**< Adaptation type. This can be any combination of polynomial (p) or mesh (h) adaptation.
	                       *   	Options: 0 (None), 1 (p), 2 (h), 3 (hp). */
	          p,      ///< Polynomial order to be used for the simulation when p adaptation is disabled.
	          ml,     ///< Mesh level to be used for the simulation when h adaptation is disabled.
	          p_max,  ///< Maximum polynomial order to be used for the simulation when p adaptation is enabled.
	          ml_max; ///< Maximum mesh level to be used for the simulation when h adaptation is enabled.

	const struct const_Intrusive_List*const elements; ///< Pointer to the head of the Element list.
	struct Intrusive_List* volumes;                   ///< Pointer to the head of the Volume  list.
	struct Intrusive_List* faces;                     ///< Pointer to the head of the Face    list.

// ToBeMoved to the solver context.
const int pde_index; ///< Index corresponding to \ref pde_name.

	const int n_var, ///< Number of variables in the PDE under consideration.
	          n_eq;  ///< Number of equations in the PDE under consideration.
};

/** \brief Constructor for \ref Simulation.
 *	\return Standard.
 *
 *	As the struct contains many data members, only the memory allocation is performed as part of the constructor. The
 *	appropriate setter functions must be called to define the members.
 */
struct Simulation* constructor_Simulation ();

/// \brief Destructor for \ref Simulation.
void destructor_Simulation
	(struct Simulation* sim ///< Standard.
	);

/** \brief Set core parameters for the simulation as specified in the control file.
 *
 * 	Requires `sim` to be dynamically allocated. This allows for the definition of `const` members after the declaration
 *	which would otherwise be undefined behaviour.
 */
void set_simulation_core
	(struct Simulation*const sim, ///< Standard.
	 const char*const ctrl_name   ///< Control file name (excluding the file extension).
	);

/// \brief Set several \ref Simulation flags.
void set_Simulation_flags
	(struct Simulation*const sim, ///< Standard.
	 const bool collocated        ///< See \ref Simulation.
	);

/// \brief Set several \ref Simulation parameters.
void set_Simulation_parameters
	(struct Simulation*const sim, ///< Standard.
	 const int d,                 ///< See \ref Simulation.
	 const int n_var,             ///< See \ref Simulation.
	 const int n_eq               ///< See \ref Simulation.
	);

/// \brief Set \ref Simulation::elements.
void set_Simulation_elements
	(struct Simulation*const sim,          ///< Standard.
	 struct const_Intrusive_List* elements ///< See \ref Simulation.
	);

#endif // DPG__Simulation_h__INCLUDED
