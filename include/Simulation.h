// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__Simulation_h__INCLUDED
#define DPG__Simulation_h__INCLUDED
/**	\file
 *	\brief Provides the struct for top-level simulation related information and its associated functions.
 */

/**\{ \name Alternate notation
 *	Provide an alternate notation to be deleted after eventual code refactoring.
 *	\todo Fix the notation.
 */
#define Element S_ELEMENT
#define Volume  S_VOLUME
#define Face    S_FACE
///\}

/// \brief Struct holding data related to the simulation.
struct Simulation {
	const unsigned int d,     ///< Dimension.
	                   n_var, ///< Number of variables in the PDE under consideration.
	                   n_eq;  ///< Number of equations in the PDE under consideration.

	const struct Element*const element_head; ///< Pointer to the head of the Element list.
	      struct Volume*       volume_head;  ///< Pointer to the head of the Volume  list.
	      struct Face*         face_head;    ///< Pointer to the head of the Face    list.
};

/** \brief Constructor for \ref Simulation.
 *
 *	As the struct contains many data members, only the memory allocation is performed as part of the constructor. The
 *	appropriate setter functions must be called to define the members.
 */
struct Simulation* constructor_Simulation ();

/// \brief Destructor for \ref Simulation.
void destructor_Simulation
	(struct Simulation* sim ///< Standard.
	);

/// \brief Set several \ref Simulation parameters.
void set_Simulation_parameters
	(struct Simulation*const sim, ///< Standard.
	 unsigned int d,              ///< Defined in \ref Simulation.
	 unsigned int n_var,          ///< Defined in \ref Simulation.
	 unsigned int n_eq            ///< Defined in \ref Simulation.
	);

/// \brief Set \ref Simulation::element_head.
void set_Simulation_element
	(struct Simulation*const sim,      ///< Standard.
	 const struct Element*const e_head ///< Defined in \ref Simulation.
	);

/// \brief Set \ref Simulation::volume_head.
void set_Simulation_volume
	(struct Simulation*const sim, ///< Standard.
	 struct Volume* v_head        ///< Defined in \ref Simulation.
	);

/// \brief Set \ref Simulation::face_head.
void set_Simulation_face
	(struct Simulation*const sim, ///< Standard.
	 struct Face* f_head          ///< Defined in \ref Simulation.
	);

#endif // DPG__Simulation_h__INCLUDED
