// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__geometry_h__INCLUDED
#define DPG__geometry_h__INCLUDED
/**	\file
 *	\brief Provides the interface to functions used for geometry processing.
 */

struct Simulation;
struct Intrusive_List;

/** \brief Set up the geometry for the simulation. Computes:
 *	- \ref Volume::geom_coef.
 */
void set_up_geometry
	(struct Simulation* sim,        ///< \ref Simulation.
	 struct Intrusive_List* volumes ///< The volumes for which to set up the geometry.
	);

/** \brief Set up the geometry for the solve. Computes:
 *	- \ref Solver_Volume::metrics_vg;
 *	- \ref Solver_Volume::metrics_vc;
 *	- \todo [ref here] Solver_Face::metrics_fc;
 */
void set_up_geometry_solver
	(struct Simulation* sim,        ///< \ref Simulation.
	 struct Intrusive_List* volumes ///< The volumes for which to set up the solver geometry.
	);

#endif // DPG__geometry_h__INCLUDED
