// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__geometry_h__INCLUDED
#define DPG__geometry_h__INCLUDED
/**	\file
 *	\brief Provides the interface to functions used for geometry processing.
 */

struct Simulation;
struct Intrusive_List;
struct Solver_Volume;

/**	\brief Set up the geometry for the simulation. Computes:
 *	- \ref Volume::geom_coef.
 */
void set_up_geometry
	(struct Simulation* sim,        ///< \ref Simulation.
	 struct Intrusive_List* volumes ///< The volumes for which to set up the geometry.
	);

/**	\brief Set up the solver geometry:
 *	- \ref Solver_Volume::metrics_vg;
 *	- \ref Solver_Volume::metrics_vc;
 *	- \ref Solver_Volume::jacobian_det_vc;
 *	- \todo [ref here] Solver_Face::normals_fc;
 *	- \todo [ref here] Solver_Face::jacobian_det_fc;
 *
 *	Requires that:
 *	- \ref Simulation::volumes points to a list of \ref Solver_Volume\*s;
 *	- \ref Simulation::faces   points to a list of \todo [ref here] Solver_Face\*s.
 */
void set_up_solver_geometry
	(struct Simulation* sim ///< \ref Simulation.
	);

#endif // DPG__geometry_h__INCLUDED
