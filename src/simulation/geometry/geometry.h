// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__geometry_h__INCLUDED
#define DPG__geometry_h__INCLUDED
/**	\file
 *	\brief Provides the interface to functions used for geometry processing.
 */

struct Simulation;

/** \brief Set up the geometry for the simulation. Computes:
 *	- \ref Volume::geom_coef.
 */
void set_up_geometry
	(struct Simulation* sim ///< \ref Simulation.
	);


#endif // DPG__geometry_h__INCLUDED
