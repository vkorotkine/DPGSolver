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

#ifndef DPG__functionals_h__INCLUDED
#define DPG__functionals_h__INCLUDED

#include <complex.h>

struct Simulation;
struct Intrusive_List;


/** \brief Pointer to objective functions used for the optimization
 *
 *  \param sim    \ref Simulation.
 */

typedef double (*functional_fptr)
	(const struct Simulation* sim ///< Standard. \ref Simulation
	);


/** \brief Pointer to objective functions (complex version) used for the optimization
 *
 *  \param sim    \ref Simulation.
 */

typedef double complex (*functional_fptr_c)
	(const struct Simulation* sim ///< Standard. \ref Simulation
	);


/** \brief A functional to compute the lift coefficient.
 *
 * \return The lift coefficient (Cl) value
 */
double functional_cl(
	const struct Simulation* sim ///< Standard. \ref Simulation
	);


/** \brief Complex version of functional_cl
 *
 * \return The complex value of the lift coefficient (Cl) value
 */
double complex functional_cl_c(
	const struct Simulation* sim_c ///< Standard. \ref Simulation
	);


/** \brief A functional which compute the total volume of the mesh.
 *
 *	\return The total volume of the mesh
 */
double functional_mesh_volume(
	const struct Simulation* sim ///< Standard. \ref Simulation
	);


/** \brief Complex version of functional_mesh_volume.
 *
 *	\return The total volume of the mesh (complex)
 */
double complex functional_mesh_volume_c(
	const struct Simulation* sim ///< Standard. \ref Simulation
	);


/** \brief Computes the functional defined by
 *
 *		I = cm_le
 *
 *	where cm_le is the moment coefficient about the leading edge.
 * 
 * 	\return Value of the functional
 */
double functional_cm_le(
	const struct Simulation* sim ///< Standard. \ref Simulation
	); 


/** \brief Complex version of functional_cm_le. 
 * \todo: Will remove this function after templating is complete.
 * 
 * \return Complex value of the functional
 */
double complex functional_cm_le_c(
	const struct Simulation* sim ///< Standard. \ref Simulation
	); 


/** \brief Computes the functional defined by
 *
 *		I = 0.5 * (cl - cl_target)^2
 *
 *	where Cl is the computed lift coefficient and Cl_target is the target
 *	lift coefficient.
 * 
 * 	\return Value of the functional
 */
double functional_target_cl(
	const struct Simulation* sim ///< Standard. \ref Simulation
	); 


/** \brief Complex version of functional_target_cl. 
 * \todo: Will remove this function after templating is complete.
 * 
 * \return Complex value of the functional
 */
double complex functional_target_cl_c(
	const struct Simulation* sim ///< Standard. \ref Simulation
	); 


/** \brief Computes the functional defined by
 * 		
 *		I = V/V0
 *	
 *	where V is the current mesh volume and V0 is the initial mesh volume.
 *
 * 	NOTE: For now, V0 is computed using the sim->volume_initial instance variable. 
 *		Find a more elegant way to do this later perhaps.
 *
 * \return Value of the functional
 */
double functional_mesh_volume_fractional_change(
	const struct Simulation* sim ///< Standard. \ref Simulation
	);


/** \brief Complex version of functional_mesh_volume_fractional_change
 * \todo: Will remove this function after templating is complete.
 *
 * \return Complex value of the functional
 */
double complex functional_mesh_volume_fractional_change_c(
	const struct Simulation* sim_c ///< Standard. \ref Simulation
	);


#endif // DPG__functionals_h__INCLUDED


