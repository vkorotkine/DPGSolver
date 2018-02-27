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

#ifndef DPG__adaptation_h__INCLUDED
#define DPG__adaptation_h__INCLUDED
/** \file
 *  \brief Provides the interface to functions used to perform hp adaption.
 */

struct Simulation;

/// \brief Container for data relating to element adaptation in the domain.
struct Adaptation_Data {
	/// Coordinates of vertices for which adjacent elements should be refined.
	const struct const_Multiarray_d* xyz_ve_refine;

	const struct const_Vector_i* xyz_ve_ml; ///< Mesh levels to which xyz_ve_refine should be refined.
};

/// \brief Adapt the computational elements using the input adaptation strategy.
void adapt_hp
	(struct Simulation* sim,                       ///< \ref Simulation.
	 const int adapt_strategy,                     /**< The adaptation strategy. Options:
	                                                *   see \ref definitions_adaptation.h. */
	 const struct Adaptation_Data*const adapt_data ///< \ref Adaptation_Data.
	);

#endif // DPG__adaptation_h__INCLUDED
