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
 *  \todo update includes.
 */

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void increment_lhs_boundary_face_T
	(struct Matrix_T* lhs, const struct Numerical_Flux_T* num_flux, const struct Solver_Face* s_face,
	 const struct Simulation* sim)
{
	UNUSED(sim);
	assert(((struct Face*)s_face)->boundary);

	struct Matrix_T* lhs_ll = constructor_lhs_f_1_T((int[]){0,0},num_flux,s_face); // destructed

	set_block_Matrix_T(lhs,(struct const_Matrix_T*)lhs_ll,0,0,'a');
	destructor_Matrix_T(lhs_ll);
}
