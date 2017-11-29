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
 *  \brief Provides the common static templated boundary condition functions.
 */

#include "boundary.h"

/** \brief Get the compute_member array for the current boundary value computation.
 *  \return Statically allocated array of flags. */
static bool* get_compute_member
	(const struct Boundary_Value_Input_T* bv_i ///< \ref Boundary_Value_Input_T.
	)
{
	static bool c_m[] = { true, false, false, };

#ifdef TYPE_REAL
	if (1) {
#elif TYPE_COMPLEX
	if (bv_i->has_complex_J) {
#endif
		struct Boundary_Value_Input_R* bv_i_r = (struct Boundary_Value_Input_R*) bv_i;
		for (int i = 1; i < MAX_FLUX_OUT; ++i)
			c_m[i] = bv_i_r->compute_member[i];
	}

	return c_m;
}

/// \brief Constructor for the solution on the boundary.
static const struct const_Multiarray_T* constructor_sol_bv
	(const struct const_Multiarray_d* xyz, ///< xyz coordinates.
	 const struct Simulation* sim          ///< \ref Simulation.
	)
{
	const struct const_Multiarray_d* s = sim->test_case->constructor_sol(xyz,sim); // keep/destructed
#ifdef TYPE_REAL
	return s;
#elif TYPE_COMPLEX
	const struct const_Multiarray_c* s_c = constructor_copy_const_Multiarray_c_Multiarray_d(s); // keep
	destructor_const_Multiarray_d(s);
	return s_c;
#endif
}
