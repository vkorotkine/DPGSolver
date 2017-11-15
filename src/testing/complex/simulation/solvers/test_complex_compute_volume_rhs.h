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

#ifndef DPG__test_complex_compute_volume_rhs_h__INCLUDED
#define DPG__test_complex_compute_volume_rhs_h__INCLUDED
/** \file
 *  \brief Provides `complex` versions of functions defined in \ref compute_volume_rlhs.h.
 */

struct Flux_Input_c;
struct Simulation;
struct Solver_Volume;

/** \brief `complex` version of \ref constructor_sol_vc_fptr.
 *  \return Standard.
 *
 *  \param s_vol The current volume.
 *  \param sim   \ref Simulation.
 */
typedef const struct const_Multiarray_c* (*constructor_sol_vc_c_fptr)
	(const struct Solver_Volume* s_vol,
	 const struct Simulation* sim
	);

/// \brief `complex` version of \ref S_Params_Volume_Structor.
struct S_Params_Volume_Structor_c {
	constructor_sol_vc_c_fptr constructor_sol_vc; ///< Pointer to the appropriate function.
};

/// \brief `complex` version of \ref Flux_Ref.
struct Flux_Ref_c {
	const struct const_Multiarray_c* fr; ///< See brief.
};

// Interface functions ********************************************************************************************** //

/** \brief `complex` version of \ref constructor_Flux_Ref_vol.
 *  \return See brief. */
struct Flux_Ref_c* constructor_Flux_Ref_vol_c
	(const struct S_Params_Volume_Structor_c* spvs, ///< \ref S_Params_Volume_Structor.
	 struct Flux_Input_c* flux_i,                   ///< \ref Flux_Input.
	 const struct Solver_Volume* s_vol,             ///< \ref Solver_Volume.
	 const struct Simulation* sim                   ///< \ref Simulation.
	);

/// \brief `complex` version of \ref destructor_Flux_Ref.
void destructor_Flux_Ref_c
	(struct Flux_Ref_c* flux_ref ///< See brief.
	);

#endif // DPG__test_complex_compute_volume_rhs_h__INCLUDED
