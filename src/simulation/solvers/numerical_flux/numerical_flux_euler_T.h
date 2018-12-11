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
 *  \brief Provides templated functions relating to euler numerical fluxes.
 */

#include "def_templates_numerical_flux.h"

struct Numerical_Flux_Input_T;
struct mutable_Numerical_Flux_T;

/** \brief Version of \ref compute_Numerical_Flux_fptr_T computing the numerical fluxes using the Lax-Friedrichs method.
 *  The implementation was copied from that of [Hesthaven et al.'s Nodal DG code][hest_lf].
 *
 *  <!-- References: -->
 *  [hest_lf]: https://github.cppom/tcew/nodal-dg/blob/master/Codes1.1/CFD2D/EulerLF2D.m
 */
void compute_Numerical_Flux_T_euler_lax_friedrichs
	(const struct Numerical_Flux_Input_T* num_flux_i, ///< See brief.
	 struct mutable_Numerical_Flux_T* num_flux        ///< See brief.
	);

/** \brief Version of \ref compute_Numerical_Flux_fptr_T computing the numerical fluxes and Jacobians using the
 *         Lax-Friedrichs method.
 *  See comments for \ref compute_Numerical_Flux_T_euler_lax_friedrichs. */
void compute_Numerical_Flux_T_euler_lax_friedrichs_jacobian
	(const struct Numerical_Flux_Input_T* num_flux_i, ///< See brief.
	 struct mutable_Numerical_Flux_T* num_flux        ///< See brief.
		);

/** \brief Version of \ref compute_Numerical_Flux_fptr_T computing the numerical fluxes using the Roe-Pike method.
 *  The implementation is based off of that explained in (Ch. 11.3, \cite Toro2009). */
void compute_Numerical_Flux_T_euler_roe_pike
	(const struct Numerical_Flux_Input_T* num_flux_i, ///< See brief.
	 struct mutable_Numerical_Flux_T* num_flux        ///< See brief.
	);

/** \brief Version of \ref compute_Numerical_Flux_fptr_T computing the numerical fluxes and Jacobians using the Roe-Pike
 *         method.
 *  See comments for \ref compute_Numerical_Flux_T_euler_roe_pike. */
void compute_Numerical_Flux_T_euler_roe_pike_jacobian
	(const struct Numerical_Flux_Input_T* num_flux_i, ///< See brief.
	 struct mutable_Numerical_Flux_T* num_flux        ///< See brief.
	);

#include "undef_templates_numerical_flux.h"
