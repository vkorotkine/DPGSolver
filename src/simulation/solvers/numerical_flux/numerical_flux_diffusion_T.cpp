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
 *  \brief Provides the templated diffusion numerical flux functions.
 */

#include <assert.h>
#include <stddef.h>

#include "macros.h"
#include "definitions_core.h"


#include "def_templates_numerical_flux.h"

// Static function declarations ************************************************************************************* //

#define NEQ 1 ///< Number of equations.
#define NVR 1 ///< Number of variables.

// Interface functions ********************************************************************************************** //

#include "numerical_flux_central_T.cpp"

void compute_Numerical_Flux_T_diffusion_central
	(const struct Numerical_Flux_Input_T* num_flux_i, struct mutable_Numerical_Flux_T* num_flux)
{
	compute_Numerical_Flux_T_central(num_flux_i,num_flux);
}

void compute_Numerical_Flux_T_diffusion_central_jacobian
	(const struct Numerical_Flux_Input_T* num_flux_i, struct mutable_Numerical_Flux_T* num_flux)
{
	compute_Numerical_Flux_T_central_jacobian(num_flux_i,num_flux);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

#include "undef_templates_numerical_flux.h"
