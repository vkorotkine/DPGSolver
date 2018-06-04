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

#ifndef DPG__numerical_flux_h__INCLUDED
#define DPG__numerical_flux_h__INCLUDED
/** \file
 *  \brief Provides real containers and functions relating to the supported numerical fluxes.
 */

#include "boundary.h"

#include "def_templates_type_d.h"
#include "def_templates_numerical_flux.h"
#include "def_templates_multiarray.h"
#include "def_templates_boundary.h"
#include "def_templates_flux.h"
#include "numerical_flux_T.h"
#include "undef_templates_type.h"
#include "undef_templates_numerical_flux.h"
#include "undef_templates_multiarray.h"
#include "undef_templates_boundary.h"
#include "undef_templates_flux.h"

/** \brief Return a statically allocated array of `int*` holding the values of the numerical flux indices.
 *  \return See brief.
 *
 *  Passing a non-NULL input for `new_vals` sets the statically allocated array to the contained values.
 */
const int* get_set_ind_num_flux
	(const int*const new_vals ///< New values.
	);

#endif // DPG__numerical_flux_h__INCLUDED
