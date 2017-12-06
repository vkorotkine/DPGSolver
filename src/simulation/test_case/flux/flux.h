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

#ifndef DPG__flux_h__INCLUDED
#define DPG__flux_h__INCLUDED
/** \file
 *  \brief Provides real containers and functions relating to fluxes of the supported PDEs.
 */

///\{ \name The maximum number of outputs from the flux functions.
#define MAX_FLUX_OUT 3 ///< See the members of \ref Flux.
///\}

#include "def_templates_type_d.h"
#include "def_templates_multiarray.h"
#include "def_templates_flux.h"
#include "flux_T.h"
#include "undef_templates_type.h"
#include "undef_templates_multiarray.h"
#include "undef_templates_flux.h"

#endif // DPG__flux_h__INCLUDED
