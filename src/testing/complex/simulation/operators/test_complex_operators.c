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
 */

#include "test_complex_operators.h"

#include <assert.h>

#include "macros.h"

#include "complex_matrix.h"
#include "complex_multiarray.h"
#include "matrix.h"
#include "multiarray.h"

#include "multiarray_operator.h"
#include "operator.h"

// Templated functions ********************************************************************************************** //

#include "def_templates_type_dc.h"
#include "def_templates_multiarray.h"
#include "def_templates_operators_c.h"
#include "operator_T.c"
#include "undef_templates_type.h"
#include "undef_templates_multiarray.h"
#include "undef_templates_operators.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
