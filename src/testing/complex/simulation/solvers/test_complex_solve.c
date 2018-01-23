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

#include "test_complex_solve.h"

#include "complex_matrix.h"
#include "complex_multiarray.h"
#include "complex_vector.h"

#include "test_complex_face_solver.h"
#include "test_complex_volume_solver.h"

#include "test_complex_solve_dg.h"
#include "test_complex_solve_dpg.h"
#include "test_complex_test_case.h"


#include "solve.h"

#include "multiarray.h"
#include "vector.h"

#include "element_solver.h"

#include "intrusive.h"
#include "multiarray_operator.h"
#include "simulation.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

#include "def_templates_type_dc.h"
#include "solve_T.c"

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
