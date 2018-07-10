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


#include "optimize.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>
#include <time.h>

#include "face.h"
#include "volume.h"

#include "test_integration.h"

#include "definitions_adaptation.h"
#include "definitions_core.h"
#include "definitions_tol.h"

#include "macros.h"

#include "matrix.h"
#include "matrix_constructors.h"

#include "multiarray.h"
#include "multiarray_constructors.h"
#include "multiarray_math.h"

#include "intrusive.h"
#include "definitions_intrusive.h"

#include "geometry.h"
#include "geometry_parametric.h"
#include "solve_implicit.h"

#include "simulation.h"

#include "volume_solver.h"
#include "face_solver.h"
#include "visualization.h"
#include "definitions_visualization.h"


#include "optimization_case.h"
#include "optimizer_line_search.h"
#include "optimizer_NLPQLP.h"


// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void optimize(struct Simulation* sim){

	struct Optimization_Case *optimization_case = constructor_Optimization_Case(sim);

	// Optimization routine
	if (strstr(optimization_case->optimizer_spec, "LINE_SEARCH"))
		optimizer_line_search(optimization_case);
	else if(strstr(optimization_case->optimizer_spec, "NLPQLP"))
		optimizer_NLPQLP(optimization_case);
	else
		EXIT_UNSUPPORTED;

	// Clear allocated structures:
	destructor_Optimization_Case(optimization_case);

}



