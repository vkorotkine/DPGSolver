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

#include <stdlib.h>
#include <stdio.h>

#include "test_base.h"
#include "macros.h"

#include "test_integration_advection.h"
#include "test_integration_convergence.h"

void test_integration_advection (struct Test_Info*const test_info)
{
	if (test_info->t_int.equivalence_real_complex) {
		test_print_warning(test_info,"Advection equivalence real/complex testing not yet implemented");
	} else {
		test_print_warning(test_info,"Advection equivalence real/complex testing currently disabled");
	}

	if (test_info->t_int.equivalence_algorithms) {
		test_print_warning(test_info,"Advection equivalence algorithms testing not yet implemented");
	} else {
		test_print_warning(test_info,"Advection equivalence algorithms testing currently disabled");
	}

	if (test_info->t_int.linearization) {
		test_print_warning(test_info,"Advection linearization testing not yet implemented");
	} else {
		test_print_warning(test_info,"Advection linearization testing currently disabled");
	}

	if (test_info->t_int.conv_order) {
		test_integration_convergence(test_info,"advection/TEST_Advection_Peterson_TRI__ml0__p1");
	} else {
		test_print_warning(test_info,"Advection convergence order testing currently disabled");
	}
}
