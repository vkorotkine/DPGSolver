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

#include "test_integration_euler.h"
#include "test_integration_linearization.h"
#include "test_integration_convergence.h"

void test_integration_euler (struct Test_Info*const test_info)
{
	if (test_info->t_int.equivalence_real_complex) {
		test_print_warning(test_info,"Euler equivalence real/complex testing not yet implemented");
	} else {
		test_print_warning(test_info,"Euler equivalence real/complex testing currently disabled");
	}

	if (test_info->t_int.equivalence_algorithms) {
		test_print_warning(test_info,"Euler equivalence algorithms testing not yet implemented");
	} else {
		test_print_warning(test_info,"Euler equivalence algorithms testing currently disabled");
	}

	if (test_info->t_int.linearization) {
		test_print_warning(test_info,"Euler linearization testing not yet implemented");
//		test_integration_linearization("test/euler/Test_Euler_SupersonicVortex_ToBeCurvedMIXED2D");
//		test_integration_linearization("test/euler/Test_Euler_SupersonicVortex_CurvedMIXED2D");
//		test_integration_linearization("test/euler/Test_Euler_PeriodicVortex_Stationary_QUAD");
	} else {
		test_print_warning(test_info,"Euler linearization testing currently disabled");
	}

	if (test_info->t_int.conv_order) {
// remove 'stationary' from the name here as it is not the case.
		test_integration_convergence(test_info,"euler/TEST_Euler_PeriodicVortex_Stationary_QUAD__ml0__p2");
	} else {
		test_print_warning(test_info,"Euler convergence order testing currently disabled");
	}
}
