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

void test_integration_euler (struct Test_Info*const test_info)
{
	if (test_info->t_int.equivalence_real_complex) {
		;
	} else {
		test_print_warning(test_info,"Euler equivalence real/complex testing currently disabled");
	}

	if (test_info->t_int.equivalence_algorithms) {
		;
	} else {
		test_print_warning(test_info,"Euler equivalence algorithms testing currently disabled");
	}

	if (test_info->t_int.linearization) {
#if 1
 #if 0
		test_integration_linearization("test/euler/Test_Euler_SupersonicVortex_ToBeCurvedMIXED2D");
 #else
		test_integration_linearization("test/euler/Test_Euler_SupersonicVortex_CurvedMIXED2D");
 #endif
#else
		test_integration_linearization("test/euler/Test_Euler_PeriodicVortex_Stationary_QUAD");
#endif
	} else {
		test_print_warning(test_info,"Euler linearization testing currently disabled");
	}

	if (test_info->t_int.conv_order) {
		;
	} else {
		test_print_warning(test_info,"Euler convergence order testing currently disabled");
	}
}
