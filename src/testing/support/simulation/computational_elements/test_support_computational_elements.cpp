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

#include "test_support_computational_elements.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "macros.h"

#include "intrusive.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

struct Intrusive_Link* constructor_copied_Intrusive_Link
	(struct Intrusive_Link* base, const size_t sizeof_base, const size_t sizeof_derived)
{
	struct Intrusive_Link* derived = calloc(1,sizeof_derived); // returned
	memcpy(derived,base,sizeof_base); // shallow copy of the base.

	assert(base->derived == NULL);

	return derived;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
