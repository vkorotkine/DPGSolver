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

#include "const_cast.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

/// \brief Cast a single variable of input type to a variable of `const` input type.
void const_cast_T
	(const Type* dest, ///< Destination.
	 const Type src    ///< Source.
	)
{
	*(Type*) dest = src;
}

/// \brief Cast `n` variables of input type to a variable of `const` input type.
void const_cast_T_n
	(const Type* dest, ///< Destination.
	 const Type* src,  ///< Source.
	 const int n       ///< The number of variables.
	)
{
	for (int i = 0; i < n; ++i)
		const_cast_T(&dest[i],src[i]);
}

/// \brief Cast a single pointer of input type to a `const` pointer of `const` input type.
void const_cast_T1
	(const Type*const* dest, ///< Destination.
	 const Type*const src    ///< Source.
	)
{
	*(Type**) dest = (Type*) src;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
