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

#ifndef DPG__core_h__INCLUDED
#define DPG__core_h__INCLUDED
/** \file
 *  Provides core functions.
 */

/** \brief Return a statically allocated `char*` holding the name of the PETSc options file with full path.
 *  \return See brief. */
const char* set_petsc_options_name
	(const char* petsc_options_name ///< The name of the PETSc options file (excluding the path).
	);

#endif // DPG__core_h__INCLUDED
