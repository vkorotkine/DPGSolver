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
 *  \brief Provides the macro definitions used for c-style templating related to the `double` \ref DG_Solver_Volume
 *         containers/functions.
 */

#if TYPE_RC == TYPE_REAL

///\{ \name Data types
#define DG_Solver_Volume_T DG_Solver_Volume
///\}

///\{ \name Function names
#define constructor_derived_DG_Solver_Volume_T constructor_derived_DG_Solver_Volume
#define destructor_derived_DG_Solver_Volume_T  destructor_derived_DG_Solver_Volume
///\}

#elif TYPE_RC == TYPE_COMPLEX
#endif
