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
 *  \brief Provides the macro definitions used for c-style templating related to the \ref Solver_Face
 *         containers/functions.
 */

#if TYPE_RC == TYPE_REAL

///\{ \name Data types
#define Solver_Face_T Solver_Face
///\}

///\{ \name Function names
#define constructor_derived_Solver_Face_T constructor_derived_Solver_Face
#define destructor_derived_Solver_Face_T  destructor_derived_Solver_Face
///\}

#elif TYPE_RC == TYPE_COMPLEX

///\{ \name Data types
#define Solver_Face_T Solver_Face_c
///\}

///\{ \name Function names
#define constructor_derived_Solver_Face_T constructor_derived_Solver_Face_c
#define destructor_derived_Solver_Face_T  destructor_derived_Solver_Face_c
///\}

#endif
