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
 *  \brief Provides the macro definitions used for c-style templating related to the `double` \ref OPG_Solver_Volume_T
 *         containers/functions.
 */

#include "def_templates_volume_solver.h"

#if TYPE_RC == TYPE_REAL

///\{ \name Data types
#define OPG_Solver_Volume_T OPG_Solver_Volume
///\}

///\{ \name Function names
#define constructor_derived_OPG_Solver_Volume_T constructor_derived_OPG_Solver_Volume
#define destructor_derived_OPG_Solver_Volume_T  destructor_derived_OPG_Solver_Volume
#define get_operator__vc0_vs_vs_T get_operator__vc0_vs_vs_d
#define get_operator__cv1_vt_vs_T get_operator__cv1_vt_vs_d
///\}

#elif TYPE_RC == TYPE_COMPLEX

///\{ \name Data types
#define OPG_Solver_Volume_T OPG_Solver_Volume_c
///\}

///\{ \name Function names
#define constructor_derived_OPG_Solver_Volume_T constructor_derived_OPG_Solver_Volume_c
#define destructor_derived_OPG_Solver_Volume_T  destructor_derived_OPG_Solver_Volume_c
#define get_operator__vc0_vs_vs_T get_operator__vc0_vs_vs_c
#define get_operator__cv1_vt_vs_T get_operator__cv1_vt_vs_c
///\}

#endif
