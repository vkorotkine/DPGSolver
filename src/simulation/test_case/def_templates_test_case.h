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
 *  \brief Provides the macro definitions used for c-style templating related to the \ref Test_Case_T functions.
 */

#if TYPE_RC == TYPE_REAL

///\{ \name Data types
#define Test_Case_T Test_Case
///\}

///\{ \name Function names
#define set_string_associations_T              set_string_associations
#define set_pde_related_T                      set_pde_related
#define set_function_pointers_T                set_function_pointers
#define read_test_case_parameters_T            read_test_case_parameters
#define set_string_associations_test_case_T    set_string_associations_test_case
#define correct_invalid_test_case_parameters_T correct_invalid_test_case_parameters

#define constructor_Test_Case_T constructor_Test_Case
#define destructor_Test_Case_T  destructor_Test_Case
///\}

#elif TYPE_RC == TYPE_COMPLEX

#endif
