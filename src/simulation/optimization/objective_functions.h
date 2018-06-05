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

#ifndef DPG__objective_functions_h__INCLUDED
#define DPG__objective_functions_h__INCLUDED

#include <complex.h>

struct Simulation;
struct Intrusive_List;

/** \brief Pointer to objective functions (complex version) used for the optimization
 *
 *  \param sim    \ref Simulation.
 */

typedef double complex (*objective_function_fptr_c)
	(const struct Simulation* sim
	);


/** \brief Pointer to objective functions used for the optimization
 *
 *  \param sim    \ref Simulation.
 */

typedef double (*objective_function_fptr)
	(const struct Simulation* sim
	);


// Real and complex versions of the target CL objective function
double objective_function_target_Cl(const struct Simulation* sim); 
double complex objective_function_target_Cl_c(const struct Simulation* sim); 


#endif // DPG__objective_functions_h__INCLUDED


