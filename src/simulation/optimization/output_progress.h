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


#ifndef DPG__output_progress_h__INCLUDED
#define DPG__output_progress_h__INCLUDED

#include <stdio.h>
#include "stdbool.h"

#include "optimization_case.h"
#include "gradient.h"


/** \brief Setup the optimization progress file. This file will store the 
 * convergence information for the optimization case. The file will be placed in 
 * the Paraview output directory for the given case for now and will be titled
 * optimization_convergence.txt.
 *
 * The header of the file data will also be set within this function.
 *
 * \return File pointer to the progress file
 */
FILE* constructor_optimization_progress_file(
	struct Optimization_Case *optimization_case ///< Consult optimization_case.h (holds the case information)
	);


/** \brief Close the optimization progress file.
 */
void destructor_optimization_progress_file(FILE *fp);


/** \brief Add the information for a given design iteration into the progress file. If 
 * output_to_stdout is set to true, then the information will also be printed to standard output
 * to monitor the progress.
 */
void progress_file_add_information(
	FILE *fp, ///< The 
	struct Optimization_Case *optimization_case, ///< Consult optimization_case.h (holds the case information)
	double L2_grad, ///< The L2 norm of the gradient
	double objective_func_value, ///< The objective function value
	double time_elapsed, ///< Total CPU time elapsed when this design iteration finished
	int design_iteration, ///< The current design iteration
	bool output_to_stdout ///< Boolean for whether the progress should also be printed to the screen
	);


/** \brief Output the NURBS patch information into a text file (Optimized_NURBS_Patch.txt)
 * 	which will be located in the paraview output directory for the case. 
 */
void output_NURBS_patch_information(
	struct Optimization_Case* optimization_case ///< Consult optimization_case.h (holds the case information)
	);


/** \brief Output the gradient of the objective function with resepct to the design 
 *	variables into a file titled Gradient.txt.
 */
void output_gradient(
	struct Optimization_Case* optimization_case, ///< Consult optimization_case.h (holds the case information)
	double *gradient ///< Array of doubles with the gradient of the objective function. Length is found using optimization_case
	);

#endif // DPG__output_progress_h__INCLUDED
