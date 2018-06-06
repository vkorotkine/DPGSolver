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

#ifndef DPG__solve_implicit_h__INCLUDED
#define DPG__solve_implicit_h__INCLUDED
/** \file
 *  \brief Provides the interface to functions used to solve for the solution using implicit procedures.
 *
 *  An outline of exactly what terms are stored for the implicit solve follows. The code was originally implemented for
 *  explicit schemes where the unsteady term is to be updated according to:
 *  \f[
 *  	\frac{\partial s_{coef}}{\partial t} = \text{rhs}.
 *  \f]
 *
 *  For implicit schemes, the linearization was then performed such that:
 *  \f[
 *  	\text{lhs} = \frac{\partial \text{rhs}}{\partial s_{coef}}.
 *  \f]
 *
 *  As the full Newton method obtains the solution updates using the following equation:
 *  \f{eqnarray*}{
 *  	0 = \text{rhs}(s_{coef,exact}) &\approx& \text{rhs}(s_{coef})
 *  	                                + \frac{\partial \text{rhs}}{\partial s_{coef}}|_{s_{coef}}\ \Delta(s_{coef})\\
 *  	                               &\approx& \text{rhs}(s_{coef}) + \text{lhs}(s_{coef})\ \Delta(s_{coef}),
 *  \f}
 *
 *  the corresponding linear system is:
 *  \f[
 *  	\text{lhs}(s_{coef}) \Delta(s_{coef}) = -\text{rhs}(s_{coef}).
 *  \f]
 *
 *  Interpretting the above in the general form of a linear system Ax = b, we thus have:
 *  - A = \f$ \text{lhs}(s_{coef}) \f$;
 *  - x = \f$ \Delta(s_{coef}) \f$; and
 *  - b = \f$ -\text{rhs}(s_{coef}) \f$.
 */

#include <stddef.h>
#include <stdbool.h>
#include "petscmat.h"

struct Simulation;
struct Solver_Storage_Implicit;
struct Vector_i;

/// \brief Solve for the solution using an implicit solver.
void solve_implicit
	(struct Simulation* sim ///< \ref Simulation.
	);

/** \brief Check whether the matrix under consideration is symmetric based on \ref Simulation::method and
 *         \ref Test_Case_T::pde_index.
 *  \return `true` if symmetric; `false` otherwise. */
bool check_symmetric
	(const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Output a PETSc Mat to the file of given input name.
void output_petsc_mat
	(Mat A,                ///< The PETSc Mat.
	 const char* file_name ///< The file name.
	);

#endif // DPG__solve_implicit_h__INCLUDED
