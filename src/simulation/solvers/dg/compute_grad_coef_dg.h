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

#ifndef DPG__compute_grad_coef_dg_h__INCLUDED
#define DPG__compute_grad_coef_dg_h__INCLUDED
/** \file
 *  \brief Provides functions used for computing the weak gradient (auxiliary variable) terms required for the
 *         computation of the DG viscous fluxes.
 *
 *  For all 2nd order equations, the auxiliary variable,
 *
 *  \f[
 *  \text{grad_coef} = \text{grad_coef_v} + \sum_{f} \text{grad_coef_f},
 *  \f]
 *
 *  is first
 *  computed here and then substituted into the primary equation such that a primal formulation of the DG scheme is
 *  obtained. Thus, despite the presence of the two equations, a mixed formulation is not being used to solve the
 *  system; in that case, both `sol_coef` and `grad_coef` would be global unknowns. Based on the mathematical analysis
 *  in Arnold et al. \cite Arnold2002, restructuring the code to solve the system directly in the primal form would
 *  potentially have advantages in terms of:
 *  	1. Consistency with the mathematical theory;
 *  	2. Possible elimination of redundant storage/operations:
 *  		- storage of grad_coef, grad_coef_v, grad_coef_f;
 *  		- computation of the same symmetric operators twice;
 *  		- computation of both contributions of the lifting operator (for CDG2);
 *  		- reordering operation in the penalty terms.
 */

struct Simulation;

/** \brief Compute the weak gradient terms for the DG scheme for PDEs which have 2nd order terms.
 *
 *  Computes: \todo Update the list when working for both explicit/implicit.
 *  - \ref Solver_Volume_T::grad_coef;
 *  - \ref DG_Solver_Volume_T::grad_coef_v;
 *  - \ref DG_Solver_Face_T::grad_coef_f.
 */
void compute_grad_coef_dg
	(const struct Simulation* sim ///< \ref Simulation.
	);

#endif // DPG__compute_grad_coef_dg_h__INCLUDED
