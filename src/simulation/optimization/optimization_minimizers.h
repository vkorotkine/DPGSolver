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


#ifndef DPG__optimization_minimizers_h__INCLUDED
#define DPG__optimization_minimizers_h__INCLUDED

struct Optimization_Case;

void preprocessor_minimizer(struct Optimization_Case* optimization_case);

void postprocessor_minimizer(struct Optimization_Case* optimization_case);

void gradient_descent(struct Optimization_Case *optimization_case, int design_iteration);

void BFGS_minimizer(struct Optimization_Case *optimization_case, int design_iteration);

#endif // DPG__optimization_minimizers_h__INCLUDED
