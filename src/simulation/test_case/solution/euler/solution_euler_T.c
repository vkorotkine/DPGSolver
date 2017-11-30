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
 *  \brief Provides the templated Euler solution functions.
 */

#include <assert.h>

#include "macros.h"
#include "definitions_core.h"
#include "definitions_test_case.h"

// Static function declarations ************************************************************************************* //

#define NEQ  NEQ_EULER  ///< Number of equations.
#define NVAR NVAR_EULER ///< Number of variables.

// Interface functions ********************************************************************************************** //

void convert_variables_T (struct Multiarray_T* vars, const char type_i, const char type_o)
{
	assert(type_i != type_o);
	assert(vars->layout == 'C');

	const ptrdiff_t ext_0 = vars->extents[0];

	assert(vars->extents[1] == NVAR);

	switch (type_i) {
	case 'p': {
		Type* rho = get_col_Multiarray_T(0,vars),
		    * p   = get_col_Multiarray_T(NVAR-1,vars),
		    * E   = p;

		Type* uvw[DMAX] = {            get_col_Multiarray_T(1,vars),
		                    (DIM > 1 ? get_col_Multiarray_T(2,vars) : NULL),
		                    (DIM > 2 ? get_col_Multiarray_T(3,vars) : NULL), };

		Type* rhouvw[DMAX] = { uvw[0], uvw[1], uvw[2], };
		switch (type_o) {
		case 'c':
			for (ptrdiff_t i = 0; i < ext_0; ++i) {
				Type V2 = 0.0;
				for (int d = 0; d < DIM; ++d) {
					V2 += uvw[d][i]*uvw[d][i];
					rhouvw[d][i] = rho[i]*uvw[d][i];
				}
				E[i] = p[i]/GM1 + 0.5*rho[i]*V2;
			}
			break;
		case 'p':
			return;
			break;
		default:
			EXIT_ERROR("Unsupported: %c\n",type_o);
			break;
		}
		break;
	} case 'c': {
		Type* rho = get_col_Multiarray_T(0,vars),
		    * p   = get_col_Multiarray_T(NVAR-1,vars),
		    * E   = p;

		Type* uvw[DMAX] = {            get_col_Multiarray_T(1,vars),
		                    (DIM > 1 ? get_col_Multiarray_T(2,vars) : NULL),
		                    (DIM > 2 ? get_col_Multiarray_T(3,vars) : NULL), };

		Type* rhouvw[DMAX] = { uvw[0], uvw[1], uvw[2], };
		switch (type_o) {
		case 'p':
			for (ptrdiff_t i = 0; i < ext_0; ++i) {
				const Type rho_inv = 1.0/rho[i];
				Type rho2V2 = 0.0;
				for (int d = 0; d < DIM; ++d) {
					rho2V2 += rhouvw[d][i]*rhouvw[d][i];
					uvw[d][i] = rhouvw[d][i]*rho_inv;
				}
				p[i] = GM1*(E[i]-0.5*rho2V2*rho_inv);
			}
			break;
		case 'c':
			return;
			break;
		default:
			EXIT_ERROR("Unsupported: %c\n",type_o);
			break;
		}
		break;
	} default:
		EXIT_ERROR("Unsupported: %c\n",type_i);
		break;
	}
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
