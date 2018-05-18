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
 */

#include "geometry_NURBS_parametric.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "macros.h"

#include "multiarray.h"

#include "bases.h"

void test(){

	// Create a simple B spline test here to see that the basis function
	// is evaluating correctly

	int n_knots = 9;  // The length of the knot vector
	int P = 2;

	struct Multiarray_d* xi_vals = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){4,1} ); 
	double *xi_vals_i = get_col_Multiarray_d(0, xi_vals);

	xi_vals_i[0] = -2.0;
	xi_vals_i[1] = -0.5;
	xi_vals_i[2] = 0.25;
	xi_vals_i[3] = 2.0;

	struct Multiarray_d* weights = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){6,1} ); 
	double *weights_i = get_col_Multiarray_d(0, weights);

	weights_i[0] = 0.5;
	weights_i[1] = 1.75;
	weights_i[2] = 0.8;
	weights_i[3] = 1.8;
	weights_i[4] = 1.75;
	weights_i[5] = 1.5;

	struct Multiarray_d* knots = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){n_knots,1}); 
	double *knots_i = get_col_Multiarray_d(0, knots);
	
	knots_i[0] = -2.0;
	knots_i[1] = -2.0;
	knots_i[2] = -2.0;
	knots_i[3] = -1.0;
	knots_i[4] =  0.0;
	knots_i[5] =  1.0;
	knots_i[6] =  2.0;
	knots_i[7] =  2.0;
	knots_i[8] =  2.0;

	const struct const_Multiarray_d *B_Spline_basis_values;
	B_Spline_basis_values = B_Spline_Basis_p(P, (struct const_Multiarray_d*)xi_vals, 
		(struct const_Multiarray_d*)knots);

	printf("NURBS values: \n");
	print_const_Multiarray_d_tol(NURBS_Basis_p(B_Spline_basis_values, (struct const_Multiarray_d*)weights), 1E-4);

	//printf(" NUBRS Basis Val : %f \n", NURBS_Basis_ip(1, xi_val_basis_values, (struct const_Multiarray_d*)weights)); 

	//printf(" B Spline Basis Val : %f \n", B_Spline_Basis_ip(0, 2, xi_val, (struct const_Multiarray_d*)knots)); 

}

