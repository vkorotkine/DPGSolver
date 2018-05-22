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

const struct const_Multiarray_d *grad_xyz_NURBS_patch_mapping(
	const struct const_Multiarray_d* xi_eta_i, int P, int Q, 
	const struct const_Multiarray_d* knots_xi, 
	const struct const_Multiarray_d* knots_eta, 
	const struct const_Multiarray_d* control_points_and_weights,
	const struct const_Multiarray_i* control_point_connectivity){

	/*
	TODO: Need to template this function in order to take into consideration complex 
	perturbations

	TODO: Put this function in geometry_parametric_T.c

	Compute the gradient of the mapping using the positions xi_eta_i on the
	parametric domain.

	Arguments:
		xi_eta_i = The values on the parametric domain (knot domain). Is a mutliarray
			of dimension [num_points x DIM]. For now, only consider 2D cases.
		P = The order of the basis functions in the xi direction
		Q = The order of the basis functions in the eta direction
		knots_xi = The knot vector in the xi direction
		knots_eta = The knot vector in the eta direction
		control_points_and_weights = The multiarray holding the list of control points
			and their associated weights. The multiarray is of dimension [num_pts x (DIM + 1)]
		control_point_connectivity = The multiarray holding the connectivity information for
			the control points and weights. Is of dimension [num_I x num_J], where num_I and 
			num_J are the number of control points / basis functions in either coordinate 
			direction. This multiarray is in row major form
	
	Return:
		A multiarray of dimension [num_points x (DIM^2)] that holds the partial derivatives
		of the mapping. In the 2D case (what has been implemented so far), the ith 
		row of the multiarray holds data of the form [x_xi, y_xi, x_eta, y_eta], where 
		x_xi, for instance, is the partial of x with respect to the xi knot domain variable.
	*/


	// Preprocessing:

	// Create multiarray structures to hold the control points and weights. The
	// control point structure will be of dimension [num_I x num_J] and there will be two
	// separate structures for the x and y values. The weight structure is of order 
	// [num_I x num_J].

	int i, j, k, numI, numJ, control_pt_index;
	numI = (int) control_point_connectivity->extents[0];
	numJ = (int) control_point_connectivity->extents[1];

	struct Multiarray_d *control_pts_x = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){numI, numJ});  // free
	struct Multiarray_d *control_pts_y = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){numI, numJ});  // free
	struct Multiarray_d *control_pts_w = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){numI, numJ});  // free

	const double *const control_pts_x_list = get_col_const_Multiarray_d(0, control_points_and_weights);
	const double *const control_pts_y_list = get_col_const_Multiarray_d(1, control_points_and_weights);
	const double *const control_pts_w_list = get_col_const_Multiarray_d(2, control_points_and_weights);
	
	for (i = 0; i < numI; i++){
		for (j = 0; j < numJ; j++){

			// The index of the control point in the list
			control_pt_index = get_row_const_Multiarray_i(i, control_point_connectivity)[j]; 

			get_col_Multiarray_d(j, control_pts_x)[i] = control_pts_x_list[control_pt_index];
			get_col_Multiarray_d(j, control_pts_y)[i] = control_pts_y_list[control_pt_index];
			get_col_Multiarray_d(j, control_pts_w)[i] = control_pts_w_list[control_pt_index];

		}
	}


	// Compute the partial derivative terms:
	int num_pts = (int) xi_eta_i->extents[0];
	struct Multiarray_d *grad_xyz = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){num_pts, 4});  // return

	// The values on the parametric domain to obtain partial derivatives at
	const double *const xi_vals_list = get_col_const_Multiarray_d(0, xi_eta_i);
	const double *const eta_vals_list = get_col_const_Multiarray_d(1, xi_eta_i);

	const struct const_Multiarray_d* grad_vals;
	double x_xi, y_xi, x_eta, y_eta;

	for (k = 0; k < num_pts; k++){
		// Loop over all the points

		const double xi = xi_vals_list[k];
		const double eta = eta_vals_list[k];

		x_xi = 0;
		y_xi = 0;
		x_eta = 0;
		y_eta = 0;

		for (i = 0; i < numI; i++){
			for (j = 0; j < numJ; j++){

				// Evaluate the i,j NURBS basis function gradient
				grad_vals = grad_NURBS_basis_ij_pq(i, j, P, Q, xi, eta, 
					knots_xi, knots_eta, (const struct const_Multiarray_d*)control_pts_w);

				// partial with respect to xi:
				x_xi += get_col_const_Multiarray_d(0, grad_vals)[0] * get_col_Multiarray_d(j, control_pts_x)[i];
				y_xi += get_col_const_Multiarray_d(0, grad_vals)[0] * get_col_Multiarray_d(j, control_pts_y)[i];

				// partial with respect to eta:
				x_eta += get_col_const_Multiarray_d(0, grad_vals)[1] * get_col_Multiarray_d(j, control_pts_x)[i];
				y_eta += get_col_const_Multiarray_d(0, grad_vals)[1] * get_col_Multiarray_d(j, control_pts_y)[i];

			}
		}

		// Store the partial terms
		get_col_Multiarray_d(0, grad_xyz)[k] = x_xi;
		get_col_Multiarray_d(1, grad_xyz)[k] = y_xi;
		get_col_Multiarray_d(2, grad_xyz)[k] = x_eta;
		get_col_Multiarray_d(3, grad_xyz)[k] = y_eta;

	}

	// Free allocated arrays
	destructor_Multiarray_d(control_pts_x);
	destructor_Multiarray_d(control_pts_y);
	destructor_Multiarray_d(control_pts_w);

	return (const struct const_Multiarray_d*)grad_xyz;

}

const struct const_Multiarray_d *xyz_NURBS_patch_mapping(
	const struct const_Multiarray_d* xi_eta_i, int P, int Q, 
	const struct const_Multiarray_d* knots_xi, 
	const struct const_Multiarray_d* knots_eta, 
	const struct const_Multiarray_d* control_points_and_weights,
	const struct const_Multiarray_i* control_point_connectivity){

	/*
	TODO: Need to template this function in order to take into consideration complex 
	perturbations

	TODO: Put this function in geometry_parametric_T.c

	Compute the values xyz on the physical domain using the positions xi_eta_i on the 
	parametric domain.

	Arguments:
		xi_eta_i = The values on the parametric domain (knot domain). Is a mutliarray
			of dimension [num_points x DIM]. For now, only consider 2D cases.
		P = The order of the basis functions in the xi direction
		Q = The order of the basis functions in the eta direction
		knots_xi = The knot vector in the xi direction
		knots_eta = The knot vector in the eta direction
		control_points_and_weights = The multiarray holding the list of control points
			and their associated weights. The multiarray is of dimension [num_pts x (DIM + 1)]
		control_point_connectivity = The multiarray holding the connectivity information for
			the control points and weights. Is of dimension [num_I x num_J], where num_I and 
			num_J are the number of control points / basis functions in either coordinate 
			direction. This multiarray is in row major form

	Return:
		A multiarray of dimension [num_points x DIM] that holds the location of
		the points on the physical domain. 
	*/

	// Preprocessing:

	// Create multiarray structures to hold the control points and weights. The
	// control point structure will be of dimension [num_I x num_J] and there will be two
	// separate structures for the x and y values. The weight structure is of order 
	// [num_I x num_J].

	int i, j, k, numI, numJ, control_pt_index;
	numI = (int) control_point_connectivity->extents[0];
	numJ = (int) control_point_connectivity->extents[1];

	struct Multiarray_d *control_pts_x = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){numI, numJ});  // free
	struct Multiarray_d *control_pts_y = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){numI, numJ});  // free
	struct Multiarray_d *control_pts_w = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){numI, numJ});  // free

	const double *const control_pts_x_list = get_col_const_Multiarray_d(0, control_points_and_weights);
	const double *const control_pts_y_list = get_col_const_Multiarray_d(1, control_points_and_weights);
	const double *const control_pts_w_list = get_col_const_Multiarray_d(2, control_points_and_weights);
	
	for (i = 0; i < numI; i++){
		for (j = 0; j < numJ; j++){

			// The index of the control point in the list
			control_pt_index = get_row_const_Multiarray_i(i, control_point_connectivity)[j]; 

			get_col_Multiarray_d(j, control_pts_x)[i] = control_pts_x_list[control_pt_index];
			get_col_Multiarray_d(j, control_pts_y)[i] = control_pts_y_list[control_pt_index];
			get_col_Multiarray_d(j, control_pts_w)[i] = control_pts_w_list[control_pt_index];

		}
	}

	// Perform the Mapping:

	int num_pts = (int) xi_eta_i->extents[0];

	struct Multiarray_d *xyz = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){num_pts, 2});  // return

	const double *const xi_vals_list = get_col_const_Multiarray_d(0, xi_eta_i);
	const double *const eta_vals_list = get_col_const_Multiarray_d(1, xi_eta_i);

	double x, y, basis_val;
	for (k = 0; k < num_pts; k++){
		// Loop over all the points

		const double xi = xi_vals_list[k];
		const double eta = eta_vals_list[k];

		x = 0;
		y = 0;

		for (i = 0; i < numI; i++){
			for (j = 0; j < numJ; j++){

				// Evaluate the i,j NURBS basis function
				basis_val = NURBS_basis_ij_pq(i, j, P, Q, xi, eta, 
					knots_xi, knots_eta, (const struct const_Multiarray_d*)control_pts_w);

				x += basis_val * get_col_Multiarray_d(j, control_pts_x)[i];
				y += basis_val * get_col_Multiarray_d(j, control_pts_y)[i];

			}
		}

		// Store the mapped values
		get_col_Multiarray_d(0, xyz)[k] = x;
		get_col_Multiarray_d(1, xyz)[k] = y;

	}

	// Free allocated arrays
	destructor_Multiarray_d(control_pts_x);
	destructor_Multiarray_d(control_pts_y);
	destructor_Multiarray_d(control_pts_w);

	return (struct const_Multiarray_d*) xyz;

}


void test_mapping(){

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

