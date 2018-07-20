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


// Start with this file for the optimization objective function to test
// functionals out. This must be eventually templated

#include "functionals.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include "macros.h"
#include "simulation.h"

#include "volume.h"
#include "volume_solver.h"
#include "face_solver.h"
#include "element_solver.h"
#include "element_solution.h"

#include "multiarray.h"
#include "vector.h"

#include "intrusive.h"
#include "file_processing.h"
#include "operator.h"
#include "compute_face_rlhs.h"
#include "solution_euler.h"
#include "boundary.h"
#include "definitions_nodes.h"
#include "nodes.h"

#include "definitions_physics.h"
#include "definitions_math.h"


// Static function declarations ************************************************************************************* //

/// \brief Container for solution data relating to 'c'oefficients of 'd'rag, 'l'ift and 'm'oment
struct Sol_Data__c_dlm {
	double rho,      ///< The free stream density.
	       p,        ///< The free stream pressure.
	       mach,     ///< The free stream 'mach' number.
	       theta,    ///< The free stream flow angle in the xy-plane (in radians).
	       area_ref, ///< The reference area used to scale the coefficients.
	       cm_le_x,  ///< The leading edge x location (for moment computations about the LE)
	       cm_le_y;  ///< The leading edge y location (for moment computations about the LE)
};


/** \brief Container for the functional data required for the various functionals defined
 */
struct Functional_Data {
	
	double target_cl; ///< Target lift coefficient value for the target cl functional

	char target_pressure_distribution_file_abs_path[STRLEN_MAX]; ///< The absolute path to the file with the pressure distribution data
	struct Multiarray_d *target_pressure_distribution; ///< The xi and pressure values for the target
};


/** \brief Read all the functional data information from the optimization.data file.
 *
 * \return The functional data structure with the data from the optimization.data file
 */
static struct Functional_Data get_functional_data(
	);


/** \brief Read the file provided containing the target pressure distribution. Load the values
 * 	into the target_pressure_distribution Multiarray in functional_data. 
 *
 *	NOTE: If no file has been given, the function will exit.
 */
static void read_target_pressure_distribution_file(
	struct Functional_Data *functional_data ///< Standard. \ref Functional_Data
	);


/** \brief Compute the lift, drag and moment coefficient of the airfoil (2D)
 *
 * \return Multiarray with the values of the lift coefficients. Returned
 * 	multiarray is in column major form and stores the data in the
 *	form [cl, cd, cm].
 */
static const struct const_Multiarray_d* compute_cl_cd_cm(
	const struct Simulation* sim ///< Standard. \ref Simulation
	);


/** \brief Complex version of compute_cl_cd. Will remove this function when the templating is
 * 	complete.
 *
 * \return Refer to compute_cl_cd
 */
static const struct const_Multiarray_c* compute_cl_cd_cm_c(
	const struct Simulation* sim ///< Standard. \ref Simulation. NOTE: Simulation must hold complex structures.
	);


/** \brief Return the statically allocated \ref Sol_Data__c_dlm container.
 *  \return See brief. */
static struct Sol_Data__c_dlm get_sol_data__c_dlm( 
	);


/** \brief Read the relevant physical data for the Cd, Cl and CM computation
 */
static void read_data_c_dlm (
	struct Sol_Data__c_dlm*const sol_data ///< Container to hold solution data
	);


/** \brief Compute the volume of a given solver volume (svol). 
 *
 * \return The value of the volume of the given solver volume.
 */
static double compute_volume_of_s_vol(
	const struct Solver_Volume* s_vol ///< Standard. \ref Solver_Volume_T.
	);


/** \brief Complex version of compute_volume_of_s_vol. Uses the complex
 *	metric terms to compute the complex volume of the solver volume.
 *
 *	\return The complex value of the volume.
 */
static double complex compute_volume_of_s_vol_c(
	const struct Solver_Volume_c* s_vol_c ///< Standard. \ref Solver_Volume_T.
	);


// Interface functions ********************************************************************************************** //


double functional_cl(const struct Simulation* sim){

	const struct const_Multiarray_d* cl_cd_cm = compute_cl_cd_cm(sim);  // free
	double cl = cl_cd_cm->data[0];
	destructor_const_Multiarray_d(cl_cd_cm);

	return cl;
}


double complex functional_cl_c(const struct Simulation* sim_c){

	const struct const_Multiarray_c* cl_cd_cm = compute_cl_cd_cm_c(sim_c);  // free
	double complex cl = cl_cd_cm->data[0];
	destructor_const_Multiarray_c(cl_cd_cm);

	return cl;
}


double functional_cm_le(const struct Simulation* sim){

	const struct const_Multiarray_d* cl_cd_cm = compute_cl_cd_cm(sim);  // free
	double cm = cl_cd_cm->data[2];
	destructor_const_Multiarray_d(cl_cd_cm);

	return cm;
}


double complex functional_cm_le_c(const struct Simulation* sim_c){

	const struct const_Multiarray_c* cl_cd_cm = compute_cl_cd_cm_c(sim_c);  // free
	double complex cm = cl_cd_cm->data[2];
	destructor_const_Multiarray_c(cl_cd_cm);

	return cm;
}


double functional_mesh_volume(const struct Simulation* sim){

	double mesh_volume = 0;

	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next){
		const struct Solver_Volume*const s_vol = (struct Solver_Volume*) curr;
		mesh_volume += compute_volume_of_s_vol(s_vol);
	}

	return mesh_volume;
}


double complex functional_mesh_volume_c(const struct Simulation* sim_c){

	double complex mesh_volume = 0;

	for (struct Intrusive_Link* curr = sim_c->volumes->first; curr; curr = curr->next){
		const struct Solver_Volume_c*const s_vol_c = (struct Solver_Volume_c*) curr;
		mesh_volume += compute_volume_of_s_vol_c(s_vol_c);
	}

	return mesh_volume;
}


double functional_target_cl(const struct Simulation* sim){

	struct Functional_Data functional_data = get_functional_data();
	double target_cl = functional_data.target_cl;
	double cl = functional_cl(sim);

	return 0.5 * (cl - target_cl) * (cl - target_cl);
}


double complex functional_target_cl_c(const struct Simulation* sim_c){

	struct Functional_Data functional_data = get_functional_data();
	double target_cl = functional_data.target_cl;
	double complex cl = functional_cl_c(sim_c);

	return 0.5 * (cl - target_cl) * (cl - target_cl);
}


double functional_mesh_volume_fractional_change(const struct Simulation* sim){

	double V = functional_mesh_volume(sim);
	double V0 = sim->mesh_volume_initial;

	return V/V0;
}


double complex functional_mesh_volume_fractional_change_c(const struct Simulation* sim_c){
	
	double complex V = functional_mesh_volume_c(sim_c);
	double V0 = sim_c->mesh_volume_initial;

	return V/V0;
}


double functional_inverse_pressure_design(const struct Simulation *sim){

	struct Functional_Data functional_data = get_functional_data();

	double inverse_pressure_design_functional_value = 0.0;

	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {

		// Loop through the faces that are on the wall boundary
		const struct Face*const face = (struct Face*) curr;
		const struct Solver_Face*const s_face = (struct Solver_Face*) curr;
		if (!is_face_wall_boundary(face))
			continue;

		const struct Volume*const vol = (struct Volume*) face->neigh_info[0].volume;

		const struct const_Element*const el_vol      = vol->element;
		const struct const_Multiarray_d*const xyz_ve = vol->xyz_ve;

		const struct const_Vector_i*const f_ve_lf = el_vol->f_ve->data[face->neigh_info[0].ind_lf];

		const double*const xyz_ve_a[] = { get_row_const_Multiarray_d(f_ve_lf->data[0],xyz_ve),
		                                  get_row_const_Multiarray_d(f_ve_lf->data[1],xyz_ve), };


		// Get the xi limits for this face
		double 	xi_val1  = xyz_ve_a[0][0],
				xi_val2  = xyz_ve_a[1][0],
				eta_val1 = xyz_ve_a[0][1],
				eta_val2 = xyz_ve_a[1][1];

		double xi_range[2], eta_range[2];
		xi_range[0]  = xi_val1  > xi_val2  ? xi_val2  : xi_val1;
		xi_range[1]  = xi_val1  > xi_val2  ? xi_val1  : xi_val2;
		eta_range[0] = eta_val1 > eta_val2 ? eta_val2 : eta_val1;
		eta_range[1] = eta_val1 > eta_val2 ? eta_val1 : eta_val2;


		// Obtain the face cubature (fc) nodes on the parametric knot domain
		const struct const_Nodes *nodes_c = constructor_const_Nodes_tp(face->element->d, s_face->p_ref, NODES_GL);
		const struct const_Multiarray_d *rst_i = constructor_move_const_Multiarray_d_Matrix_d(nodes_c->rst);  // free
		const ptrdiff_t n_fc = rst_i->extents[0];


		// \todo MSB: Map these values using some operator and be able to generalize to any type of
		//	element (not just QUADS)

		// Map the node point locations to knot domain:
		struct Multiarray_d *rst_knots_i = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){n_fc, 2}); // destructed
	
		// Flag is set to true if this is an eta = -1 face
		bool eta_min_1_face = false; 

		switch (face->neigh_info[0].ind_lf){
			case 0:
				// Face 0 : xi is fixed (left face)
				for (ptrdiff_t i = 0; i < n_fc; i++){
					get_col_Multiarray_d(0, rst_knots_i)[i] = xi_range[0];
					get_col_Multiarray_d(1, rst_knots_i)[i] = get_col_const_Multiarray_d(0, rst_i)[i]*(eta_range[1] - eta_range[0])*0.5 
																+ 0.5*(eta_range[1] + eta_range[0]);
				}
				break;

			case 1:
				// Face 1 : xi is fixed (right face)
				for (ptrdiff_t i = 0; i < n_fc; i++){
					get_col_Multiarray_d(0, rst_knots_i)[i] = xi_range[1];
					get_col_Multiarray_d(1, rst_knots_i)[i] = get_col_const_Multiarray_d(0, rst_i)[i]*(eta_range[1] - eta_range[0])*0.5 
																+ 0.5*(eta_range[1] + eta_range[0]);
				}
				break;

			case 2:
				// Face 2 : eta is fixed (bottom face)
				for (ptrdiff_t i = 0; i < n_fc; i++){
					get_col_Multiarray_d(0, rst_knots_i)[i] = get_col_const_Multiarray_d(0, rst_i)[i]*(xi_range[1] - xi_range[0])*0.5 
																+ 0.5*(xi_range[1] + xi_range[0]);
					get_col_Multiarray_d(1, rst_knots_i)[i] = eta_range[0];
				}
				eta_min_1_face = true;
				break;

			case 3:
				// Face 3 : eta is fixed (top face)
				for (ptrdiff_t i = 0; i < n_fc; i++){
					get_col_Multiarray_d(0, rst_knots_i)[i] = get_col_const_Multiarray_d(0, rst_i)[i]*(xi_range[1] - xi_range[0])*0.5 
																+ 0.5*(xi_range[1] + xi_range[0]);
					get_col_Multiarray_d(1, rst_knots_i)[i] = eta_range[1];
				}
				break;

			default:
				EXIT_ERROR("Unsupported face index"); break;
		}

		// Get the solution (s) at the face cubature (fc) nodes and convert it to primitive variables (pv)
		const struct const_Multiarray_d* pv_fc = constructor_s_fc_interp_d(0, s_face);
		convert_variables((struct Multiarray_d*)pv_fc,'c','p');

		// jacobian determinant values at the face cubature (fc) nodes
		const struct const_Multiarray_d* jac_det_fc = s_face->jacobian_det_fc;

		// Pressure (p) at face cubature (fc) nodes as an array
		const double *p_fc_i = get_col_const_Multiarray_d(pv_fc->extents[1]-1,pv_fc);

		// Cubature weights at the face cubature (fc) nodes
		const struct const_Vector_d*const w_fc = get_operator__w_fc__s_e(s_face);

		// The target values. s_values correspond to the location on the knot domain (xi coordinate)
		// of the given target value and p_values stores the target pressure values.
		double 	*target_s_values = get_col_Multiarray_d(0, functional_data.target_pressure_distribution),
				*target_p_values = get_col_Multiarray_d(1, functional_data.target_pressure_distribution);
		int num_target_values = (int) functional_data.target_pressure_distribution->extents[0];

		double p_value, s_value, p_target_value, p1, p2, s1, s2;

		if (eta_min_1_face){
		for (int i = 0; i < n_fc; i++){

			p_value = p_fc_i[i];
			s_value = get_col_Multiarray_d(0, rst_knots_i)[i];

			// Compute the target pressure at the ith cubature point using either interpolation or extrapolation
			
			if (target_s_values[0] > s_value){
				// Extrapolation Case 1
				
				s1 = target_s_values[0];
				p1 = target_p_values[0];

				s2 = target_s_values[1];
				p2 = target_p_values[1];

			} else if(target_s_values[num_target_values-1] < s_value){
				// Extrapolation Case 2

				s1 = target_s_values[num_target_values-2];
				p1 = target_p_values[num_target_values-2];

				s2 = target_s_values[num_target_values-1];
				p2 = target_p_values[num_target_values-1];

			} else{
				// Interpolation Case

				for (int j = 0; j < num_target_values-1; j++){

					if (target_s_values[j] < s_value && target_s_values[j+1] > s_value){
						s1 = target_s_values[j];
						p1 = target_p_values[j];

						s2 = target_s_values[j+1];
						p2 = target_p_values[j+1];

						break;
					}
				}
			}

			p_target_value = (s_value - s1) * ((p2 - p1)/ (s2 - s1)) + p1;

			inverse_pressure_design_functional_value += (p_value - p_target_value) * 
				(p_value - p_target_value) * jac_det_fc->data[i] * w_fc->data[i];
		}
		}

		// Free allocated memory
		destructor_const_Multiarray_d(rst_i);
		destructor_Multiarray_d(rst_knots_i);
		destructor_const_Nodes(nodes_c);
		destructor_const_Multiarray_d(pv_fc);
	}

	return inverse_pressure_design_functional_value;
}


double complex functional_inverse_pressure_design_c(const struct Simulation *sim_c){

	struct Functional_Data functional_data = get_functional_data();

	double complex inverse_pressure_design_functional_value = 0.0;

	for (struct Intrusive_Link* curr = sim_c->faces->first; curr; curr = curr->next) {

		// Loop through the faces that are on the wall boundary
		const struct Face*const face = (struct Face*) curr;
		const struct Solver_Face_c*const s_face_c = (struct Solver_Face_c*) curr;
		if (!is_face_wall_boundary(face))
			continue;

		// Face vertex points on the knot domain are always real
		const struct Volume*const vol = (struct Volume*) face->neigh_info[0].volume;
		const struct const_Element*const el_vol      = vol->element;
		const struct const_Multiarray_d*const xyz_ve = vol->xyz_ve;
		const struct const_Vector_i*const f_ve_lf = el_vol->f_ve->data[face->neigh_info[0].ind_lf];
		const double*const xyz_ve_a[] = { get_row_const_Multiarray_d(f_ve_lf->data[0],xyz_ve),
		                                  get_row_const_Multiarray_d(f_ve_lf->data[1],xyz_ve), };

		// Get the xi limits for this face
		double 	xi_val1  = xyz_ve_a[0][0],
				xi_val2  = xyz_ve_a[1][0],
				eta_val1 = xyz_ve_a[0][1],
				eta_val2 = xyz_ve_a[1][1];

		double xi_range[2], eta_range[2];
		xi_range[0]  = xi_val1  > xi_val2  ? xi_val2  : xi_val1;
		xi_range[1]  = xi_val1  > xi_val2  ? xi_val1  : xi_val2;
		eta_range[0] = eta_val1 > eta_val2 ? eta_val2 : eta_val1;
		eta_range[1] = eta_val1 > eta_val2 ? eta_val1 : eta_val2;


		// Obtain the face cubature (fc) nodes on the parametric knot domain
		const struct const_Nodes *nodes_c = constructor_const_Nodes_tp(face->element->d, s_face_c->p_ref, NODES_GL);
		const struct const_Multiarray_d *rst_i = constructor_move_const_Multiarray_d_Matrix_d(nodes_c->rst);  // free
		const ptrdiff_t n_fc = rst_i->extents[0];


		// \todo MSB: Map these values using some operator and be able to generalize to any type of
		//	element (not just QUADS)

		// Map the node point locations to knot domain:
		struct Multiarray_d *rst_knots_i = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){n_fc, 2}); // destructed
	
		// Flag is set to true if this is an eta = -1 face
		bool eta_min_1_face = false; 

		switch (face->neigh_info[0].ind_lf){
			case 0:
				// Face 0 : xi is fixed (left face)
				for (ptrdiff_t i = 0; i < n_fc; i++){
					get_col_Multiarray_d(0, rst_knots_i)[i] = xi_range[0];
					get_col_Multiarray_d(1, rst_knots_i)[i] = get_col_const_Multiarray_d(0, rst_i)[i]*(eta_range[1] - eta_range[0])*0.5 
																+ 0.5*(eta_range[1] + eta_range[0]);
				}
				break;

			case 1:
				// Face 1 : xi is fixed (right face)
				for (ptrdiff_t i = 0; i < n_fc; i++){
					get_col_Multiarray_d(0, rst_knots_i)[i] = xi_range[1];
					get_col_Multiarray_d(1, rst_knots_i)[i] = get_col_const_Multiarray_d(0, rst_i)[i]*(eta_range[1] - eta_range[0])*0.5 
																+ 0.5*(eta_range[1] + eta_range[0]);
				}
				break;

			case 2:
				// Face 2 : eta is fixed (bottom face)
				for (ptrdiff_t i = 0; i < n_fc; i++){
					get_col_Multiarray_d(0, rst_knots_i)[i] = get_col_const_Multiarray_d(0, rst_i)[i]*(xi_range[1] - xi_range[0])*0.5 
																+ 0.5*(xi_range[1] + xi_range[0]);
					get_col_Multiarray_d(1, rst_knots_i)[i] = eta_range[0];
				}
				eta_min_1_face = true;
				break;

			case 3:
				// Face 3 : eta is fixed (top face)
				for (ptrdiff_t i = 0; i < n_fc; i++){
					get_col_Multiarray_d(0, rst_knots_i)[i] = get_col_const_Multiarray_d(0, rst_i)[i]*(xi_range[1] - xi_range[0])*0.5 
																+ 0.5*(xi_range[1] + xi_range[0]);
					get_col_Multiarray_d(1, rst_knots_i)[i] = eta_range[1];
				}
				break;

			default:
				EXIT_ERROR("Unsupported face index"); break;
		}

		// Get the solution (s) at the face cubature (fc) nodes and convert it to primitive variables (pv)
		const struct const_Multiarray_c* pv_fc_c = constructor_s_fc_interp_c(0, s_face_c);
		convert_variables_c((struct Multiarray_c*)pv_fc_c,'c','p');

		// jacobian determinant values at the face cubature (fc) nodes
		const struct const_Multiarray_c* jac_det_fc_c = s_face_c->jacobian_det_fc;

		// Pressure (p) at face cubature (fc) nodes as an array
		const double complex *p_fc_c_i = get_col_const_Multiarray_c(pv_fc_c->extents[1]-1,pv_fc_c);

		// Cubature weights at the face cubature (fc) nodes
		const struct const_Vector_d*const w_fc = get_operator__w_fc__s_e_c(s_face_c);

		// The target values. s_values correspond to the location on the knot domain (xi coordinate)
		// of the given target value and p_values stores the target pressure values.
		double 	*target_s_values = get_col_Multiarray_d(0, functional_data.target_pressure_distribution),
				*target_p_values = get_col_Multiarray_d(1, functional_data.target_pressure_distribution);
		int num_target_values = (int) functional_data.target_pressure_distribution->extents[0];

		double complex p_value;
		double s_value, p_target_value, p1, p2, s1, s2;

		if(eta_min_1_face){
		for (int i = 0; i < n_fc; i++){

			p_value = p_fc_c_i[i];
			s_value = get_col_Multiarray_d(0, rst_knots_i)[i];

			// Compute the target pressure at the ith cubature point using either interpolation or extrapolation
			
			if (target_s_values[0] > s_value){
				// Extrapolation Case 1
				
				s1 = target_s_values[0];
				p1 = target_p_values[0];

				s2 = target_s_values[1];
				p2 = target_p_values[1];

			} else if(target_s_values[num_target_values-1] < s_value){
				// Extrapolation Case 2

				s1 = target_s_values[num_target_values-2];
				p1 = target_p_values[num_target_values-2];

				s2 = target_s_values[num_target_values-1];
				p2 = target_p_values[num_target_values-1];

			} else{
				// Interpolation Case

				for (int j = 0; j < num_target_values-1; j++){

					if (target_s_values[j] < s_value && target_s_values[j+1] > s_value){
						s1 = target_s_values[j];
						p1 = target_p_values[j];

						s2 = target_s_values[j+1];
						p2 = target_p_values[j+1];

						break;
					}
				}
			}

			p_target_value = (s_value - s1) * ((p2 - p1)/ (s2 - s1)) + p1;

			inverse_pressure_design_functional_value += (p_value - p_target_value) * 
				(p_value - p_target_value) * jac_det_fc_c->data[i] * w_fc->data[i];
		}
		}

		// Free allocated memory
		destructor_const_Multiarray_d(rst_i);
		destructor_Multiarray_d(rst_knots_i);
		destructor_const_Nodes(nodes_c);
		destructor_const_Multiarray_c(pv_fc_c);
	}

	return inverse_pressure_design_functional_value;
}


// Static functions ************************************************************************************************* //


static const struct const_Multiarray_d* compute_cl_cd_cm(const struct Simulation* sim){

	int i, num_fc;

	// Get the freestream flow data and reference area
	const struct Sol_Data__c_dlm sol_data = get_sol_data__c_dlm();

	const double rho_fs   = sol_data.rho,
	             p_fs     = sol_data.p,
	             c_fs     = sqrt(GAMMA*p_fs/rho_fs),
	             V_fs     = sol_data.mach*c_fs,
	             denom    = 0.5*rho_fs*V_fs*V_fs*sol_data.area_ref,
	             theta_fs = sol_data.theta,
	             cm_le_x  = sol_data.cm_le_x,
	             cm_le_y  = sol_data.cm_le_y;

	// The force components along the Normal and Axial directions and the moment about the leading edge
	double	Force_N = 0.0, 
			Force_A = 0.0,
			Moment  = 0.0;

	double delta_x, delta_y, x_fc, y_fc;

	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
		// Loop through the faces that are on the wall boundary
		const struct Face*const face = (struct Face*) curr;
		const struct Solver_Face*const s_face = (struct Solver_Face*) curr;
		if (!is_face_wall_boundary(face))
			continue;

		// Get the solution (s) at the face cubature (fc) nodes and convert it to 
		// primitive variables (pv)
		const struct const_Multiarray_d* pv_fc = constructor_s_fc_interp_d(0, s_face);
		convert_variables((struct Multiarray_d*)pv_fc,'c','p');
		
		// Get the normals (n) at the face cubature (fc) nodes and jacobian values
		const struct const_Multiarray_d* n_fc 		= s_face->normals_fc;
		const struct const_Multiarray_d* jac_det_fc = s_face->jacobian_det_fc;

		// Get the xyz values of the face cubature (fc) nodes
		const struct const_Multiarray_d* xyz_fc 	= s_face->xyz_fc;

		num_fc = (int) pv_fc->extents[0]; // number of face cubature (fc) nodes

		// Pressure (p) at face cubature (fc) nodes as an array
		const double *p_fc_i = get_col_const_Multiarray_d(pv_fc->extents[1]-1,pv_fc);

		// Integrate the pressure multiplied by the the normal component (into airfoil face).
		// Compute the integral using quadrature
		const struct const_Vector_d*const w_fc = get_operator__w_fc__s_e(s_face);

		for (i = 0; i < num_fc; i++){
			const double*const n = get_row_const_Multiarray_d(i,n_fc);

			// Force Computations
			Force_A += 1.0*p_fc_i[i]*n[0] * jac_det_fc->data[i] * w_fc->data[i];
			Force_N += 1.0*p_fc_i[i]*n[1] * jac_det_fc->data[i] * w_fc->data[i];


			// Moment Computations
			x_fc = get_col_const_Multiarray_d(0, xyz_fc)[i];
			y_fc = get_col_const_Multiarray_d(1, xyz_fc)[i];

			delta_x = x_fc - cm_le_x;
			delta_y = y_fc - cm_le_y;

			Moment += (delta_x*p_fc_i[i]*n[1] - delta_y*p_fc_i[i]*n[0]) * jac_det_fc->data[i] * w_fc->data[i];
		}

		// Destroy the allocated vectors
		destructor_const_Multiarray_d(pv_fc);
	}

	// Get the lift and drag by taking into account the angle of attack 
	double 	L = Force_N*cos(theta_fs) - Force_A*sin(theta_fs),
			D = Force_N*sin(theta_fs) + Force_A*cos(theta_fs);

	double 	cl = L/denom,
			cd = D/denom,
			cm = -1.0*Moment/(denom*sol_data.area_ref);

 	struct Multiarray_d *aero_coeffs = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){3, 1}); // destructed
	get_col_Multiarray_d(0, aero_coeffs)[0] = cl;
	get_col_Multiarray_d(0, aero_coeffs)[1] = cd;
	get_col_Multiarray_d(0, aero_coeffs)[2] = cm;

	return (const struct const_Multiarray_d*) aero_coeffs;
}


static const struct const_Multiarray_c* compute_cl_cd_cm_c(const struct Simulation* sim_c){

	int i, num_fc;

	// Get the freestream flow data and reference area
	const struct Sol_Data__c_dlm sol_data = get_sol_data__c_dlm();

	const double rho_fs   = sol_data.rho,
	             p_fs     = sol_data.p,
	             c_fs     = sqrt(GAMMA*p_fs/rho_fs),
	             V_fs     = sol_data.mach*c_fs,
	             denom    = 0.5*rho_fs*V_fs*V_fs*sol_data.area_ref,
	             theta_fs = sol_data.theta,
	             cm_le_x  = sol_data.cm_le_x,
	             cm_le_y  = sol_data.cm_le_y;

	// The force components along the Normal and Axial directions and the moment about the leading edge
	double complex	Force_N = 0.0,
					Force_A = 0.0,
					Moment  = 0.0;

	double complex delta_x, delta_y, x_fc, y_fc;

	for (struct Intrusive_Link* curr = sim_c->faces->first; curr; curr = curr->next) {
		// Loop through the faces that are on the wall boundary
		const struct Face*const face = (struct Face*) curr;
		const struct Solver_Face_c*const s_face = (struct Solver_Face_c*) curr;
		if (!is_face_wall_boundary(face))
			continue;

		// Get the solution (s) at the face cubature (fc) nodes and convert it to 
		// primitive variables (pv)
		const struct const_Multiarray_c* pv_fc_c = constructor_s_fc_interp_c(0, s_face);
		convert_variables_c((struct Multiarray_c*)pv_fc_c,'c','p');

		// Get the normals (n) at the face cubature (fc) nodes and jacobian values
		const struct const_Multiarray_c* n_fc 		= s_face->normals_fc;
		const struct const_Multiarray_c* jac_det_fc = s_face->jacobian_det_fc;

		// Get the xyz values of the face cubature (fc) nodes
		const struct const_Multiarray_c* xyz_fc 	= s_face->xyz_fc;

		num_fc = (int) pv_fc_c->extents[0]; // number of face cubature (fc) nodes

		// Pressure (p) at face cubature (fc) nodes as an array
		const complex double *p_fc_c_i = get_col_const_Multiarray_c(pv_fc_c->extents[1]-1,pv_fc_c);

		// Integrate the pressure multiplied by the negative of the normal component.
		// Compute the integral using quadrature
		const struct const_Vector_d*const w_fc = get_operator__w_fc__s_e_c(s_face);

		for (i = 0; i < num_fc; i++){
			const double complex *const n = get_row_const_Multiarray_c(i,n_fc);

			// Force Computations:
			Force_A += 1.0*p_fc_c_i[i]*n[0] * jac_det_fc->data[i] * w_fc->data[i];
			Force_N += 1.0*p_fc_c_i[i]*n[1] * jac_det_fc->data[i] * w_fc->data[i];

			// Moment Computations:
			x_fc = get_col_const_Multiarray_c(0, xyz_fc)[i];
			y_fc = get_col_const_Multiarray_c(1, xyz_fc)[i];

			delta_x = x_fc - cm_le_x;
			delta_y = y_fc - cm_le_y;

			Moment += (delta_x*p_fc_c_i[i]*n[1] - delta_y*p_fc_c_i[i]*n[0]) * jac_det_fc->data[i] * w_fc->data[i];
		}

		// Destroy the allocated vectors
		destructor_const_Multiarray_c(pv_fc_c);
	}

	// Get the lift and drag by taking into account the angle of attack 
	complex double 	L = Force_N*cos(theta_fs) - Force_A*sin(theta_fs),
					D = Force_N*sin(theta_fs) + Force_A*cos(theta_fs);

	complex double 	cl = L/denom,
					cd = D/denom,
					cm = -1.0*Moment/(denom*sol_data.area_ref);

 	struct Multiarray_c *aero_coeffs = constructor_empty_Multiarray_c('C',2,(ptrdiff_t[]){3, 1}); // destructed
	get_col_Multiarray_c(0, aero_coeffs)[0] = cl;
	get_col_Multiarray_c(0, aero_coeffs)[1] = cd;
	get_col_Multiarray_c(0, aero_coeffs)[2] = cm;

	return (const struct const_Multiarray_c*) aero_coeffs;
}


static struct Sol_Data__c_dlm get_sol_data__c_dlm(){
	static bool need_input = true;

	static struct Sol_Data__c_dlm sol_data;
	if (need_input) {
		need_input = false;
		read_data_c_dlm(&sol_data);
	}
	return sol_data;
}


static void read_data_c_dlm (struct Sol_Data__c_dlm*const sol_data){
	const int count_to_find = 7;
	int count_found = 0;
	char line[STRLEN_MAX];
	FILE* input_file = NULL;

	input_file = fopen_input('s',NULL,NULL); // closed
	while (fgets(line,sizeof(line),input_file)) {
		read_skip_string_count_d("density",  &count_found,line,&sol_data->rho);
		read_skip_string_count_d("pressure", &count_found,line,&sol_data->p);
		read_skip_string_count_d("mach",     &count_found,line,&sol_data->mach);
		read_skip_string_count_d("theta_deg",&count_found,line,&sol_data->theta);
	}
	fclose(input_file);

	sol_data->theta *= PI/180.0;

	input_file = fopen_input('g',NULL,NULL); // closed
	while (fgets(line,sizeof(line),input_file)) {
		read_skip_string_count_c_style_d("area_ref",&count_found,line,&sol_data->area_ref);
		read_skip_string_count_c_style_d("cm_le_x",&count_found,line,&sol_data->cm_le_x);
		read_skip_string_count_c_style_d("cm_le_y",&count_found,line,&sol_data->cm_le_y);
	}
	fclose(input_file);

	if (count_found != count_to_find)
		EXIT_ERROR("Did not find the required number of variables");
}


static struct Functional_Data get_functional_data(){

	// Make functional_data a static variable so that, once declared, it is saved in 
	// memory. Therefore, no need to read the file again in future function calls.
	static bool need_input = true;
	static struct Functional_Data functional_data;

	if(need_input){
		need_input = false;

		FILE* input_file = fopen_input('o',NULL,NULL); // closed
		char line[STRLEN_MAX];

		// Read in the information from the file
		while (fgets(line,sizeof(line),input_file)) {

			// Any line with a comment flag should be skipped
			if (strstr(line, "//"))
				continue;

			// Read specific functional data
			if (strstr(line, "FUNCTIONAL_TARGET_CL_VALUE")) read_skip_d_1(line, 1, &functional_data.target_cl, 1);
			if (strstr(line, "FUNCTIONAL_INVERSE_PRESSURE_DESIGN_TARGET_FILE_ABS_PATH")) read_skip_c(line, functional_data.target_pressure_distribution_file_abs_path);		
		}
		
		read_target_pressure_distribution_file(&functional_data);
		fclose(input_file);
	}

	return functional_data;
}


static void read_target_pressure_distribution_file(struct Functional_Data *functional_data){

	// If no file name was given, then there is no target pressure data to read
	if(strlen(functional_data->target_pressure_distribution_file_abs_path) == 0)
		return;

	// A valid file name was given. Read the data and load it into the Functional_Data struct
	FILE* fp;
	if ((fp = fopen(functional_data->target_pressure_distribution_file_abs_path,"r")) == NULL)
		printf("Error: Target Pressure Distribution file %s did not open.\n", 
			functional_data->target_pressure_distribution_file_abs_path), exit(1);


	char line[STRLEN_MAX];

	// 1) Get the number of pts
	int num_pts;
	fgets(line,sizeof(line),fp);
	sscanf(line, "%d", &num_pts);

	// 2) Read the point data 
	functional_data->target_pressure_distribution = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){num_pts,2});
	double 	*s = get_col_Multiarray_d(0, functional_data->target_pressure_distribution),
			*p = get_col_Multiarray_d(1, functional_data->target_pressure_distribution);


	for(int i = 0; i < num_pts; i++){
		fgets(line,sizeof(line),fp);
		sscanf(line, "%lf %lf", &s[i],&p[i]);
	}

	// Sort the pressure distribution data according to the knot domain location
	// Use an insertion sort

	int i = 1;
	int j;

	while (i < num_pts){
		j = i;
		while (j > 0 && (s[j-1] > s[j])){

			// Swap j and j-1
			double 	temp_s = s[j],
					temp_p = p[j];

			s[j] = s[j-1];
			p[j] = p[j-1];

			s[j-1] = temp_s;
			p[j-1] = temp_p;

			j--;
		}
		i++;
	}

	fclose(fp);
}


static double compute_volume_of_s_vol(const struct Solver_Volume* s_vol){


	struct Volume* vol = (struct Volume*)s_vol;
	const struct Solution_Element* s_e = &((struct Solver_Element*)vol->element)->s_e;

	const int curved = vol->curved,
	          p      = s_vol->p_ref;

	const struct const_Vector_d* w_vc = get_const_Multiarray_Vector_d(s_e->w_vc[curved],(ptrdiff_t[]){0,0,p,p});
	const struct const_Vector_d jacobian_det_vc = interpret_const_Multiarray_as_Vector_d(s_vol->jacobian_det_vc);

	const ptrdiff_t ext_0 = w_vc->ext_0;

	double volume = 0.0;
	for (int i = 0; i < ext_0; ++i)
		volume += w_vc->data[i]*jacobian_det_vc.data[i];

	return volume;

}


static double complex compute_volume_of_s_vol_c(const struct Solver_Volume_c* s_vol_c){
	

	struct Volume* vol = (struct Volume*)s_vol_c;
	const struct Solution_Element* s_e = &((struct Solver_Element*)vol->element)->s_e;

	const int curved = vol->curved,
	          p      = s_vol_c->p_ref;

	const struct const_Vector_d* w_vc = get_const_Multiarray_Vector_d(s_e->w_vc[curved],(ptrdiff_t[]){0,0,p,p});
	const struct const_Vector_c jacobian_det_vc = interpret_const_Multiarray_as_Vector_c(s_vol_c->jacobian_det_vc);

	const ptrdiff_t ext_0 = w_vc->ext_0;

	double complex volume = 0.0;
	for (int i = 0; i < ext_0; ++i)
		volume += w_vc->data[i]*jacobian_det_vc.data[i];

	return volume;

}


