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

#ifndef DPG__optimization_case_h__INCLUDED
#define DPG__optimization_case_h__INCLUDED

//#include "objective_functions.h"

#include "petscvec.h"
#include "petscmat.h"
#include "definitions_alloc.h"
#include "objective_functions.h"

struct Multiarray_d;
struct Matrix_d;

/** \brief Container for test case specific information.
 *
 *  This container is used to hold optimization case specific variables as well as function 
 *  pointers to various functions such that the control flow is not broken during run-time.
 *  
 */
struct Optimization_Case {

	/** The number of design control points degrees of freedom. For example, if there are two
	 * control points, one which can move in the x and y direction and the other that can only
	 * move in the y, the value of num_design_pts_dofs is 3.
	 */
	int num_design_pts_dofs;

	functional_fptr 	objective_function;  ///< Function pointer to the real objective function
	functional_fptr_c 	objective_function_c;  ///< Function pointer to the complex objective function

	struct Simulation *sim; ///< Pointer to the real sim object
	struct Simulation *sim_c; ///< Pointer to the complex sim object. Complex counterpart is needed for the complex step.

	/** Specifies the type of optimizer. Options are:
	 * LINE_SEARCH_STEEPEST_DESCENT = Uses gradient descent to minimize the objective function.
	 * LINE_SEARCH_BFGS = Uses the BFGS algorithm to minimize the objective function.
	 * NLPQLP = Uses NLPQLP to minimize the objective function
	 */
	char optimizer_spec[STRLEN_MIN];

	/** Specifies the type of optimization to take place. Options are:
	 * TARGET_CL = Optimize to attain a target lift coefficient
	 */
	char optimization_type_spec[STRLEN_MIN];

	// =================================
	//      Geometry Data Structures
	// =================================

	/// \brief Container for information relating to the geometry
	struct Geo_Data {

		int P, ///< Consult Geo_Data in geoemtry_parameteric_T.c 
			Q; ///< Consult Geo_Data in geoemtry_parameteric_T.c 

		struct Multiarray_d *knots_xi, ///< Consult Geo_Data in geoemtry_parameteric_T.c 
							*knots_eta; ///< Consult Geo_Data in geoemtry_parameteric_T.c 

		struct Multiarray_d *control_points; ///< Control points for the geometry
		struct Multiarray_c *control_points_c; ///< Complex version of the control points needed for the complex step

		/** Multiarray holding the control point weights. Stored as
		 * a multiarray of dimension [num_points x 1]
		 */
		struct Multiarray_d *control_weights;

		/** Multiarray holding the connectivity information for the control
		 * points (pt) and weights (wt) in the patch.
		 * The xi control points correspond to the i index and eta control 
		 * points to the j index. Data stored in row major form.
		 */
		struct Multiarray_i *control_pt_wt_connectivity;

		/** The multiarray used to store which design points are used for the optimization.
		 * Each row holds data of the form (control_pt_index, x dof bool, y dof bool). The
		 * bools are 0 if the point cannot be moved in the given coordinate direction and 
		 * 1 if it can.
		 */
		struct Multiarray_i *control_points_optimization;

		/** The limits for the dofs (ordering is the same as control_points_optimization)
		 * NOTE: For now, we only consider y displacements here.
		 * NOTE: Read this information from maybe the optimization.data file
		 */
		struct Multiarray_d *control_points_optimization_lims;

	} geo_data;	
};


/** \brief Create the optimization case structure. Fill it with data relevant to the optimization
 * and allocate memory (except for that used by the optimizer or adjoint solve)
 *
 * \return Optimization_Case data structure with the information for the optimization
 */
struct Optimization_Case* constructor_Optimization_Case (
	struct Simulation* sim ///< Standard. The simulation data structure
	);


/** \brief Free the optimization_case structure and clear allocated memory in the constructor
 */
void destructor_Optimization_Case (
	struct Optimization_Case* optimization_case ///< The optimization_case structure to be freed
	);


/** \brief 	Copy the data from the real sim object to the complex sim object. This method will 
 *	transfer all the face solver and volume solver data into the complex sim object.
 *
 *	NOTE: Since the same mesh and control file is read (no adaptation was done) we can 
 *		assume that the ordering of the volumes and faces in the sim and sim_c structure
 *		is the same.
 */
void copy_data_r_to_c_sim(
	struct Simulation* sim, ///< The sim data structure with real Volumes and Faces
	struct Simulation* sim_c ///< The sim data structure with Volumes and Faces that hold complex data.
	);

#endif // DPG__optimization_case_h__INCLUDED
