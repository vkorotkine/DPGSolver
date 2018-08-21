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
 *  \brief Provides the interface to templated functions used for parametric geometry processing.
 */

#include "def_templates_geometry.h"
#include "def_templates_volume_solver.h"
#include "def_templates_multiarray.h"

struct Solver_Volume_T;
struct Simulation;
struct const_Multiarray_i;

/** \brief Version of \ref constructor_xyz_fptr_T for the n-cube scaled and translated to an arbitrary fixed cube.
 *  \return See brief. */
const struct const_Multiarray_T* constructor_xyz_fixed_cube_parametric_T
	(const char n_type,                      ///< See brief.
	 const struct const_Multiarray_T* xyz_i, ///< See brief.
	 const struct Solver_Volume_T* s_vol,    ///< See brief.
	 const struct Simulation* sim            ///< See brief.
	);

/** \brief Version of \ref constructor_xyz_fptr_T for the parametric cylinder.
 *  \return See brief.
 *
 *  Uses \f$ r-\theta \f$ parametrization to transform from square to circular sections.
 */
const struct const_Multiarray_T* constructor_xyz_cylinder_parametric_T
	(const char n_type,                      ///< See brief.
	 const struct const_Multiarray_T* xyz_i, ///< See brief.
	 const struct Solver_Volume_T* s_vol,    ///< See brief.
	 const struct Simulation* sim            ///< See brief.
	);

/** \brief Version of \ref constructor_xyz_fptr_T for the n-cube perturbed with trigonometric functions.
 *  \return See brief.
 */
const struct const_Multiarray_T* constructor_xyz_trigonometric_cube_parametric_T
	(const char n_type,                      ///< See brief.
	 const struct const_Multiarray_T* xyz_i, ///< See brief.
	 const struct Solver_Volume_T* s_vol,    ///< See brief.
	 const struct Simulation* sim            ///< See brief.
	);

/** \brief Version of \ref constructor_xyz_trigonometric_cube_parametric_T where the only perturbed boundary is that
 *         corresponding to the leftmost x-coordinate.
 *  \return See brief.
 */
const struct const_Multiarray_T* constructor_xyz_trigonometric_cube_parametric_xl_T
	(const char n_type,                      ///< See brief.
	 const struct const_Multiarray_T* xyz_i, ///< See brief.
	 const struct Solver_Volume_T* s_vol,    ///< See brief.
	 const struct Simulation* sim            ///< See brief.
	);

/** \brief Version of \ref constructor_xyz_trigonometric_cube_parametric_xl_T with additional shift and scaling to use
 *         only the first octant.
 *  \return See brief. */
const struct const_Multiarray_T* constructor_xyz_trigonometric_cube_parametric_xl_oct1_T
	(const char n_type,                      ///< See brief.
	 const struct const_Multiarray_T* xyz_i, ///< See brief.
	 const struct Solver_Volume_T* s_vol,    ///< See brief.
	 const struct Simulation* sim            ///< See brief.
	);

/** \brief Version of \ref constructor_xyz_fptr_T for a parametric Joukowski to circle blended domain.
 *  \return See brief.
 *
 *  The procedure for the transformation is as follows:
 *  - Find the point corresponding to each input coordinate on the Joukowski and the external cylinder surfaces based on
 *    the reference \f$ r \f$ coordinate.
 *  - Find the final physical coordinate as a linearly blended contribution of the two surface coordinates using
 *    \f$ s \f$ as the blending parameter.
 */
const struct const_Multiarray_T* constructor_xyz_joukowski_parametric_T
	(const char n_type,                      ///< See brief.
	 const struct const_Multiarray_T* xyz_i, ///< See brief.
	 const struct Solver_Volume_T* s_vol,    ///< See brief.
	 const struct Simulation* sim            ///< See brief.
	);

/** \brief Version of \ref constructor_xyz_fptr_T for a parametric Gaussian bump to plane blended domain.
 *  \return See brief. */
const struct const_Multiarray_T* constructor_xyz_gaussian_bump_parametric_T
	(const char n_type,                      ///< See brief.
	 const struct const_Multiarray_T* xyz_i, ///< See brief.
	 const struct Solver_Volume_T* s_vol,    ///< See brief.
	 const struct Simulation* sim            ///< See brief.
	);

/** \brief Version of \ref constructor_xyz_fptr_T for a parametric NURBS patch domain.
 *  Computes the parametric NURBS mapping from the parametric domain to the physical domain.
 *  \return See brief. */
const struct const_Multiarray_T* constructor_xyz_NURBS_parametric_T
	(const char n_type,                      ///< See brief.
	 const struct const_Multiarray_T* xyz_i, ///< See brief.
	 const struct Solver_Volume_T* s_vol,    ///< See brief.
	 const struct Simulation* sim            ///< See brief.
	);

/** \brief Similar of \ref constructor_xyz_fptr_T. This function is unique in that it computes 
 *	the gradient terms of the parametric NURBS mapping from the NURBS parametric domain
 *	to the physical domain.
 *
 *  \return The multiarray of the gradient values using the exact NURBS mapping. 
 *	Consult \ref grad_xyz_NURBS_patch_mapping_T for details of multiarray ordering.
 */
const struct const_Multiarray_T* constructor_grad_xyz_NURBS_parametric_T
	(const char n_type, 						///< See brief.
	 const struct const_Multiarray_R* xyz_i,  	///< See brief.
	 const struct Solver_Volume_T* s_vol,  		///< See brief.
	 const struct Simulation* sim 				///< See brief.
	);

/** \brief 	Update the NURBS data held in the geo_data structure. This method
 *	is used with the optimization routines for updating the information 
 *	about the NURBS patch by adjusting the location of the control points to 
 *	allow for shape optimization to take place.
 */
void update_geo_data_NURBS_parametric_T
	(const struct const_Multiarray_T* control_points ///< The new values for the control point data
	);

/** \brief Computes the gradient terms of the NURBS mapping at the specified points on the 
 *	knot domain (xi_eta_i). 
 *
 * 	NOTE: Only the 2D case has been implemented for now.
 *	NOTE: (xi, eta) form the standard NURBS parametric domain variables (eta = 0 is the horizontal axis)
 *
 *	\return A multiarray of dimension [num_points x (DIM^2)] that holds the partial derivatives
 *	of the mapping. The ith row of the multiarray holds the partials at the ith point 
 *	and stores the data in the form [x_xi, y_xi, x_eta, y_eta], where x_xi, for instance, is the 
 *	partial of x with respect to the xi variable.
 */
const struct const_Multiarray_T *grad_xyz_NURBS_patch_mapping_T(

	/** The values on the parametric domain (knot domain). Is a multiarray
	 * of dimension [num_points x DIM]. For now, only consider 2D cases.
	 */
	const struct const_Multiarray_d* xi_eta_i,

	int P, 		///< The order of the basis functions in the xi direction
	int Q,		///< The order of the basis functions in the eta direction

	const struct const_Multiarray_d* knots_xi, 	///< The knot vector in the xi direction
	const struct const_Multiarray_d* knots_eta, ///< The knot vector in the eta direction

	/** The multiarray holding the list of control points. 
	 * The multiarray is of dimension [num_ctrl_pts x (DIM)]
	 */
	const struct const_Multiarray_T* control_points,
	
	/** The multiarray holding the list of weights associated to each point. 
	 * The multiarray is of dimension [num_ctrl_pts x 1]
	 */
	const struct const_Multiarray_d* control_weights, 
	
	/** The multiarray holding the connectivity information for
	 * the control points (pt) and weights (wt). Is of dimension [num_I x num_J], where num_I and
	 * num_J are the number of control points / basis functions in either coordinate 
	 * direction. This multiarray is in row major form
	 */
	const struct const_Multiarray_i* control_pt_wt_connectivity
	);


/** \brief Computes the mapped points from the knot domain onto the physical domain using the
 * 	NURBS mapping. 
 *
 * 	NOTE: Only the 2D case has been implemented for now.
 *	NOTE: (xi, eta) form the standard NURBS parametric domain variables (eta = 0 is the horizontal axis)
 *
 *	\return A multiarray of dimension [num_points x (DIM)] that holds the mapped values.
 */
const struct const_Multiarray_T *xyz_NURBS_patch_mapping_T(

	/** The values on the parametric domain (knot domain). Is a multiarray
	 * of dimension [num_points x DIM]. For now, only consider 2D cases.
	 */
	const struct const_Multiarray_d* xi_eta_i, 

	int P, 	///< The order of the basis functions in the xi direction
	int Q,	///< The order of the basis functions in the eta direction

	const struct const_Multiarray_d* knots_xi, 	///< The knot vector in the xi direction
	const struct const_Multiarray_d* knots_eta, ///< The knot vector in the eta direction

	/** The multiarray holding the list of control points. 
	 * The multiarray is of dimension [num_ctrl_pts x (DIM)]
	 */
	const struct const_Multiarray_T* control_points, 
	
	/** The multiarray holding the list of weights associated to each point. 
	 * The multiarray is of dimension [num_ctrl_pts x 1]
	 */
	const struct const_Multiarray_d* control_weights, 
	
	/** The multiarray holding the connectivity information for
	 * the control points (pt) and weights (wt). Is of dimension [num_I x num_J], where num_I and
	 * num_J are the number of control points / basis functions in either coordinate 
	 * direction. This multiarray is in row major form
	 */
	const struct const_Multiarray_i* control_pt_wt_connectivity 
	);


#include "undef_templates_geometry.h"
#include "undef_templates_volume_solver.h"
#include "undef_templates_multiarray.h"
