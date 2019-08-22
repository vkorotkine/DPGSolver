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

#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "gsl/gsl_math.h"

#include "macros.h"
#include "definitions_alloc.h"
#include "definitions_core.h"
#include "definitions_math.h"
#include "definitions_tol.h"

#include "def_templates_geometry_parametric.h"
#include "def_templates_volume_solver.h"
#include "def_templates_multiarray.h"
#include "def_templates_math_functions.h"
#include "simulation.h"

#include "bases.h"


// Static function declarations ************************************************************************************* //


struct NURBS_Patch_Geo_Data{
	// For multipatch, multiple NURBS patch structs need to be stored in Geo_Data.
	// This struct holds their information. 
	// TODO: Properly integrate with Manmeet's NURBS.  
	// NURBS Parametric Domain parameters: (2D Domain for now only)
	
	int P, ///< Order of the patch in the xi direction 
		Q;  ///< Order of the patch in the eta direction
	
	struct Multiarray_R *knots_xi, ///< Knot vector in the xi parametric direction
						*knots_eta; ///< Knot vector in the eta parametric direction

	/** Multiarray holding the control point data. Stored as
	 * a multiarray of dimension [num_points x DIM]. Is templated to allow
	 * for complex control points so that sensitivities can be calculated using the complex step.
	 */
	struct Multiarray_T *control_points;
	
	/** Multiarray holding the control point weights. Stored as
	 * a multiarray of dimension [num_points x 1]
	 */
	struct Multiarray_R *control_weights;

	/** Multiarray holding the connectivity information for the control
	 * points (pt) and weights (wt) in the patch.
	 * The xi control points correspond to the i index and eta control 
	 * points to the j index. Data stored in row major form.
	 */
	struct Multiarray_i *control_pt_wt_connectivity;
};

/// \brief Container for the geometric data required for the supported parametric cases.
struct Geo_Data {
	Real x_scale,  ///< Multiplicative scaling for the x-coordinates.
	     r_i,      ///< Radius corresponding to the 'i'nner surface.
	     r_o,      ///< Radius corresponding to the 'o'uter surface.
	     s_offset; ///< The offset for the reference domain 's'-coordinate.

	Real a,     ///< Geometry parameter.
	     b,     ///< Geometry parameter.
	     c,     ///< Geometry parameter.
	     d,     ///< Geometry parameter.
	     h,     ///< Geometry parameter ('h'eight generally).
	     x_max; ///< Maximum range for the x-coordinates.

	Real xyz_l, ///< Left  xyz-coordinate.
	     xyz_r; ///< Right xyz-coordinate.


	// NURBS Parametric Domain parameters: (2D Domain for now only)
	
	int P, ///< Order of the patch in the xi direction 
		Q;  ///< Order of the patch in the eta direction
	
	struct Multiarray_R *knots_xi, ///< Knot vector in the xi parametric direction
						*knots_eta; ///< Knot vector in the eta parametric direction

	/** Multiarray holding the control point data. Stored as
	 * a multiarray of dimension [num_points x DIM]. Is templated to allow
	 * for complex control points so that sensitivities can be calculated using the complex step.
	 */
	struct Multiarray_T *control_points;
	
	/** Multiarray holding the control point weights. Stored as
	 * a multiarray of dimension [num_points x 1]
	 */
	struct Multiarray_R *control_weights;

	/** Multiarray holding the connectivity information for the control
	 * points (pt) and weights (wt) in the patch.
	 * The xi control points correspond to the i index and eta control 
	 * points to the j index. Data stored in row major form.
	 */
	struct Multiarray_i *control_pt_wt_connectivity;

	struct NURBS_Patch_Geo_Data * NURBS_Patch_Data;
};


/** \brief Return the statically allocated \ref Geo_Data container.
 *  \return See brief. */
static struct Geo_Data get_geo_data
	(const char*const geo_name ///< The name corresponding to the type of geometry for which to obtain data.
	);

/** \brief Return the estimate of the definite integral of the input function over the input range computed using
 *         Romberg's method.
 *  \return See brief.
 *
 *  The code was copied from the [Wikipedia implementation of Romberg's method][wikipedia_romberg].
 *
 *  \note As curves being integrated here will likely always be \f$ C^{\infty} \f$, using a high-order quadrature may be
 *        more efficient. This can be tested if relevant after profiling.
 *
 *  <!-- References: -->
 *  [wikipedia_romberg]: https://en.wikipedia.org/wiki/Romberg%27s_method#Implementation
 */
static Real romberg
	(Real (*f)(const Real x, const struct Function_Data_GP*const f_data), ///< Function to integration.
	 const Real a,                                                        ///< Lower limit.
	 const Real b,                                                        ///< Upper limit.
	 const Real acc,                                                      ///< Accuracy tolerance.
	 const struct Function_Data_GP*const f_data                           ///< \ref Function_Data_GP.
	);

/** \brief Return the value of the arc length integrand for the Gaussian Bump function.
 *  \return See brief. */
static Real f_al_gaussian_bump
	(const Real x,                              ///< x-coordinate.
	 const struct Function_Data_GP*const f_data ///< \ref Function_Data_GP.
	);

// Interface functions ********************************************************************************************** //

const struct const_Multiarray_T* constructor_xyz_fixed_cube_parametric_T
	(const char n_type, const struct const_Multiarray_T* xyz_i, const struct Solver_Volume_T* s_vol,
	 const struct Simulation* sim)
{
	UNUSED(n_type);
	UNUSED(s_vol);
	UNUSED(sim);
	assert(DIM == xyz_i->extents[1]);

	const ptrdiff_t n_n = xyz_i->extents[0];

	struct Multiarray_T* xyz = constructor_empty_Multiarray_T('C',2,(ptrdiff_t[]){n_n,DIM}); // returned

	const Type*const x_i = get_col_const_Multiarray_T(0,xyz_i),
	          *const y_i = ( DIM > 1 ? get_col_const_Multiarray_T(1,xyz_i) : NULL ),
	          *const z_i = ( DIM > 2 ? get_col_const_Multiarray_T(2,xyz_i) : NULL );

	Type*const x = get_col_Multiarray_T(0,xyz),
	    *const y = ( DIM > 1 ? get_col_Multiarray_T(1,xyz) : NULL ),
	    *const z = ( DIM > 2 ? get_col_Multiarray_T(2,xyz) : NULL );

	struct Geo_Data geo_data = get_geo_data("fixed_cube");

	const Real xyz_l = geo_data.xyz_l,
	           xyz_r = geo_data.xyz_r;
	const Real dxyz  = xyz_r-xyz_l;

	for (int n = 0; n < n_n; ++n) {
		const double dx = 0.5*(real_T(x_i[n])+1.0);
		x[n] = xyz_l+dx*dxyz;

		if (DIM > 1)
			EXIT_ADD_SUPPORT; UNUSED(y); UNUSED(z); UNUSED(y_i); UNUSED(z_i);
	}
	return (struct const_Multiarray_T*) xyz;
}

const struct const_Multiarray_T* constructor_xyz_cylinder_parametric_T
	(const char n_type, const struct const_Multiarray_T* xyz_i, const struct Solver_Volume_T* s_vol,
	 const struct Simulation* sim)
{
	UNUSED(n_type);
	UNUSED(s_vol);
	UNUSED(sim);
	assert(DIM >= 2);
	assert(DIM == xyz_i->extents[1]);

	const ptrdiff_t n_n = xyz_i->extents[0];

	struct Multiarray_T* xyz = constructor_empty_Multiarray_T('C',2,(ptrdiff_t[]){n_n,DIM}); // returned

	const Type*const x_i = get_col_const_Multiarray_T(0,xyz_i),
	          *const y_i = get_col_const_Multiarray_T(1,xyz_i),
	          *const z_i = ( DIM > 2 ? get_col_const_Multiarray_T(2,xyz_i) : NULL );

	Type*const x = get_col_Multiarray_T(0,xyz),
	    *const y = get_col_Multiarray_T(1,xyz),
	    *const z = ( DIM > 2 ? get_col_Multiarray_T(2,xyz) : NULL );

	for (int n = 0; n < n_n; ++n) {
		const Real xy[2] = { real_T(x_i[n]), real_T(y_i[n]), };

		const Real r = GSL_MAX(fabs(xy[0]),fabs(xy[1]));

		Real t = atan2(xy[1],xy[0]);
		if (t >= -1.0*PI_OVER_4 && t < 1.0*PI_OVER_4)
			t = xy[1]/r*PI_OVER_4;
		else if ( t >=  1.0*PI_OVER_4 && t < 3.0*PI_OVER_4)
			t = 0.5*PI - xy[0]/r*PI_OVER_4;
		else if ((t >=  3.0*PI_OVER_4 && t <= 4.0*PI_OVER_4) || (t >= -4.0*PI_OVER_4 && t < -3.0*PI_OVER_4))
			t = PI - xy[1]/r*PI_OVER_4;
		else if (t >= -3.0*PI_OVER_4 && t < -1.0*PI_OVER_4)
			t = 1.5*PI + xy[0]/r*PI_OVER_4;
		else
			EXIT_ERROR("Unsupported: %f %f %f\n",xy[0],xy[1],t);

		x[n] = r*cos(t);
		y[n] = r*sin(t);
		if (DIM > 2)
			z[n] = z_i[n];
	}
	return (struct const_Multiarray_T*) xyz;
}

const struct const_Multiarray_T* constructor_xyz_trigonometric_cube_parametric_T
	(const char n_type, const struct const_Multiarray_T* xyz_i, const struct Solver_Volume_T* s_vol,
	 const struct Simulation* sim)
{
	UNUSED(n_type);
	UNUSED(s_vol);
	UNUSED(sim);
	assert(DIM >= 2);
	assert(DIM == xyz_i->extents[1]);

	const ptrdiff_t n_n = xyz_i->extents[0];

	struct Multiarray_T* xyz = constructor_empty_Multiarray_T('C',2,(ptrdiff_t[]){n_n,DIM}); // returned

	const Type*const x_i = get_col_const_Multiarray_T(0,xyz_i),
	          *const y_i = ( DIM > 1 ? get_col_const_Multiarray_T(1,xyz_i) : NULL ),
	          *const z_i = ( DIM > 2 ? get_col_const_Multiarray_T(2,xyz_i) : NULL );

	Type*const x = get_col_Multiarray_T(0,xyz),
	    *const y = ( DIM > 1 ? get_col_Multiarray_T(1,xyz) : NULL ),
	    *const z = ( DIM > 2 ? get_col_Multiarray_T(2,xyz) : NULL );

	const Real dxyz = 0.05;
	for (int n = 0; n < n_n; ++n) {
		if (DIM == 2) {
			x[n] = x_i[n] + dxyz*sin(PI*real_T(y_i[n]));
			y[n] = y_i[n];
		} else if (DIM == 3) {
			EXIT_ADD_SUPPORT; UNUSED(z_i); UNUSED(z);
		} else {
			EXIT_UNSUPPORTED;
		}
	}
	return (struct const_Multiarray_T*) xyz;
}

const struct const_Multiarray_T* constructor_xyz_trigonometric_cube_parametric_xl_T
	(const char n_type, const struct const_Multiarray_T* xyz_i, const struct Solver_Volume_T* s_vol,
	 const struct Simulation* sim)
{
	UNUSED(n_type);
	UNUSED(s_vol);
	UNUSED(sim);
	assert(DIM >= 2);
	assert(DIM == xyz_i->extents[1]);

	const ptrdiff_t n_n = xyz_i->extents[0];

	struct Multiarray_T* xyz = constructor_empty_Multiarray_T('C',2,(ptrdiff_t[]){n_n,DIM}); // returned

	const Type*const x_i = get_col_const_Multiarray_T(0,xyz_i),
	          *const y_i = ( DIM > 1 ? get_col_const_Multiarray_T(1,xyz_i) : NULL ),
	          *const z_i = ( DIM > 2 ? get_col_const_Multiarray_T(2,xyz_i) : NULL );

	Type*const x = get_col_Multiarray_T(0,xyz),
	    *const y = ( DIM > 1 ? get_col_Multiarray_T(1,xyz) : NULL ),
	    *const z = ( DIM > 2 ? get_col_Multiarray_T(2,xyz) : NULL );

	const Real dxyz = 0.15;
	for (int n = 0; n < n_n; ++n) {
		const Real x_r = real_T(x_i[n]),
		           y_r = real_T(y_i[n]);
		if (DIM == 2) {
			x[n] = x_i[n] + (1.0-x_r)/2.0*dxyz*cos(PI*x_r+0.2)*cos(0.5*PI*y_r+0.3)*sin(PI*y_r);
			y[n] = y_i[n] + 2.0*dxyz*cos(PI/2.0*y_r)*sin(0.25*PI*y_r-0.1)*cos(PI/2.0*x_r-0.2);
		} else if (DIM == 3) {
			EXIT_ADD_SUPPORT; UNUSED(z_i); UNUSED(z);
		} else {
			EXIT_UNSUPPORTED;
		}
	}
	return (struct const_Multiarray_T*) xyz;
}

const struct const_Multiarray_T* constructor_xyz_trigonometric_cube_parametric_xl_oct1_T
	(const char n_type, const struct const_Multiarray_T* xyz_i, const struct Solver_Volume_T* s_vol,
	 const struct Simulation* sim)
{
	struct Multiarray_T*const xyz = (struct Multiarray_T*)
		constructor_xyz_trigonometric_cube_parametric_xl_T(n_type,xyz_i,s_vol,sim); // returned

	Type*const xyz_a[DIM] = ARRAY_DIM( get_col_Multiarray_T(0,xyz),
	                                   get_col_Multiarray_T(1,xyz),
	                                   get_col_Multiarray_T(2,xyz) );

	const ptrdiff_t n_n = xyz->extents[0];
	for (int n = 0; n < n_n; ++n) {
	for (int d = 0; d < DIM; ++d) {
		xyz_a[d][n] = 0.5*(xyz_a[d][n]+1.0)+0.5;
	}}
	return (struct const_Multiarray_T*) xyz;
}

const struct const_Multiarray_T* constructor_xyz_joukowski_parametric_T
	(const char n_type, const struct const_Multiarray_T* xyz_i, const struct Solver_Volume_T* s_vol,
	 const struct Simulation* sim)
{
	UNUSED(n_type);
	UNUSED(s_vol);
	UNUSED(sim);
	assert(DIM >= 2);
	assert(DIM == xyz_i->extents[1]);

	const ptrdiff_t n_n = xyz_i->extents[0];

	struct Multiarray_T* xyz = constructor_empty_Multiarray_T('C',2,(ptrdiff_t[]){n_n,DIM}); // returned

	const Type*const x_i = get_col_const_Multiarray_T(0,xyz_i),
	          *const y_i = get_col_const_Multiarray_T(1,xyz_i),
	          *const z_i = ( DIM > 2 ? get_col_const_Multiarray_T(2,xyz_i) : NULL );

	Type*const x = get_col_Multiarray_T(0,xyz),
	    *const y = get_col_Multiarray_T(1,xyz),
	    *const z = ( DIM > 2 ? get_col_Multiarray_T(2,xyz) : NULL );

	struct Geo_Data geo_data = get_geo_data("joukowski");

	const Real x_scale  = geo_data.x_scale,
	           r_i      = geo_data.r_i,
	           r_o      = geo_data.r_o,
	           s_offset = geo_data.s_offset,
	           center_x = 1.0-r_i; // The Joukowski transformation requires that the point (0,1) is included.

	const Real center_x_c  = -x_scale,
	           zeta_tail   = center_x+r_i, // theta = 0.0 in equation for zeta below.
	           x_mid_p     = (zeta_tail+1.0/zeta_tail)*x_scale,
	           y_mid_p     = sqrt(pow(r_o,2.0)-pow(x_mid_p-center_x_c,2.0)),
	           min_radians = atan2(y_mid_p,x_mid_p-center_x_c);

	for (int n = 0; n < n_n; ++n) {
		assert(y_i[n] != 0.0);

		const int sign_y = ( real_T(y_i[n]) > 0 ? 1 : -1 );
		const Real r = real_T(x_i[n]),
		           s = fabs(real_T(y_i[n]))-s_offset;

		if (r <= 0.0) {
			const Real b_r[] = { -r, 1.0+r, },
			           b_s[] = { 1.0-s, s, };

			// Joukowski portion
			const Real t_j = PI*(-r);

			const Complex zeta = center_x+r_i*cos(t_j) + I*(r_i*sin(t_j));
			const Complex z_j = zeta + 1.0/zeta;

			const Real x_j = creal(z_j)*x_scale,
			           y_j = cimag(z_j);

			// Cylinder portion
			const Real t_c = PI*b_r[0]+min_radians*b_r[1],
			           x_c = r_o*cos(t_c) + center_x_c,
			           y_c = r_o*sin(t_c);

			// Blended contribution
			x[n] =         x_j*b_s[0] + x_c*b_s[1];
			y[n] = sign_y*(y_j*b_s[0] + y_c*b_s[1]);
		} else {
			const Real b_r[] = { 1.0-r, r, },
			           b_s[] = { 1.0-s, s, };
			x[n] =         x_mid_p*b_r[0] + (x_mid_p+r_o)*b_r[1];
			y[n] = sign_y*(0.0    *b_s[0] +  y_mid_p*b_s[1]);
		}
		if (DIM > 2)
			z[n] = z_i[n];
	}
	return (struct const_Multiarray_T*) xyz;
}

void update_geo_data_NURBS_parametric_T(const struct const_Multiarray_T* control_points){

	struct Geo_Data geo_data = get_geo_data("NURBS");  // get static geo_data struct
	struct Multiarray_T *dest = geo_data.control_points;

	for (int i = 0; i < (int)dest->extents[0]; i++){
		for(int j = 0; j < (int)dest->extents[1]; j++){
			get_col_Multiarray_T(j, dest)[i] = get_col_const_Multiarray_T(j, control_points)[i];
		}
	}
}


const struct const_Multiarray_T* constructor_grad_xyz_NURBS_parametric_T
(const char n_type, const struct const_Multiarray_R* xyz_i, const struct Solver_Volume_T* s_vol,
	const struct Simulation* sim){

	UNUSED(n_type);
	UNUSED(s_vol);

	// Read the geometric data for the NURBS patch
	struct Geo_Data geo_data;
	struct NURBS_Patch_Geo_Data geo_patch_data;
	struct Multiarray_R *xyz_i_real = constructor_empty_Multiarray_R('C',2,xyz_i->extents);  // free
	
	int k = 0;
	for (int i = 0; i < xyz_i->extents[0]; i++){
		for (int j = 0; j < xyz_i->extents[1]; j++){
				
#if TYPE_RC == TYPE_REAL
			xyz_i_real->data[k] = xyz_i->data[k];
#elif TYPE_RC == TYPE_COMPLEX
			xyz_i_real->data[k] = creal(xyz_i->data[k]);
#endif

			k++;
		}
	}

	const Real*const x_i = get_col_Multiarray_d(0,xyz_i_real),
	          *const y_i = get_col_Multiarray_d(1,xyz_i_real);

	
	if (sim->nurbs_multipatch){
		geo_data = get_geo_data("Multipatch");
		int patch_index;
		if (*x_i>0. && *y_i<=0.){
			patch_index=0;
		}
		else if (*x_i<=0. && *y_i<=0.){
			patch_index=1;
		}
		else if (*x_i<=0. && *y_i>0.){
			patch_index=2;
		}
		else if (*x_i>0. && *y_i>0.){
			patch_index=3;
		};
		
		geo_patch_data=geo_data.NURBS_Patch_Data[patch_index];
	}
	else{
		geo_data = get_geo_data("NURBS");
	}
		
	const struct const_Multiarray_T *grad_xyz;
	if (sim->nurbs_multipatch){
	grad_xyz = grad_xyz_NURBS_patch_mapping_T(
		(const struct const_Multiarray_d*)xyz_i, geo_patch_data.P, geo_patch_data.Q, 
		(const struct const_Multiarray_d*)geo_patch_data.knots_xi, 
		(const struct const_Multiarray_d*)geo_patch_data.knots_eta,
		(const struct const_Multiarray_T*)geo_patch_data.control_points,
		(const struct const_Multiarray_d*)geo_patch_data.control_weights,
		(const struct const_Multiarray_i*)geo_patch_data.control_pt_wt_connectivity);
	}else{
	// Compute the gradient values at the desired points on the knot domain (xyz_i) using the mapping
	grad_xyz = grad_xyz_NURBS_patch_mapping_T(
		(const struct const_Multiarray_d*)xyz_i, geo_data.P, geo_data.Q, 
		(const struct const_Multiarray_d*)geo_data.knots_xi, 
		(const struct const_Multiarray_d*)geo_data.knots_eta,
		(const struct const_Multiarray_T*)geo_data.control_points,
		(const struct const_Multiarray_d*)geo_data.control_weights,
		(const struct const_Multiarray_i*)geo_data.control_pt_wt_connectivity);
	}


	return (const struct const_Multiarray_T*)grad_xyz;
}

const struct const_Multiarray_T* constructor_xyz_Multipatch_parametric_T
(const char n_type, const struct const_Multiarray_T* xyz_i, const struct Solver_Volume_T* s_vol,
	const struct Simulation* sim){

	UNUSED(n_type);
	UNUSED(s_vol);
	UNUSED(sim);

	assert(DIM == 2); // Add support for 3D if required.
	assert(DIM == xyz_i->extents[1]);

	// Read the geometric data for the NURBS patch
	const struct Geo_Data geo_data = get_geo_data("Multipatch");

	// xyz_i are the initial values for the xyz coordinates. These are, in 
	// this function, the location on the NURBS parametric domain (or knot domain). Map these
	// values onto the physical domain
 
	// NOTE: For NURBS mappings, the coordinates on the parametric domain (xyz_i) will need to be real
	struct Multiarray_R *xyz_i_real = constructor_empty_Multiarray_R('C',2,xyz_i->extents);  // free
	
	int k = 0;
	for (int i = 0; i < xyz_i->extents[0]; i++){
		for (int j = 0; j < xyz_i->extents[1]; j++){
				
#if TYPE_RC == TYPE_REAL
			xyz_i_real->data[k] = xyz_i->data[k];
#elif TYPE_RC == TYPE_COMPLEX
			xyz_i_real->data[k] = creal(xyz_i->data[k]);
#endif

			k++;
		}
	}

	const Real*const x_i = get_col_Multiarray_d(0,xyz_i_real),
	          *const y_i = get_col_Multiarray_d(1,xyz_i_real);


	int patch_index;
	if (*x_i>0. && *y_i<=0.){
		patch_index=0;
	}
	else if (*x_i<=0. && *y_i<=0.){
		patch_index=1;
	}
	else if (*x_i<=0. && *y_i>0.){
		patch_index=2;
	}
	else if (*x_i>0. && *y_i>0.){
		patch_index=3;
	};
	
	struct NURBS_Patch_Geo_Data geo_patch_data=geo_data.NURBS_Patch_Data[patch_index];
	const struct const_Multiarray_T *xyz = xyz_NURBS_patch_mapping_T(
		(const struct const_Multiarray_d*)xyz_i_real, geo_patch_data.P, geo_patch_data.Q, 
		(const struct const_Multiarray_d*)geo_patch_data.knots_xi, 
		(const struct const_Multiarray_d*)geo_patch_data.knots_eta,
		(const struct const_Multiarray_T*)geo_patch_data.control_points,
		(const struct const_Multiarray_d*)geo_patch_data.control_weights,
		(const struct const_Multiarray_i*)geo_patch_data.control_pt_wt_connectivity);

	destructor_Multiarray_d(xyz_i_real);
	
	return (const struct const_Multiarray_T*)xyz;
}

const struct const_Multiarray_T* constructor_xyz_NURBS_parametric_T
(const char n_type, const struct const_Multiarray_T* xyz_i, const struct Solver_Volume_T* s_vol,
	const struct Simulation* sim){

	UNUSED(n_type);
	UNUSED(s_vol);
	UNUSED(sim);

	assert(DIM == 2); // Add support for 3D if required.
	assert(DIM == xyz_i->extents[1]);

	// Read the geometric data for the NURBS patch
	const struct Geo_Data geo_data = get_geo_data("NURBS");

	// xyz_i are the initial values for the xyz coordinates. These are, in 
	// this function, the location on the NURBS parametric domain (or knot domain). Map these
	// values onto the physical domain
 
	// NOTE: For NURBS mappings, the coordinates on the parametric domain (xyz_i) will need to be real
	struct Multiarray_R *xyz_i_real = constructor_empty_Multiarray_R('C',2,xyz_i->extents);  // free
	
	int k = 0;
	for (int i = 0; i < xyz_i->extents[0]; i++){
		for (int j = 0; j < xyz_i->extents[1]; j++){
				
#if TYPE_RC == TYPE_REAL
			xyz_i_real->data[k] = xyz_i->data[k];
#elif TYPE_RC == TYPE_COMPLEX
			xyz_i_real->data[k] = creal(xyz_i->data[k]);
#endif

			k++;
		}
	}

	const struct const_Multiarray_T *xyz = xyz_NURBS_patch_mapping_T(
		(const struct const_Multiarray_d*)xyz_i_real, geo_data.P, geo_data.Q, 
		(const struct const_Multiarray_d*)geo_data.knots_xi, 
		(const struct const_Multiarray_d*)geo_data.knots_eta,
		(const struct const_Multiarray_T*)geo_data.control_points,
		(const struct const_Multiarray_d*)geo_data.control_weights,
		(const struct const_Multiarray_i*)geo_data.control_pt_wt_connectivity);

	destructor_Multiarray_d(xyz_i_real);
	
	return (const struct const_Multiarray_T*)xyz;
}
const struct const_Multiarray_T* constructor_xyz_gaussian_bump_parametric_T
	(const char n_type, const struct const_Multiarray_T* xyz_i, const struct Solver_Volume_T* s_vol,
	 const struct Simulation* sim)
{
	UNUSED(n_type);
	UNUSED(s_vol);
	UNUSED(sim);
	assert(DIM == 2); // Add support for 3D if required.
	assert(DIM == xyz_i->extents[1]);

	const ptrdiff_t n_n = xyz_i->extents[0];

	struct Multiarray_T* xyz = constructor_empty_Multiarray_T('C',2,(ptrdiff_t[]){n_n,DIM}); // returned

	const Type*const x_i = get_col_const_Multiarray_T(0,xyz_i),
	          *const y_i = get_col_const_Multiarray_T(1,xyz_i);

	Type*const x = get_col_Multiarray_T(0,xyz),
	    *const y = get_col_Multiarray_T(1,xyz);

	const struct Geo_Data geo_data = get_geo_data("gaussian_bump");
	const Real h     = geo_data.h,
	           x_max = geo_data.x_max;

	struct Function_Data_GP f_data = { .scale = 1.0, };

	const Real tol = 1e2*EPS;
	const Real x_l     = -x_max,
	           x_total = 2.0*x_max;
	const Real al_total = romberg(f_al_gaussian_bump,x_l,x_max,tol,&f_data);

	for (int n = 0; n < n_n; ++n) {
		const Real r = real_T(x_i[n]),
		           s = real_T(y_i[n]);

		// Plane contribution
		const Real x_p = x_max*r,
		           y_p = h;

		// Gaussian bump contribution
		Real x_g = x_p;

		int count = 0;
		for (Real newton_f = 1.0; fabs(newton_f) > tol; ) {
			const Real al = romberg(f_al_gaussian_bump,x_l,x_g,tol,&f_data);
			           newton_f     = 1.0/al_total*al - (x_p-x_l)/x_total;
			const Real newton_df_dx = 1.0/al_total*f_al_gaussian_bump(x_g,&f_data);
			x_g -= newton_f/newton_df_dx;

			++count;
			if (count > 10)
				EXIT_ERROR("Newton's method not converging. Error: % .3e.",newton_f);
		}
		const Real y_g = f_gaussian_bump(x_g,0,&f_data);

		// Blended contribution
		const Real b_s[] = { 0.5*(1.0-s), 0.5*(1.0+s), };
		x[n] = x_g*b_s[0] + x_p*b_s[1];
		y[n] = y_g*b_s[0] + y_p*b_s[1];
	}
	return (struct const_Multiarray_T*) xyz;
}


const struct const_Multiarray_T *grad_xyz_NURBS_patch_mapping_T(
	const struct const_Multiarray_d* xi_eta_i, int P, int Q, 
	const struct const_Multiarray_d* knots_xi, 
	const struct const_Multiarray_d* knots_eta, 
	const struct const_Multiarray_T* control_points,
	const struct const_Multiarray_d* control_weights,
	const struct const_Multiarray_i* control_pt_wt_connectivity){

	// Preprocessing:

	// Create multiarray structures to hold the control points and weights. The
	// control point structure will be of dimension [num_I x num_J] and there will be two
	// separate structures for the x and y values. The weight structure is of order 
	// [num_I x num_J].

	int i, j, k, numI, numJ, control_pt_index;
	numI = (int) control_pt_wt_connectivity->extents[0];
	numJ = (int) control_pt_wt_connectivity->extents[1];

	struct Multiarray_T *control_pts_x = constructor_empty_Multiarray_T('C',2,(ptrdiff_t[]){numI, numJ});  // free
	struct Multiarray_T *control_pts_y = constructor_empty_Multiarray_T('C',2,(ptrdiff_t[]){numI, numJ});  // free
	
	struct Multiarray_d *control_pts_w = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){numI, numJ});  // free

	const Type *const control_pts_x_list = get_col_const_Multiarray_T(0, control_points);
	const Type *const control_pts_y_list = get_col_const_Multiarray_T(1, control_points);
	
	const double *const control_pts_w_list = get_col_const_Multiarray_d(0, control_weights);
	
	for (i = 0; i < numI; i++){
		for (j = 0; j < numJ; j++){

			// The index of the control point in the list
			control_pt_index = get_row_const_Multiarray_i(i, control_pt_wt_connectivity)[j]; 

			get_col_Multiarray_T(j, control_pts_x)[i] = control_pts_x_list[control_pt_index];
			get_col_Multiarray_T(j, control_pts_y)[i] = control_pts_y_list[control_pt_index];
			
			get_col_Multiarray_d(j, control_pts_w)[i] = control_pts_w_list[control_pt_index];

		}
	}

	// Compute the partial derivative terms:
	int num_pts = (int) xi_eta_i->extents[0];
	struct Multiarray_T *grad_xyz = constructor_empty_Multiarray_T('C',2,(ptrdiff_t[]){num_pts, 4});  // returned

	// The values on the parametric domain to obtain partial derivatives at (xi and eta)
	struct Multiarray_d *xi_vals = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){num_pts,1});  // free
	struct Multiarray_d *eta_vals = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){num_pts,1});  // free
	
	for (i = 0; i < num_pts; i++){
		get_col_Multiarray_d(0, xi_vals)[i]  = get_col_const_Multiarray_d(0, xi_eta_i)[i];
		get_col_Multiarray_d(0, eta_vals)[i] = get_col_const_Multiarray_d(1, xi_eta_i)[i];
	}

	// Get the B Spline Basis (and derivative) function values (1D) at the given xi points
	const struct const_Multiarray_d *N_p_xi  = B_Spline_Basis_p(P, (const struct const_Multiarray_d*) xi_vals, knots_xi);  // free
	const struct const_Multiarray_d *dN_p_xi = derivative_B_Spline_Basis_p(P, (const struct const_Multiarray_d*) xi_vals, knots_xi);  // free

	// Get the B Spline Basis (and derivative) function values (1D) at the given eta points
	const struct const_Multiarray_d *N_q_eta  = B_Spline_Basis_p(Q, (const struct const_Multiarray_d*) eta_vals, knots_eta);  // free
	const struct const_Multiarray_d *dN_q_eta = derivative_B_Spline_Basis_p(Q, (const struct const_Multiarray_d*) eta_vals, knots_eta);  // free

	const struct const_Multiarray_d* grad_vals;
	grad_vals = grad_NURBS_basis_pq(N_p_xi, dN_p_xi, N_q_eta, dN_q_eta, (const struct const_Multiarray_d*)control_pts_w);

	const double *const grad_vals_i = grad_vals->data;
	double R_ij_xi, R_ij_eta;  // The partials of the NURBS basis function

	Type x_xi, y_xi, x_eta, y_eta;

	for (k = 0; k < num_pts; k++){
		// Loop over all the (xi,eta) points

		x_xi = 0;
		y_xi = 0;
		x_eta = 0;
		y_eta = 0;

		for (i = 0; i < numI; i++){
			for (j = 0; j < numJ; j++){

				// Get the i,j NURBS basis function gradient
				R_ij_xi  = grad_vals_i[i + j*numI + k*numI*numJ + 0*numI*numJ*num_pts];
				R_ij_eta = grad_vals_i[i + j*numI + k*numI*numJ + 1*numI*numJ*num_pts];

				// partial with respect to xi:
				x_xi += R_ij_xi * get_col_Multiarray_T(j, control_pts_x)[i];
				y_xi += R_ij_xi * get_col_Multiarray_T(j, control_pts_y)[i];

				// partial with respect to eta:
				x_eta += R_ij_eta * get_col_Multiarray_T(j, control_pts_x)[i];
				y_eta += R_ij_eta * get_col_Multiarray_T(j, control_pts_y)[i];
			}
		}

		// Store the partial terms
		get_col_Multiarray_T(0, grad_xyz)[k] = x_xi;
		get_col_Multiarray_T(1, grad_xyz)[k] = y_xi;
		get_col_Multiarray_T(2, grad_xyz)[k] = x_eta;
		get_col_Multiarray_T(3, grad_xyz)[k] = y_eta;

	}

	// Free allocated arrays
	destructor_Multiarray_T(control_pts_x);
	destructor_Multiarray_T(control_pts_y);
	destructor_Multiarray_d(control_pts_w);
	destructor_Multiarray_d(xi_vals);
	destructor_Multiarray_d(eta_vals);
	destructor_const_Multiarray_d(N_p_xi);
	destructor_const_Multiarray_d(dN_p_xi);
	destructor_const_Multiarray_d(N_q_eta);	
	destructor_const_Multiarray_d(dN_q_eta);
	destructor_const_Multiarray_d(grad_vals);

	return (const struct const_Multiarray_T*)grad_xyz;
}


const struct const_Multiarray_T *xyz_NURBS_patch_mapping_T(
	const struct const_Multiarray_d* xi_eta_i, int P, int Q, 
	const struct const_Multiarray_d* knots_xi, 
	const struct const_Multiarray_d* knots_eta, 
	const struct const_Multiarray_T* control_points,
	const struct const_Multiarray_d* control_weights,
	const struct const_Multiarray_i* control_pt_wt_connectivity){

	// Preprocessing:

	// Create multiarray structures to hold the control points and weights. The
	// control point structure will be of dimension [num_I x num_J] and there will be two
	// separate structures for the x and y values. The weight structure is of order 
	// [num_I x num_J].

	int i, j, k, numI, numJ, control_pt_index;
	numI = (int) control_pt_wt_connectivity->extents[0];
	numJ = (int) control_pt_wt_connectivity->extents[1];

	struct Multiarray_T *control_pts_x = constructor_empty_Multiarray_T('C',2,(ptrdiff_t[]){numI, numJ});  // free
	struct Multiarray_T *control_pts_y = constructor_empty_Multiarray_T('C',2,(ptrdiff_t[]){numI, numJ});  // free
	
	const Type *const control_pts_x_list = get_col_const_Multiarray_T(0, control_points);
	const Type *const control_pts_y_list = get_col_const_Multiarray_T(1, control_points);
	
	struct Multiarray_d *control_pts_w = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){numI, numJ});  // free
	
	const double *const control_pts_w_list = get_col_const_Multiarray_d(0, control_weights);
	
	for (i = 0; i < numI; i++){
		for (j = 0; j < numJ; j++){

			// The index of the control point in the list
			control_pt_index = get_row_const_Multiarray_i(i, control_pt_wt_connectivity)[j]; 

			get_col_Multiarray_T(j, control_pts_x)[i] = control_pts_x_list[control_pt_index];
			get_col_Multiarray_T(j, control_pts_y)[i] = control_pts_y_list[control_pt_index];
			
			get_col_Multiarray_d(j, control_pts_w)[i] = control_pts_w_list[control_pt_index];
		}
	}

	// Perform the Mapping:

	int num_pts = (int) xi_eta_i->extents[0];
	struct Multiarray_T *xyz = constructor_empty_Multiarray_T('C',2, (ptrdiff_t[]){num_pts, 2});  // return

	// The values on the parametric domain to obtain partial derivatives at (xi and eta)
	struct Multiarray_d *xi_vals  = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){num_pts,1});  // free
	struct Multiarray_d *eta_vals = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){num_pts,1});  // free
	
	for (i = 0; i < num_pts; i++){
		get_col_Multiarray_d(0, xi_vals)[i]  = get_col_const_Multiarray_d(0, xi_eta_i)[i];
		get_col_Multiarray_d(0, eta_vals)[i] = get_col_const_Multiarray_d(1, xi_eta_i)[i];
	}

	// Get the B Spline Basis function values (1D) at the given xi points
	const struct const_Multiarray_d *N_p_xi = B_Spline_Basis_p(P, (const struct const_Multiarray_d*) xi_vals, knots_xi);  // free
	
	// Get the B Spline Basis function values (1D) at the given eta points
	const struct const_Multiarray_d *N_q_eta = B_Spline_Basis_p(Q, (const struct const_Multiarray_d*) eta_vals, knots_eta);  // free
	
	// Get the NURBS basis functions evaluated at the provided points
	const struct const_Multiarray_d* basis_vals;
	basis_vals = NURBS_basis_pq(N_p_xi, N_q_eta, (const struct const_Multiarray_d*)control_pts_w);  // free

	const double *const basis_vals_i = basis_vals->data;

	Type x, y;
	double R_ij;

	for (k = 0; k < num_pts; k++){
		// Loop over all the (xi,eta) points

		x = 0;
		y = 0;

		for (i = 0; i < numI; i++){
			for (j = 0; j < numJ; j++){

				// Get the i,j NURBS basis function
				R_ij = basis_vals_i[i + j*numI + k*numI*numJ];

				x += R_ij * get_col_Multiarray_T(j, control_pts_x)[i];
				y += R_ij * get_col_Multiarray_T(j, control_pts_y)[i];
			}
		}

		// Store the mapped values
		get_col_Multiarray_T(0, xyz)[k] = x;
		get_col_Multiarray_T(1, xyz)[k] = y;
	}

	// Free allocated arrays
	destructor_Multiarray_T(control_pts_x);
	destructor_Multiarray_T(control_pts_y);
	destructor_Multiarray_d(control_pts_w);
	destructor_Multiarray_d(xi_vals);
	destructor_Multiarray_d(eta_vals);
	destructor_const_Multiarray_d(N_p_xi);
	destructor_const_Multiarray_d(N_q_eta);
	destructor_const_Multiarray_d(basis_vals);

	return (struct const_Multiarray_T*) xyz;

}


// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Read the required geometry data for the Joukowski parametric domain into the \ref Geo_Data container.
static void read_data_joukowski
	(struct Geo_Data*const geo_data ///< \ref Geo_Data.
	);

/// \brief Read the required geometry data for the Gaussian bump parametric domain into the \ref Geo_Data container.
static void read_data_gaussian_bump
	(struct Geo_Data*const geo_data ///< \ref Geo_Data.
	);

/// \brief Read the required geometry data for the fixed cube parametric domain into the \ref Geo_Data container.
static void read_data_fixed_cube
	(struct Geo_Data*const geo_data ///< \ref Geo_Data.
	);

/// \brief Read the required geometry data for the NURBS domain into the \ref Geo_Data container.
/// 
static void read_data_NURBS 
	(struct Geo_Data*const geo_data ///< \ref Geo_Data.
	);

/// \brief Read the required geometry data for the Multipatch domain into the \ref Geo_Data container.
/// 
static void read_data_Multipatch 
	(struct Geo_Data*const geo_data ///< \ref Geo_Data.
	);

static struct Geo_Data get_geo_data (const char*const geo_name)
{
	static bool need_input = true;
	static struct Geo_Data geo_data;
	if (need_input) {
		need_input = false;
		if (strcmp(geo_name,"joukowski") == 0)
			read_data_joukowski(&geo_data);
		else if (strcmp(geo_name,"gaussian_bump") == 0)
			read_data_gaussian_bump(&geo_data);
		else if (strcmp(geo_name,"fixed_cube") == 0)
			read_data_fixed_cube(&geo_data);
		else if (strcmp(geo_name, "NURBS") == 0)
			read_data_NURBS(&geo_data);
		else if (strcmp(geo_name, "Multipatch") == 0)
			read_data_Multipatch(&geo_data);
		else
			EXIT_ERROR("Unsupported: %s.\n",geo_name);
	}
	return geo_data;
}

static Real romberg
	(Real (*f)(const Real x, const struct Function_Data_GP*const f_data), const Real a, const Real b, const Real acc,
	 const struct Function_Data_GP*const f_data)
{
	const size_t max_steps = 20;
	Real R1[max_steps], R2[max_steps]; //buffers
	Real *Rp = &R1[0], *Rc = &R2[0]; //Rp is previous row, Rc is current row
	Real h = (b-a); //step size
	Rp[0] = (f(a,f_data) + f(b,f_data))*h*.5; //first trapezoidal step

	for (size_t i = 1; i < max_steps; ++i) {
		h /= 2.0;
		Real c = 0;
		const size_t ep = 1UL << (i-1); //2^(i-1)
		for (size_t j = 1; j <= ep; ++j)
			c += f(a+(double)(2*j-1)*h,f_data);
		Rc[0] = h*c + .5*Rp[0]; //R(i,0)

		for (size_t j = 1; j <= i; ++j) {
			Real n_k = pow(4.0,(double)j);
			Rc[j] = (n_k*Rc[j-1] - Rp[j-1])/(n_k-1); //compute R(i,j)
		}

		if (i > 1 && fabs(Rp[i-1]-Rc[i]) < acc)
			return Rc[i-1];

		//swap Rn and Rc as we only need the last row
		Real *rt = Rp;
		Rp = Rc;
		Rc = rt;
	}

	assert(fabs(Rp[max_steps-2]-Rp[max_steps-1]) < acc);
	return Rp[max_steps-1]; //return our best guess
}

static Real f_al_gaussian_bump (const Real x, const struct Function_Data_GP*const f_data)
{
	const Real df_dx = f_gaussian_bump(x,1,f_data);
	return sqrt(1.0 + df_dx*df_dx);
}

// Level 1 ********************************************************************************************************** //

static void read_data_joukowski (struct Geo_Data*const geo_data)
{
	const int count_to_find = 4;

	FILE* input_file = fopen_input('g',NULL,NULL); // closed

	int count_found = 0;
	char line[STRLEN_MAX];
	while (fgets(line,sizeof(line),input_file)) {
		read_skip_string_count_c_style_d("x_scale", &count_found,line,&geo_data->x_scale);
		read_skip_string_count_c_style_d("r_i",     &count_found,line,&geo_data->r_i);
		read_skip_string_count_c_style_d("r_o",     &count_found,line,&geo_data->r_o);
		read_skip_string_count_c_style_d("s_offset",&count_found,line,&geo_data->s_offset);
	}
	fclose(input_file);

	assert(count_found == count_to_find);
}

static void read_data_gaussian_bump (struct Geo_Data*const geo_data)
{
	const int count_to_find = 6;

	FILE* input_file = fopen_input('g',NULL,NULL); // closed

	int count_found = 0;
	char line[STRLEN_MAX];
	while (fgets(line,sizeof(line),input_file)) {
		read_skip_string_count_c_style_d("a",    &count_found,line,&geo_data->a);
		read_skip_string_count_c_style_d("b",    &count_found,line,&geo_data->b);
		read_skip_string_count_c_style_d("c",    &count_found,line,&geo_data->c);
		read_skip_string_count_c_style_d("d",    &count_found,line,&geo_data->d);
		read_skip_string_count_c_style_d("h",    &count_found,line,&geo_data->h);
		read_skip_string_count_c_style_d("x_max",&count_found,line,&geo_data->x_max);
	}
	fclose(input_file);

	assert(count_found == count_to_find);
}

static void read_data_NURBS (struct Geo_Data*const geo_data)
{

	// Read the NURBS data from the geometry_parameters.geo file

	// NOTE: Consider only 2D patches for now
	assert(DIM == 2);

	// The following parameters must be listed in the file.
	//	- P (order in the xi direction)
	// 	- Q (order in the eta direction)
	//	- knots_xi
	//  - knots_eta
	// 	- Control Points and Weights
	//  - Control Point Connectivity
	const int count_to_find = 6;

	// Get the file pointer to the geometry file
	FILE* input_file = fopen_input('g',NULL,NULL); // closed

	int count_found = 0;
	char line[STRLEN_MAX];

	int i, len_knots_xi, len_knots_eta;

	// Read in the information from the file
	while (fgets(line,sizeof(line),input_file)) {

		// Get the P and Q information
		read_skip_string_count_i("P(xi_order)", &count_found, line, &geo_data->P);
		read_skip_string_count_i("Q(eta_order)", &count_found, line, &geo_data->Q);

		// The knot vectors in the eta and xi directions:
		if (strstr(line, "knots_xi")){
			
			read_skip_string_count_i("knots_xi", &count_found, line, &len_knots_xi);	
			struct Multiarray_R* knots_xi = constructor_empty_Multiarray_R('C',2,(ptrdiff_t[]){len_knots_xi,1});
			
			for (i = 0; i < len_knots_xi; i++){
				fgets(line,sizeof(line),input_file);  // read the next line
				read_skip_d_1(line, 0, &get_col_Multiarray_R(0, knots_xi)[i], 1);
			}

			geo_data->knots_xi = knots_xi;
		}

		if (strstr(line, "knots_eta")){

			read_skip_string_count_i("knots_eta", &count_found, line, &len_knots_eta);
			struct Multiarray_R* knots_eta = constructor_empty_Multiarray_R('C',2,(ptrdiff_t[]){len_knots_eta,1});
			
			for (i = 0; i < len_knots_eta; i++){
				fgets(line,sizeof(line),input_file);  // read the next line
				read_skip_d_1(line, 0, &get_col_Multiarray_R(0, knots_eta)[i], 1);
			}
			
			geo_data->knots_eta = knots_eta;

		}

		// Read the Control Point Data
		if (strstr(line,"Control_Point_Data")){
			
			// Get the number of pts
			int num_control_pts;
			read_skip_string_count_i("Control_Point_Data", &count_found, line, &num_control_pts);

			// Create the multiarray structure to hold the control point data
			struct Multiarray_T* control_points  = constructor_empty_Multiarray_T('C',2,(ptrdiff_t[]){num_control_pts,DIM}); // saved
			struct Multiarray_R* control_weights = constructor_empty_Multiarray_R('C',2,(ptrdiff_t[]){num_control_pts,1}); // saved
			 
			Type *x = get_col_Multiarray_T(0, control_points),
				 *y = get_col_Multiarray_T(1, control_points);

			Real *w = get_col_Multiarray_R(0, control_weights);

			// Read the pt data line by line
			for (i = 0; i < num_control_pts; i++){
				fgets(line,sizeof(line),input_file);
				
				// Point data stored in .geo file with real values 
				Real pt[3] = {0.0, 0.0, 0.0};
				read_skip_d_1(line, 0, pt, DIM+1);

				x[i] = pt[0];
				y[i] = pt[1];
				w[i] = pt[2];
			}

			geo_data->control_points  = control_points;
			geo_data->control_weights = control_weights;
		}

		// Read the Control Point Connectivity Data
		if (strstr(line,"Control_Point_Connectivity")){

			// Get the number of xi and eta points
			count_found++;  // Connectivity information was found
			int extents[2] = {0,0};
			read_skip_const_i_1(line,1,extents,2);

			// Create the multiarray structure to hold the connectivity data
			struct Multiarray_i* control_pt_wt_connectivity = 
				constructor_empty_Multiarray_i('R',2,(ptrdiff_t[]){extents[0], extents[1]});

			// Read in the connectivity data line by line
			for (i = 0; i < extents[0]; i++){
				fgets(line,sizeof(line),input_file);

				read_skip_i_1(line,0, get_row_Multiarray_i(i, control_pt_wt_connectivity), extents[1]);
			}

			geo_data->control_pt_wt_connectivity = control_pt_wt_connectivity;

		}
	}
	fclose(input_file);

	assert(count_found == count_to_find);
}
 
static void read_data_Multipatch (struct Geo_Data*const geo_data)
{

	// Read the NURBS data from the geometry_parameters.geo file

	// NOTE: Consider only 2D patches for now
	assert(DIM == 2);

	// The following parameters must be listed in the file.
	//	- P (order in the xi direction)
	// 	- Q (order in the eta direction)
	//	- knots_xi
	//  - knots_eta
	// 	- Control Points and Weights
	//  - Control Point Connectivity
	const int count_to_find_global=1; //amount of variables to find, not specific to patch
	const int count_to_find_per_patch = 6; // amount of variables to find, specific to each patch


	// Get the file pointer to the geometry file
	FILE* input_file = fopen_input('g',NULL,NULL); // closed
	char line[STRLEN_MAX];
	char buffer[STRLEN_MAX];

	int global_found = 0; //compare to count_to_find_global at end, make sure we have everything
	int patches_to_find;
	int patches_found=0;
	

	int i, len_knots_xi, len_knots_eta;

	// Read in the information from the file
	while (fgets(line,sizeof(line),input_file)) {
		read_skip_string_count_i("Number_patches", &global_found, line, &patches_to_find);
		}
	fclose(input_file);
	geo_data->NURBS_Patch_Data=malloc(sizeof(struct NURBS_Patch_Geo_Data)*(unsigned long)patches_to_find);

	for(int patch_index=0;patch_index<=patches_to_find;patch_index++){
		// For each patch index, open the file and search until we find the corresponding section. 
		// Read that section until the next patch starts. 
		int count_found = 0;
		input_file = fopen_input('g',NULL,NULL);
		while (fgets(line,sizeof(line),input_file)) {
			if (strstr(line, "Patch_Index")){
				sprintf(buffer, "Patch_Index %d\n", patch_index);
				if(!strcmp(line, buffer)){
					while(fgets(line,sizeof(line),input_file) && !strstr(line, "Patch_Index")){
						read_skip_string_count_i("P(xi_order)", &count_found, line, &geo_data->NURBS_Patch_Data[patch_index].P);
						read_skip_string_count_i("Q(eta_order)", &count_found, line, &geo_data->NURBS_Patch_Data[patch_index].Q);
						// The knot vectors in the eta and xi directions:
						if (strstr(line, "knots_xi")){
							read_skip_string_count_i("knots_xi", &count_found, line, &len_knots_xi);	
							struct Multiarray_R* knots_xi = constructor_empty_Multiarray_R('C',2,(ptrdiff_t[]){len_knots_xi,1});

							for (i = 0; i < len_knots_xi; i++){
								fgets(line,sizeof(line),input_file);  // read the next line
								read_skip_d_1(line, 0, &get_col_Multiarray_R(0, knots_xi)[i], 1);
							}
							geo_data->NURBS_Patch_Data[patch_index].knots_xi = knots_xi;
						}

						if (strstr(line, "knots_eta")){

							read_skip_string_count_i("knots_eta", &count_found, line, &len_knots_eta);
							struct Multiarray_R* knots_eta = constructor_empty_Multiarray_R('C',2,(ptrdiff_t[]){len_knots_eta,1});
							
							for (i = 0; i < len_knots_eta; i++){
								fgets(line,sizeof(line),input_file);  // read the next line
								read_skip_d_1(line, 0, &get_col_Multiarray_R(0, knots_eta)[i], 1);
							}
							geo_data->NURBS_Patch_Data[patch_index].knots_eta = knots_eta;

						}
								// Read the Control Point Data
						if (strstr(line,"Control_Point_Data")){
							
							// Get the number of pts
							int num_control_pts;
							read_skip_string_count_i("Control_Point_Data", &count_found, line, &num_control_pts);

							// Create the multiarray structure to hold the control point data
							struct Multiarray_T* control_points  = constructor_empty_Multiarray_T('C',2,(ptrdiff_t[]){num_control_pts,DIM}); // saved
							struct Multiarray_R* control_weights = constructor_empty_Multiarray_R('C',2,(ptrdiff_t[]){num_control_pts,1}); // saved
							
							Type *x = get_col_Multiarray_T(0, control_points),
								*y = get_col_Multiarray_T(1, control_points);

							Real *w = get_col_Multiarray_R(0, control_weights);

							// Read the pt data line by line
							for (i = 0; i < num_control_pts; i++){
								fgets(line,sizeof(line),input_file);
								
								// Point data stored in .geo file with real values 
								Real pt[3] = {0.0, 0.0, 0.0};
								read_skip_d_1(line, 0, pt, DIM+1);

								x[i] = pt[0];
								y[i] = pt[1];
								w[i] = pt[2];
							}
							geo_data->NURBS_Patch_Data[patch_index].control_points = control_points;
							geo_data->NURBS_Patch_Data[patch_index].control_weights = control_weights;
						}

						// Read the Control Point Connectivity Data
						if (strstr(line,"Control_Point_Connectivity")){

							// Get the number of xi and eta points
							count_found++;  // Connectivity information was found
							int extents[2] = {0,0};
							read_skip_const_i_1(line,1,extents,2);

							// Create the multiarray structure to hold the connectivity data
							struct Multiarray_i* control_pt_wt_connectivity = 
								constructor_empty_Multiarray_i('R',2,(ptrdiff_t[]){extents[0], extents[1]});

							// Read in the connectivity data line by line
							for (i = 0; i < extents[0]; i++){
								fgets(line,sizeof(line),input_file);

								read_skip_i_1(line,0, get_row_Multiarray_i(i, control_pt_wt_connectivity), extents[1]);
							}
							geo_data->NURBS_Patch_Data[patch_index].control_pt_wt_connectivity = control_pt_wt_connectivity;

						}
					}
					assert(count_to_find_per_patch == count_found);
					patches_found++;
					break;
				}
			}
		}
		fclose(input_file);
	}
		
	assert(patches_found=count_to_find_per_patch);
	assert(count_to_find_global == global_found);
}

static void read_data_fixed_cube (struct Geo_Data*const geo_data)
{
	const int count_to_find = 2;

	FILE* input_file = fopen_input('g',NULL,NULL); // closed

	int count_found = 0;
	char line[STRLEN_MAX];
	while (fgets(line,sizeof(line),input_file)) {
		read_skip_string_count_c_style_d("xyz_l",&count_found,line,&geo_data->xyz_l);
		read_skip_string_count_c_style_d("xyz_r",&count_found,line,&geo_data->xyz_r);
	}
	fclose(input_file);

	assert(count_found == count_to_find);
}

#include "undef_templates_geometry_parametric.h"
#include "undef_templates_volume_solver.h"
#include "undef_templates_multiarray.h"
#include "undef_templates_math_functions.h"
