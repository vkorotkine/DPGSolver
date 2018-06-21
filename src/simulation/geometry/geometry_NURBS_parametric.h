#ifndef DPG__geometry_NURBS_parametric_h__INCLUDED
#define DPG__geometry_NURBS_parametric_h__INCLUDED


struct const_Multiarray_d;
struct const_Multiarray_c;
struct const_Multiarray_i;


const struct const_Multiarray_d *grad_xyz_NURBS_patch_mapping(
	const struct const_Multiarray_d* xi_eta_i, int P, int Q, 
	const struct const_Multiarray_d* knots_xi, 
	const struct const_Multiarray_d* knots_eta, 
	const struct const_Multiarray_d* control_points_and_weights,
	const struct const_Multiarray_i* control_point_connectivity);

const struct const_Multiarray_d *xyz_NURBS_patch_mapping(
	const struct const_Multiarray_d* xi_eta_i, int P, int Q, 
	const struct const_Multiarray_d* knots_xi, 
	const struct const_Multiarray_d* knots_eta, 
	const struct const_Multiarray_d* control_points,
	const struct const_Multiarray_d* control_weights,
	const struct const_Multiarray_i* control_pt_wt_connectivity);

const struct const_Multiarray_c *xyz_NURBS_patch_mapping_c(
	const struct const_Multiarray_d* xi_eta_i, int P, int Q, 
	const struct const_Multiarray_d* knots_xi, 
	const struct const_Multiarray_d* knots_eta, 
	const struct const_Multiarray_c* control_points,
	const struct const_Multiarray_d* control_weights,
	const struct const_Multiarray_i* control_pt_wt_connectivity);

void test_mapping();

#endif // DPG__geometry_NURBS_parametric_h__INCLUDED

