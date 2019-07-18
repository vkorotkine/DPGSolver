"""
Module: Airfoil_Patch.py
------------------------------------------

Creates a NURBS patch around an airfoil. Use a least squares approach to determine the 
locations of the control points. Take as input the number of control points desired and 
spline order. 

Use a standard approach to place the slipwall on one xi/eta location, 
farfield on another, etc...

"""

import math
import numpy
import matplotlib.pyplot as plt
import Basis
import sys
import scipy.integrate
from Airfoil_Patch_Helper import *


# ==================================================
# 				Patch Parameters

# The properties of the patch in the xi direction (which
# traverses around the airfoil from the trailing edge, bottom
# surface to the leading edge and back)
CONST_P = 4
CONST_NUM_CONTROL_PTS_XI = 17
#Have to get these to match somehow
#To add patches, need to modify:
# -get_BSpline_parameters: append new spline as needed
# -CONST_SPLINE_PATCH_CONNECT: says which patch index goes to which patch
# NOTE: BOTH PATCHES AND SPLINE INDICES DEPEND ON THE ORDER THAT THEY ARE APPENDED IN IN THE FUNCTION. 
CONST_NUM_PATCHES=2
CONST_SPLINE_PATCH_CONNECT=[[0,1], [1,2]]

#CONST_NUM_PATCHES=1
#CONST_SPLINE_PATCH_CONNECT=[[1,2]]
# Properties of the patch in the eta direction (eta increases in the 
# normal direction from the airfoil surface to the farfield)
CONST_Q = 1
CONST_ETA_KNOTS = [-1, -1, 1, 1]


CONST_R_FARFIELD = 20.

CONST_CONTINUOUS_APPROXIMATION = True

# ==================================================




# Plot Parameters (for testing)
CONST_PlotXRange = [-1.0, 1.0]
CONST_PlotYRange = [-0.6, 0.6]




def get_BSpline_parameters():
	"""
	This function is added in order to be able to add more splines without 
	creating a separate function for each. Important for multiple NURBS patches. 
	Note, each still needs to be described with separate parametric equation.
	Get a return list of parameter dictionaries that will be used to parametrize the
	airfoil and farfield surfaces. This method will generate the 
	knot vector and the control point locations using a least square
	approach.

	:return: List of dictionaries. Each dictionary with the parameters for the B Spline
		return_dict = {
			'knots' : [list with the knots],
			'control_points': [[pt1_x, pt1_y], [pt2_x, pt2_y], ...],
			'P' : Order of the spline,
			'Spline_Function' : Lambda expression for the spline
		}
	"""
	# Spline parameters
	n = CONST_NUM_CONTROL_PTS_XI
	p = CONST_P
	m = n + p + 1  # number of elements in knot vector

	# Create a uniform open knot vector (use an open knot vector so first and last
	# knot are repeated p+1 times). Also, domain of the knot vector is [-1,1]
	knots = []
 
	for i in range(p+1):
		knots.append(-1.)

	num_non_end_knots = m - 2*(p+1)
	delta_knot = 2./(num_non_end_knots+1)

	for i in range(1, num_non_end_knots+1):
		knots.append(-1. + i*delta_knot)

	for i in range(p+1):
		knots.append(1.)


	# Get the B Spline basis functions
	BSpline_Basis = Basis.get_BSpline_basis_functions_1D(p, knots)
	BSpline_parameters=[]
	
	if not CONST_CONTINUOUS_APPROXIMATION:

		control_points = get_BSpline_control_points_discrete_least_square(BSpline_Basis, knots, n)
		# Get the lambda expression for the spline
		spline_function = lambda xi, BSpline_Basis=BSpline_Basis, control_points=control_points: \
			Basis.get_Spline_Function(xi, control_points, BSpline_Basis)

	else:
		control_points = get_BSpline_control_points_continuous_least_square(BSpline_Basis, knots, n, airfoil_parametric_equation)
		spline_function = lambda xi, BSpline_Basis=BSpline_Basis, control_points=control_points: \
			Basis.get_Spline_Function(xi, control_points, BSpline_Basis)

	BSpline_parameters.append({
		"knots" : knots,
		"control_points" : control_points,
		"P" : p,
		"Spline_Function" : spline_function
		})
	#Create and append for farfield:
	inner_farfield = lambda xi: (0.5*CONST_R_FARFIELD*math.cos(-1.*(xi+1)*math.pi), 
							0.5*CONST_R_FARFIELD*math.sin(-1.*(xi+1)*math.pi))

	control_points = get_BSpline_control_points_continuous_least_square(BSpline_Basis, knots, n, inner_farfield)
	spline_function = lambda xi, BSpline_Basis=BSpline_Basis, control_points=control_points: \
		Basis.get_Spline_Function(xi, control_points, BSpline_Basis)
	
	BSpline_parameters.append({
	"knots" : knots,
	"control_points" : control_points,
	"P" : p,
	"Spline_Function" : spline_function
	})

	#Create and append for farfield:
	farfield = lambda xi: (CONST_R_FARFIELD*math.cos(-1.*(xi+1)*math.pi), 
							CONST_R_FARFIELD*math.sin(-1.*(xi+1)*math.pi))

	control_points = get_BSpline_control_points_continuous_least_square(BSpline_Basis, knots, n, farfield)
	spline_function = lambda xi, BSpline_Basis=BSpline_Basis, control_points=control_points: \
		Basis.get_Spline_Function(xi, control_points, BSpline_Basis)

	BSpline_parameters.append({
	"knots" : knots,
	"control_points" : control_points,
	"P" : p,
	"Spline_Function" : spline_function
	})
	
	return BSpline_parameters




def plot_airfoil():
	# Test method for plotting airfoil
	
	s = numpy.linspace(-1, 1, 201)

	x_vals = []
	y_vals = []

	for s_val in s:
		pt = airfoil_parametric_equation(s_val)
		x_vals.append(pt[0])
		y_vals.append(pt[1])

	plt.plot(x_vals, y_vals, c='k', linestyle='--')


def plot_spline(spline_function, xi_range, control_points):

	xi_min = xi_range[0]
	xi_max = xi_range[-1]

	num_plot_pts = 101

	x_vals = []
	y_vals = []

	xi_vals = numpy.linspace(xi_min, xi_max, 201)

	for xi in xi_vals:
		pt = spline_function(xi)
		x_vals.append(pt[0])
		y_vals.append(pt[1])

	plt.plot(x_vals, y_vals, c='r')

	x_ctrl_pt = []
	y_ctrl_pt = []

	for pt in control_points:
		x_ctrl_pt.append(pt[0])
		y_ctrl_pt.append(pt[1])

	plt.scatter(x_ctrl_pt, y_ctrl_pt)


def get_optimization_pts(ControlPoints_and_Weights):

	"""
	Go through the Control Point net and set the points that will 
	act as the optimization control points. This method will combine
	the points into a list and also specify each point's degrees of 
	freedom. 

	:param ControlPoints_and_Weights: The net (matrix) of control points and
		weights.

	:return : A list with the control points for the optimization. The list
		will hold tuples of the form (cntrl_pt, dof_x_bool, dof_y_bool)
		where dof_x_bool = 0 if the x direction cannot be used as a degree of 
		freedom and is 1 if it can be (same for dof_y_bool)
	"""

	# For the Airfoil case, the trailing edge will be fixed and all other points
	# will be able to freely move in the y direction

	# The leading edge will also be fixed (although the spline does not interpolate
	# this location)

	optimization_control_pt_list = []

	num_airfoil_pts = len(ControlPoints_and_Weights) - 1  # Subtract one because first and last pt identical

	leading_edge_pt_index = int(len(ControlPoints_and_Weights)/2)

	for i in range(1,num_airfoil_pts):
		if i != leading_edge_pt_index:
			optimization_control_pt_list.append((ControlPoints_and_Weights[i][0], 0, 1))


	# Place the limits on each design variable
	# TODO: For now, only y can be adjusted so set those limits. In the future
	# implement the general case
	optimization_control_pt_limit_list = []
	
	y_max = CONST_R_FARFIELD
	y_min = -1.0*CONST_R_FARFIELD

	for i in range(len(optimization_control_pt_list)):

		pt = optimization_control_pt_list[i][0]

		# Each design point will remain in their half
		if pt[1] >= 0:
			optimization_control_pt_limit_list.append((pt, 0.0, y_max))
		else:
			optimization_control_pt_limit_list.append((pt, y_min, 0.0))

	return optimization_control_pt_list, optimization_control_pt_limit_list


def get_patch_information():

	"""
	Return the information for all patches in a list of dicts. 

	:return : patch_parameters: A list of dictionaries. patch_parameters[patch_index]
	holds following dict, corresponding to patch with index patch_index. 
		{
			xiVector: [],
			etaVector: [],
			ControlPoints_and_Weights: [],
			P : val,
			Q : val
		}
	"""
	patch_parameters=[]
	# Get the B Spline and Farfield B Splines
	BSpline_parameters=get_BSpline_parameters()

	for patch_index in range(CONST_NUM_PATCHES):
		#SPLINE 0 and 1 corresponding to patch
		SP0=CONST_SPLINE_PATCH_CONNECT[patch_index][0]
		SP1=CONST_SPLINE_PATCH_CONNECT[patch_index][1]
		Airfoil_parameters=BSpline_parameters[SP0]
		Farfield_parameters=BSpline_parameters[SP1]
		xiVector = BSpline_parameters[SP0]["knots"]
		etaVector = CONST_ETA_KNOTS
		P = CONST_P
		Q = CONST_Q

		# Build the control points and weights structure
		num_xi_pts = len(BSpline_parameters[SP0]["control_points"])
		num_eta_pts = 2

		ControlPoints_and_Weights = []

		for i in range(num_xi_pts):
			col = []
			for j in range(num_eta_pts):
				col.append([None, None, None])
			ControlPoints_and_Weights.append(col)

		for i in range(num_xi_pts):
			ControlPoints_and_Weights[i][0] = [Airfoil_parameters["control_points"][i][0],
					Airfoil_parameters["control_points"][i][1], 1.]
			ControlPoints_and_Weights[i][1] = [Farfield_parameters["control_points"][i][0],
					Farfield_parameters["control_points"][i][1], 1.]

		# Wrap the patch around by closing it. That is, make the last point on the control
		# mesh equal the first one (xi = 1 and xi = -1)
		ControlPoints_and_Weights[-1][0] = ControlPoints_and_Weights[0][0]
		ControlPoints_and_Weights[-1][1] = ControlPoints_and_Weights[0][1]
			# Get the list of optimization control points
		optimization_control_pt_list, optimization_control_pt_limit_list = get_optimization_pts(ControlPoints_and_Weights)

		patch_parameters.append( {
			"xiVector" : xiVector,
			"etaVector" : etaVector,
			"ControlPoints_and_Weights" : ControlPoints_and_Weights,
			"P" : P,
			"Q" : Q,
			"Optimization_ControlPoints_and_Weights" : optimization_control_pt_list,
			"Optimization_ControlPoints_Limits" : optimization_control_pt_limit_list,
			"area_ref" : 1.0,
			"cm_le_x"  : -0.5,
			"cm_le_y"  : 0.0
		})

	return patch_parameters


def test():
	#get_patch_information()
	BSpline_parameters=get_BSpline_parameters()
	
	plot_airfoil()
	for parameters in BSpline_parameters: 
		plot_spline(parameters["Spline_Function"],
			parameters["knots"],
			parameters["control_points"])

	plt.grid()

	#plt.gca().set_xlim(CONST_PlotXRange)
	#plt.gca().set_ylim(CONST_PlotYRange)

	plt.show(block=True)


if __name__ == "__main__":
	test()







