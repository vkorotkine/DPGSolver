"""
Module: Internal_Channel_Patch.py
------------------------------------------

Create a NURBS patch for a channel for internal flow. The patch will
be such that the middle three points will be lifted and the weights of the
two outer points (of the three middle points) will be 1/sqrt(2). This will
create the slight bump at the bottom of the channel. 

NOTE: Since the middle three points will be lifted, ensure that the bottom
	of the channel has an odd number of control points.

"""

import math
import numpy
import matplotlib.pyplot as plt
import Basis
import sys


# ==================================================
# 				Patch Parameters


# The properties of the patch in the xi direction (which
# traverses around the airfoil from the trailing edge, bottom
# surface to the leading edge and back)
CONST_P = 6
CONST_NUM_CONTROL_PTS_XI = 11


# Properties of the patch in the eta direction (eta increases in the 
# normal direction from the airfoil surface to the farfield)
CONST_Q = 1
CONST_ETA_KNOTS = [-1, -1, 1, 1]


# Domain Parameters
CONST_X_MIN = 0.0
CONST_X_MAX = 5.0

CONST_Y_MIN = 0.0
CONST_Y_MAX = 2.0

CONST_DELTA_Y = 0.5  # How much to lift middle points by

# ==================================================

def get_channel_top_parameters():

	"""
	Get the parameters for the NURBS Spline that will be used to 
	parametrize the top of the channel. This method will generate the 
	knot vector and the control point locations.

	NOTE: For now, the top of the channel has been set to be straight

	:return: Dictionary with the parameters for the B Spline
		return_dict = {
			'knots' : [list with the knots],
			'control_points': [[pt1_x, pt1_y], [pt2_x, pt2_y], ...],
			'P' : Order of the spline
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



	# Compute the control points and weights. Equally space the points along the
	# bottom channel surface
	control_points = []
	control_points_x_vals = numpy.linspace(CONST_X_MIN, CONST_X_MAX, n)
	control_points_y_vals = numpy.linspace(CONST_Y_MAX, CONST_Y_MAX, n)

	for i in range(n):
		control_points.append([control_points_x_vals[i], control_points_y_vals[i], 1.0])

	return {
		"knots" : knots,
		"control_points" : control_points,
		"P" : p
	}


def get_channel_bottom_parameters():

	"""
	Get the parameters for the NURBS Spline that will be used to 
	parametrize the bottom of the channel. This method will generate the 
	knot vector and the control point locations.

	:return: Dictionary with the parameters for the B Spline
		return_dict = {
			'knots' : [list with the knots],
			'control_points': [[pt1_x, pt1_y], [pt2_x, pt2_y], ...],
			'P' : Order of the spline
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



	# Compute the control points and weights. Equally space the points along the
	# bottom channel surface
	control_points = []
	control_points_x_vals = numpy.linspace(CONST_X_MIN, CONST_X_MAX, n)
	control_points_y_vals = numpy.linspace(CONST_Y_MIN, CONST_Y_MIN, n)

	for i in range(n):
		control_points.append([control_points_x_vals[i], control_points_y_vals[i], 1.0])
	
	# - Lift the middle three points in the y direction by CONST_DELTA_Y
	i_mid_index = int((n-1)/2.0)

	control_points[i_mid_index-1][1] += CONST_DELTA_Y
	#control_points[i_mid_index-1][2] = 1./math.sqrt(2.)

	control_points[i_mid_index][1] += CONST_DELTA_Y

	control_points[i_mid_index+1][1] += CONST_DELTA_Y
	#control_points[i_mid_index+1][2] = 1./math.sqrt(2.)

	return {
		"knots" : knots,
		"control_points" : control_points,
		"P" : p
	}


def get_patch_information():

	"""
	Return the Patch information

	:return : A dictionary holding key, value pairs
		for the patch information
		{
			xiVector: [],
			etaVector: [],
			ControlPoints_and_Weights: [],
			P : val,
			Q : val
		}
	"""

	bottom_parameters = get_channel_bottom_parameters()
	top_parameters = get_channel_top_parameters()

	xiVector = bottom_parameters["knots"]
	etaVector = CONST_ETA_KNOTS
	P = CONST_P
	Q = CONST_Q

	# Build the control points and weights structure
	num_xi_pts = len(bottom_parameters["control_points"])
	num_eta_pts = 2

	ControlPoints_and_Weights = []

	for i in range(num_xi_pts):
		col = []
		for j in range(num_eta_pts):
			col.append([None, None, None])
		ControlPoints_and_Weights.append(col)
	
	for i in range(num_xi_pts):
		ControlPoints_and_Weights[i][0] = bottom_parameters["control_points"][i]
		ControlPoints_and_Weights[i][1] = top_parameters["control_points"][i]
	
	patch_parameters = {
		"xiVector" : xiVector,
		"etaVector" : etaVector,
		"ControlPoints_and_Weights" : ControlPoints_and_Weights,
		"P" : P,
		"Q" : Q
	}

	return patch_parameters


def test():
	bottom_parameters = get_channel_bottom_parameters()


if __name__ == "__main__":
	test()



