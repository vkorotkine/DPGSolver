
"""
Module: user_defined_patch.py
------------------------------------------

Sets a user defined patch. That is, the user can manually set all 
the patch information, such as the knot vectors, weights and control
point locations. 

This module is good for testing purposes or for checking simple patch
geometries. 

IMPORTANT: Ensure that the NURBS patch parameter relations are respected
	(ex: The number of control points is dependent on the knot vector, etc...).
	Checks will be done to ensure this is the case

"""

import math

# ==============================================
#				Patch Information
# ==============================================

"""
# BUMP Patch

# Knot Vectors (xi, eta)
xiVector = [-1,-1,-1,-0.5,0,0,0.5,1,1,1]
etaVector = [-1,-1,-1,1,1,1]

# Control Points and Weights
# Order is ControlPoints[i][j]. i increases along xi and j along eta direction
#	[x, y, weight]
ControlPoints_and_Weights = [
[[0,0,1],               [0,2,1]],
[[1,0,1],               [1,2,1]],
[[2,0.5,1./math.sqrt(2)], [2,2,1]],
[[3,0.5,1],               [3,2,1]],
[[4,0.5,1./math.sqrt(2)], [4,2,1]],
[[5,0,1],               [5,2,1]],
[[6,0,1],               [6,2,1]]
]

# The order of the NURBS surface being made
# - P = xi direction order
# - Q = eta direction order
# NOTE: 
# 	n = number of control points in given xi or eta direction
#	p = order in given xi or eta direction
# 	m = length of knot vector
# 		Then, n + p + 1 = m -> Relation that must hold true
P = len(xiVector) - 1 - len(ControlPoints_and_Weights)
Q = len(etaVector) - 1 - len(ControlPoints_and_Weights[0])
"""

"""
# Test P = 1 Patch
# Knot Vectors (xi, eta)
xiVector = [-1,-1,0,1,1]
etaVector = [-1,-1,1,1]

# Control Points and Weights
# Order is ControlPoints[i][j]. i increases along xi and j along eta direction
#	[x, y, weight]
ControlPoints_and_Weights = [
[[0,0,1],               [0,2,1]],
[[1,0.5,1],             [1,2,1]],
[[2,0,1],               [2,2,1]]
]

# The order of the NURBS surface being made
# - P = xi direction order
# - Q = eta direction order
# NOTE: 
# 	n = number of control points in given xi or eta direction
#	p = order in given xi or eta direction
# 	m = length of knot vector
# 		Then, n + p + 1 = m -> Relation that must hold true
P = len(xiVector) - 1 - len(ControlPoints_and_Weights)
Q = len(etaVector) - 1 - len(ControlPoints_and_Weights[0])
"""

# BUMP Patch

# Knot Vectors (xi, eta)
xiVector = [-1,-1,-1,-0.5,0,0,0.5,1,1,1]
etaVector = [-1,-1,-1,1,1,1]

# Control Points and Weights
# Order is ControlPoints[i][j]. i increases along xi and j along eta direction
#	[x, y, weight]
ControlPoints_and_Weights = [
[[0,0,1],               [0,2,1]],
[[1,0,1],               [1,2,1]],
[[2,0.2,1./math.sqrt(2)], [2,2,1]],
[[3,0.2,1],               [3,2,1]],
[[4,0.2,1./math.sqrt(2)], [4,2,1]],
[[5,0,1],               [5,2,1]],
[[6,0,1],               [6,2,1]]
]

# The order of the NURBS surface being made
# - P = xi direction order
# - Q = eta direction order
# NOTE: 
# 	n = number of control points in given xi or eta direction
#	p = order in given xi or eta direction
# 	m = length of knot vector
# 		Then, n + p + 1 = m -> Relation that must hold true
P = len(xiVector) - 1 - len(ControlPoints_and_Weights)
Q = len(etaVector) - 1 - len(ControlPoints_and_Weights[0])

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

	patch_parameters = {
		"xiVector" : xiVector,
		"etaVector" : etaVector,
		"ControlPoints_and_Weights" : ControlPoints_and_Weights,
		"P" : P,
		"Q" : Q
	}

	return patch_parameters



