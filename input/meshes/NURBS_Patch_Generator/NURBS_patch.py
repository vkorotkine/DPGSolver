"""Module: NURBS_patch.py

Contains methods for generating the NURBS patch
"""

def NURBS_patch(xi, eta, ctrl_pts_and_weights, basis_functions):
	
	"""Compute the value of the NURBS patch at a specified xi and eta point.
	
	Args:
	    xi (float): xi location.
	    eta (float): eta location.
	    ctrl_pts_and_weights (list): List of control points and their weights.
	    basis_functions (list): List of lambda expressions for the basis functions.
	
	Returns:
	    list: List of floats for the physical point on the patch at the specified (xi, eta) point.
	"""

	val = [0, 0]

	numI = len(ctrl_pts_and_weights)
	numJ = len(ctrl_pts_and_weights[0])

	for i in range(numI):
		for j in range(numJ):
			for k in range(2):
				val[k] += basis_functions[i][j](xi, eta)*ctrl_pts_and_weights[i][j][k]

	return val


