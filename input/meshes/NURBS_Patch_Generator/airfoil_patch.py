"""
Module: airfoil_patch.py
------------------------------------------

Creates a NURBS patch around an airfoil. Use a least squares approach to determine the 
weights and the locations of the control points. Take as input the number of control points
desired. 

Use a standard approach to place the slipwall on one xi/eta location, farfield on another, etc...

"""




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


