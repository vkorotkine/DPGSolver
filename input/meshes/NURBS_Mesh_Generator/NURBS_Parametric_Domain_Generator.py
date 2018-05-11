
"""
**************************************************************
			NURBS Parametric Domain Generator (2D)
**************************************************************

This python script will be used to create a single NURBS patch to be used as a 
parametric domain in the DPGSolver. The knot domain (xi, eta) must be in [-1,1] 
(this is because the parametric domain in the DPGSolver code is in this domain). 

The setup process in the code will involve reading the knot vectors (xi and 
eta) as well as the design control point locations. The DPG code will then 
read in a GMSH file consisting of a reference domain (a square in the domain [-1,1] 
in both dimensions). This reference square will be discretized into the elements (on
the parametric domain). Using the NURBS mapping, each element's vertices and face nodes
will be found on the physical domain using their corresponding values on the 
parametric domain.

Optimization:
	- The NURBS patch will be able to be used for optimization as well. To ease the process
		and avoid the use of RBF, the patch will be such that there will be no RBF 
		volume points (no interior NURBS control points in the domain). The elements will
		not cross since the NURBS control points will be spaced apart and the parametric 
		domain will ensure the elements do not attain negative volumes.

"""

import Basis
import user_defined_patch
import numpy
import matplotlib.pyplot as plt
import sys

# The type of patch to use. Perhaps take this as a command line
# argument
CONST_Patch_Type = "user_defined_patch"
CONST_EPS = 1E-9
CONST_Output_file_name = "ChannelWithBump.nurbs_patch"

def NURBS_patch(xi,eta,BasisFunctionsList, ControlPoints_and_Weights):

	"""
	Knowing the basis functions and locations of the points,
	compute the parametric NURBS patch expression

	:param xi: The xi value to evaluate the patch at
	:param eta: The eta value to evaluate the patch at
	:param BasisFunctionsList: The list (list of lists in 2D) of basis functions
	:param ControlPoints_and_Weights: The list (list of lists in 2D) of Points. 
		The i,j index of the Point structure points to a list of the form [x,y,weight], 
		which correspond to the physical location and weight of a given control point

	:return: Value of the patch on the physical domain at the given xi,eta point on the 
		parametric domain.
		Using this function in a lambda expression will ease the plotting process
	"""

	physical_val = [0,0]

	numI = len(BasisFunctionsList)
	numJ = len(BasisFunctionsList[0])

	for i in range(numI):
		for j in range(numJ):
			for k in range(2):
				physical_val[k] = physical_val[k] + BasisFunctionsList[i][j](xi,eta) * ControlPoints_and_Weights[i][j][k]

	return physical_val


def plot_patch(patch_parameters):

	"""
	Plot the given patch using matplotlib. Display the control points as well
	as the constant knot lines. For the xi and eta
	values on the edge (boundary), subtract or add a small epsilon value to ensure
	they are in the patch domain.

	:param patch_parameters: The dictionary with the patch_parameters
	"""

	xiVector = patch_parameters["xiVector"]
	etaVector = patch_parameters["etaVector"]
	ControlPoints_and_Weights = patch_parameters["ControlPoints_and_Weights"]
	P = patch_parameters["P"]
	Q = patch_parameters["Q"]

	# Get the NURBS basis functions associated with each control point
	NURBS_basis_functions = Basis.get_NURBS_basis_functions(ControlPoints_and_Weights, P, Q, xiVector, etaVector)

	# Get the patch parametric function
	NURBS_patch_function = lambda xi,eta,BasisFunctionsList=NURBS_basis_functions, \
			ControlPoints_and_Weights=ControlPoints_and_Weights: \
			NURBS_patch(xi,eta, BasisFunctionsList, ControlPoints_and_Weights)

	# The distinct knot values along both parametric coordinates (xi and eta)
	xi_distinct_values = []
	eta_distinct_values = []

	for xi_val in xiVector:
		if xi_val not in xi_distinct_values:
			xi_distinct_values.append(xi_val)

	for eta_val in etaVector:
		if eta_val not in eta_distinct_values:
			eta_distinct_values.append(eta_val)

	# Discretize the xi and eta domains so we can plot the lines
	num_linspace_pts = 50
	eta_linspace = numpy.linspace(etaVector[0]+CONST_EPS, etaVector[-1]-CONST_EPS, num_linspace_pts)
	xi_linspace = numpy.linspace(xiVector[0]+CONST_EPS, xiVector[-1]-CONST_EPS, num_linspace_pts)

	# Plot the constant xi values:
	for xi_val in xi_distinct_values:
		if xi_distinct_values.index(xi_val) == 0:
			xi_val += CONST_EPS
		elif xi_distinct_values.index(xi_val) == (len(xi_distinct_values)-1):
			xi_val -= CONST_EPS

		xPlotPts = []
		yPlotPts = []

		for eta_val in eta_linspace:
			xy_pt = NURBS_patch_function(xi_val, eta_val)
			xPlotPts.append(xy_pt[0])
			yPlotPts.append(xy_pt[1])

		plt.plot(xPlotPts, yPlotPts, c='b')

	# Plot the constant eta values
	for eta_val in eta_distinct_values:
		if eta_distinct_values.index(eta_val) == 0:
			eta_val += CONST_EPS
		elif eta_distinct_values.index(eta_val) == (len(eta_distinct_values)-1):
			eta_val -= CONST_EPS

		xPlotPts = []
		yPlotPts = []

		for xi_val in xi_linspace:
			xy_pt = NURBS_patch_function(xi_val, eta_val)
			xPlotPts.append(xy_pt[0])
			yPlotPts.append(xy_pt[1])

		plt.plot(xPlotPts, yPlotPts, c='b')

	# Scatter the control points
	xPlotPts = []
	yPlotPts = []
	for i in range(len(ControlPoints_and_Weights)):
		for j in range(len(ControlPoints_and_Weights[0])):
			xPlotPts.append(ControlPoints_and_Weights[i][j][0])
			yPlotPts.append(ControlPoints_and_Weights[i][j][1])

	plt.scatter(xPlotPts, yPlotPts, c='r')

	plt.grid()
	plt.show(block=True)


def output_file(patch_parameters):

	"""
	Output the file with the patch parameters. The DPG code will
	read this file and use the values to perform the parametric
	mapping

	:param patch_parameters: The parameters used to define the patch
	"""

	xiVector = patch_parameters["xiVector"]
	etaVector = patch_parameters["etaVector"]
	ControlPoints_and_Weights = patch_parameters["ControlPoints_and_Weights"]
	P = patch_parameters["P"]
	Q = patch_parameters["Q"]

	num_xi_pts = len(ControlPoints_and_Weights)
	num_eta_pts = len(ControlPoints_and_Weights[0])

	with open(CONST_Output_file_name, "w") as fp:

		fp.write("**********************************\n")
		fp.write("         NURBS Patch File \n")
		fp.write("**********************************\n")
		fp.write("\n")

		fp.write("P(xi_order) %d \n" % P)
		fp.write("Q(eta_order) %d \n" % Q)
		fp.write("\n")

		fp.write("num_xi_pts %d \n"%(num_xi_pts))
		fp.write("num_eta_pts %d \n"%(num_eta_pts))
		fp.write("\n")
		
		fp.write("Control Point X Values (xi = increasing rows, eta = increasing cols)\n")
		for i in range(num_xi_pts):
			for j in range(num_eta_pts):
				fp.write(" %.14e " % ControlPoints_and_Weights[i][j][0])
			fp.write("\n")

		fp.write("Control Point Y Values (xi = increasing rows, eta = increasing cols)\n")
		for i in range(num_xi_pts):
			for j in range(num_eta_pts):
				fp.write(" %.14e " % ControlPoints_and_Weights[i][j][1])
			fp.write("\n")

		fp.write("Control Point Weights Values (xi = increasing rows, eta = increasing cols)\n")
		for i in range(num_xi_pts):
			for j in range(num_eta_pts):
				fp.write(" %.14e " % ControlPoints_and_Weights[i][j][2])
			fp.write("\n")





def main():
	
	"""
	Main method of the patch generation module
		
	Command line argument options:
		plot = Plot the given patch to see what it looks like. This will
			plot the outer lines of the patch (-1, 1 constant xi and eta lines 
			for the 2D case)
		output_file = Output the parametric NURBS patch file which will be able
			to be processed by the DPG code.
		[No Arguments] = Run the plot and output_file case (together)

	:return : -

	"""

	# Get the patch parameters
	patch_parameters = user_defined_patch.get_patch_information()

	# Parse command line arguments
	if len(sys.argv) == 1:
		# No command line arguments so plot and output the file
		plot_patch(patch_parameters)
		output_file(patch_parameters)

	elif(sys.argv[1] == "plot"):
		# Only plot the patch
		plot_patch(patch_parameters)
	
	elif(sys.argv[1] == "output_file"):
		# Only output the patch file
		output_file(patch_parameters)


if __name__ == "__main__":
	main()




