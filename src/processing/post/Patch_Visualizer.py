"""
Module: Patch_Visualizer.py
------------------------------------------

Used to visualize the NURBS patch, read from a .geo file. This script will be used
to see the progress of the optimization by plotting a given line in the patch (an 
isocurve). To keep things efficient, the function will compute the basis function
values the first time and then will simply use the updated control point locations
to find the updated profiles.

"""

import os
import sys
import numpy
import math
import Basis
import matplotlib.pyplot as plt
import matplotlib.lines as lines
import time


# Location of the optimization Patch progress file
CONST_OPTIMIZATION_FILE = "../../../build_2D/output/paraview/euler/steady/NURBS_Airfoil/Optimized_NURBS_Patch.txt"

# Location of the initial geometry file
CONST_INITIAL_PATCH_FILE = "../../../input/input_files/euler/steady/NURBS_Airfoil/geometry_parameters.geo"

# The points (on the knot domain) to plot
CONST_ISOCURVE_PTS = []
t_vals = numpy.linspace(-1., 1., 100)
for t in t_vals:
	CONST_ISOCURVE_PTS.append((t, -1.))


# Plot Parameters
CONST_PLOT_X_RANGE = [-0.7, 0.7]
CONST_PLOT_Y_RANGE = [-0.1, 0.1]


def read_Patch_file(file_path):
	# read the patch file and load the information

	patch_information = {}

	with open(file_path, "r") as fp:

		while True:

			line = fp.readline()
			if line == "":
				break;

			line = line.rstrip("\n")

			if "P(xi_order)" in line:
				P_value = int(line.split()[-1])
				patch_information["P"] = P_value

			if "Q(eta_order)" in line:
				Q_value = int(line.split()[-1])
				patch_information["Q"] = Q_value

			if "knots_xi" in line:
				num_knots = int(line.split()[-1])

				knot_vals = []
				for i in range(num_knots):
					knot_vals.append(float(fp.readline().rstrip("\n")))

				patch_information["xiVector"] = knot_vals

			if "knots_eta" in line:
				num_knots = int(line.split()[-1])

				knot_vals = []
				for i in range(num_knots):
					knot_vals.append(float(fp.readline().rstrip("\n")))

				patch_information["etaVector"] = knot_vals

			if "Control_Point_Data" in line:
				num_pts = int(line.split()[-1])

				pt_vals = []
				for i in range(num_pts):
					pt_vals.append( [float(x) for x in fp.readline().rstrip("\n").split()])

				patch_information["Control_Point_List"] = pt_vals

			if "Control_Point_Connectivity" in line:
				line = line.split()
				
				num_xi = int(line[-2])
				num_eta = int(line[-1])

				connectivity_matrix = []
				for i in range(num_xi):
					connectivity_matrix.append([int(x) for x in fp.readline().rstrip("\n").split()])

				patch_information["Connectivity_Matrix"] = connectivity_matrix


	# Create the control point net
	Control_Points_and_Weights = []

	numI = len(patch_information["Connectivity_Matrix"])
	numJ = len(patch_information["Connectivity_Matrix"][0])

	for i in range(numI):
		row_vals = []
		for j in range(numJ):
			pt_index = patch_information["Connectivity_Matrix"][i][j]
			row_vals.append(patch_information["Control_Point_List"][pt_index])
		Control_Points_and_Weights.append(row_vals)

	patch_information["Control_Points_and_Weights"] = Control_Points_and_Weights

	return patch_information





def plot_patch_points(basis_values, Control_Points_and_Weights, Control_Points_and_Weights_opt):

	"""
	Plot the patch values
	"""

	x_vals = []
	y_vals = []

	x_vals_opt = []
	y_vals_opt = []

	numI = len(Control_Points_and_Weights)
	numJ = len(Control_Points_and_Weights[0])

	for pt_k in range(len(basis_values)):

		x = 0.
		y = 0.

		x_opt = 0.
		y_opt = 0.

		for i in range(numI):
			for j in range(numJ):

				x += basis_values[pt_k][i][j] * Control_Points_and_Weights[i][j][0]
				y += basis_values[pt_k][i][j] * Control_Points_and_Weights[i][j][1]

				x_opt += basis_values[pt_k][i][j] * Control_Points_and_Weights_opt[i][j][0]
				y_opt += basis_values[pt_k][i][j] * Control_Points_and_Weights_opt[i][j][1]


		x_vals.append(x)
		y_vals.append(y)

		x_vals_opt.append(x_opt)
		y_vals_opt.append(y_opt)


	x_scatter = []
	y_scatter = []

	x_scatter_opt = []
	y_scatter_opt = []

	for i in range(numI):
		for j in range(numJ):
			x_scatter.append(Control_Points_and_Weights[i][j][0])
			y_scatter.append(Control_Points_and_Weights[i][j][1])

			x_scatter_opt.append(Control_Points_and_Weights_opt[i][j][0])
			y_scatter_opt.append(Control_Points_and_Weights_opt[i][j][1])
			

	plt.plot(x_vals, y_vals)
	plt.scatter(x_scatter, y_scatter)

	plt.plot(x_vals_opt, y_vals_opt)
	plt.scatter(x_scatter_opt, y_scatter_opt)


def check_progress(patch_information, basis_values):

	#plt.draw()
	plt.pause(0.1)

	while True:

		time.sleep(1.0)
		
		plt.cla()

		patch_information_opt = read_Patch_file(CONST_OPTIMIZATION_FILE)

		plot_patch_points(basis_values, patch_information["Control_Points_and_Weights"],
			patch_information_opt["Control_Points_and_Weights"])

		plt.grid()

		plt.gca().set_xlim(CONST_PLOT_X_RANGE)
		plt.gca().set_ylim(CONST_PLOT_Y_RANGE)


		#plt.draw()
		plt.pause(0.1)


def main():

	patch_information = read_Patch_file(CONST_INITIAL_PATCH_FILE)
	patch_information_opt = read_Patch_file(CONST_OPTIMIZATION_FILE)

	# Get the basis function lambda expressions
	R_pq = Basis.get_NURBS_basis_functions(patch_information["Control_Points_and_Weights"], 
		patch_information["P"], patch_information["Q"], patch_information["xiVector"], 
		patch_information["etaVector"])


	# Store the values of the basis functions at the required 
	# points on the knot domain. 
	numI = len(patch_information["Connectivity_Matrix"])
	numJ = len(patch_information["Connectivity_Matrix"][0])
	
	basis_values = []
	
	for k in range(len(CONST_ISOCURVE_PTS)):
		
		const_k_vals = []
		xi = CONST_ISOCURVE_PTS[k][0]
		eta = CONST_ISOCURVE_PTS[k][1]

		for i in range(numI):
			const_i_vals = []
			for j in range(numJ):
				const_i_vals.append(R_pq[i][j](xi, eta))
		
			const_k_vals.append(const_i_vals)

		basis_values.append(const_k_vals)

	print "computed basis_values"

	# plot the initial line
	plot_patch_points(basis_values, patch_information["Control_Points_and_Weights"],
		patch_information_opt["Control_Points_and_Weights"])

	plt.grid()

	plt.gca().set_xlim(CONST_PLOT_X_RANGE)
	plt.gca().set_ylim(CONST_PLOT_Y_RANGE)

	if len(sys.argv) > 1 and sys.argv[1] == "progress":
		check_progress(patch_information, basis_values)
	else:
		plt.show(block=True)


if __name__ == "__main__":
	main()
	

