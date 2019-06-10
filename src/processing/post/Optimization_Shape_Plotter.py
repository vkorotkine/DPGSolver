"""
Module: Optimization_Shape_Plotter.py

Plot the optimization initial and final configurations
"""

import os
import sys
import numpy
import math
import Basis
import matplotlib.pyplot as plt
import matplotlib.lines as lines

# The list of files and the label to associate with them when plotting them. Each tuple
# contains the name of the file first and the label second. An empty label will result 
# in no legend. The tuples are in the form:
# (file name, legend name, color, linestyle, Scatter Boolean)
# CONST_File_list = [
	
# 	# Inverse Design
# 	("geometry_parameters_NACA0012_tests_reference.geo", "Initial Profile", "k", "--", False), # Initial Profile
# 	#("geometry_parameters_NACA4412.geo", "Target Profile", "b", "--", False), # Target Profile

# 	("Optimized_NURBS_Patch.txt", "P = 1, ml = 2, NURBS Metrics", "r", "-", True),
# 	#("ml2_P1_NURBS_NIso_Optimized_NURBS_Patch.txt", "P = 1, ml = 2, Standard (Isoparameteric)", "c", "-", True),
# 	#("ml2_P1_NURBS_NSuper_Optimized_NURBS_Patch.txt", "P = 1, ml = 2, Standard (Superparameteric)", "m", "-", True),
# ]


def read_Patch_file(file_path):

	"""Read the patch file and load the NURBS information
	
	Args:
	    file_path (str): Path to file with the NURBS patch information
	
	Returns:
	    dict: Dictionary with the patch information
	"""

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


def plot_patch_points(case_patch_info, case_data_tuple, plot_points_xi_eta):

	"""
	Plot the patch values
	
	Args:
	    case_patch_info (dict): Dictionary with the information for the given patch
	    case_data_tuple (tuple): Tuple holding the information for the given patch (used to
	    	determine the line color, linestyle, etc.)
	    plot_points_xi_eta (list): List with the points (xi, eta) on the parametric 
	    	domain at which to plot the patch.
	
	Returns:
	    TYPE: -
	"""


	# Get the basis function lambda expressions
	R_pq = Basis.get_NURBS_basis_functions(case_patch_info["Control_Points_and_Weights"], 
		case_patch_info["P"], case_patch_info["Q"], case_patch_info["xiVector"], 
		case_patch_info["etaVector"])


	# Store the values of the basis functions at the required 
	# points on the knot domain. 
	numI, numJ = len(case_patch_info["Connectivity_Matrix"]), len(case_patch_info["Connectivity_Matrix"][0])
	
	basis_values = []
	
	for k in range(len(plot_points_xi_eta)):
		
		const_k_vals = []
		xi = plot_points_xi_eta[k][0]
		eta = plot_points_xi_eta[k][1]

		for i in range(numI):
			const_i_vals = []
			for j in range(numJ):
				const_i_vals.append(R_pq[i][j](xi, eta))
		
			const_k_vals.append(const_i_vals)

		basis_values.append(const_k_vals)


	x_vals, y_vals = [], []

	for pt_k in range(len(basis_values)):

		x, y = 0., 0.

		for i in range(numI):
			for j in range(numJ):

				x += basis_values[pt_k][i][j] * case_patch_info["Control_Points_and_Weights"][i][j][0]
				y += basis_values[pt_k][i][j] * case_patch_info["Control_Points_and_Weights"][i][j][1]

		x_vals.append(x)
		y_vals.append(y)

	curve_label = case_data_tuple[3]
	curve_color = case_data_tuple[4]
	line_style = case_data_tuple[5]
	scatter_bool = case_data_tuple[7]

	plt.plot(x_vals, y_vals, c=curve_color, label=curve_label, linestyle=line_style)

	if not scatter_bool:
		return

	x_scatter, y_scatter = [], []

	for i in range(numI):
		for j in range(numJ):
			x_scatter.append(case_patch_info["Control_Points_and_Weights"][i][j][0])
			y_scatter.append(case_patch_info["Control_Points_and_Weights"][i][j][1])

	plt.scatter(x_scatter, y_scatter, c=curve_color, s=15)


