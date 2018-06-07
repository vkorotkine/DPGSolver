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

#NOTE: Need the same NURBS patch parameters for all patches to plot them

# Airfoil Case
CONST_PLOT_X_RANGE = [-0.6, 0.6]
CONST_PLOT_Y_RANGE = [-0.1, 0.1]


# Absolute path to the directory with all the optimization results
CONST_OPTIMIZATION_DIR_ABS_PATH = "/Users/jm-034232/Documents/McGill/Research/DPGSolver/build_2D/output/paraview/euler/steady/NURBS_Airfoil"
append_path = ""
CONST_OPTIMIZATION_DIR_ABS_PATH = os.path.join(CONST_OPTIMIZATION_DIR_ABS_PATH, append_path)

# The list of files and the label to associate with them when plotting them. Each tuple
# contains the name of the file first and the label second. An empty label will result 
# in no legend. The tuples are in the form:
# (file name, legend name, color, linestyle, Scatter Boolean)
CONST_File_list = [
	
	# Target CL
	("NACA0012_TargetCL0.24_P2_16x10_NURBSMetricY_BFGSmaxnorm1E-2/geometry_parameters_initial.geo", "Initial", "k", "--", False),
	
	("NACA0012_TargetCL0.24_P3_20x10_NURBSMetricY_BFGSmaxnorm1E-2/Optimized_NURBS_Patch.txt", "P = 2, 16x10, NURBS Metrics", "r", "-", True),
	("NACA0012_TargetCL0.24_P3_20x10_NURBSMetricN_BFGSmaxnorm1E-2/Optimized_NURBS_Patch.txt", "P = 2, 16x10, Standard", "c", "-", True),
	("NACA0012_TargetCL0.24_P3Superparametric_20x10_NURBSMetricN_BFGSmaxnorm1E-2/Optimized_NURBS_Patch.txt", "P = 2 (Sup), 16x10, Standard", "m", "-", True)

]

# The points (on the knot domain) to plot
CONST_ISOCURVE_PTS = []
t_vals = numpy.linspace(-1., 1., 100)
for t in t_vals:
	CONST_ISOCURVE_PTS.append((t, -1.))



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


def plot_patch_points(basis_values, Control_Points_and_Weights, case_data_tuple):

	"""
	Plot the patch values
	"""

	x_vals = []
	y_vals = []

	numI = len(Control_Points_and_Weights)
	numJ = len(Control_Points_and_Weights[0])

	for pt_k in range(len(basis_values)):

		x = 0.
		y = 0.

		for i in range(numI):
			for j in range(numJ):

				x += basis_values[pt_k][i][j] * Control_Points_and_Weights[i][j][0]
				y += basis_values[pt_k][i][j] * Control_Points_and_Weights[i][j][1]

		x_vals.append(x)
		y_vals.append(y)

	curve_label = case_data_tuple[1]
	curve_color = case_data_tuple[2]
	line_style = case_data_tuple[3]
	scatter_bool = case_data_tuple[4]

	plt.plot(x_vals, y_vals, c=curve_color, label=curve_label, linestyle=line_style)

	if not scatter_bool:
		return

	x_scatter = []
	y_scatter = []

	for i in range(numI):
		for j in range(numJ):
			x_scatter.append(Control_Points_and_Weights[i][j][0])
			y_scatter.append(Control_Points_and_Weights[i][j][1])

	plt.scatter(x_scatter, y_scatter, c=curve_color, s=15)


def main():
		
	# Read the the patch information

	cases_patch_info_list = []

	for data_tuple in CONST_File_list:

		file = data_tuple[0]
		file_abs_path = os.path.join(CONST_OPTIMIZATION_DIR_ABS_PATH, file)
		patch_info = read_Patch_file(file_abs_path)

		cases_patch_info_list.append(patch_info)


	# Use the first case to find the basis function values.
	# These will be used for plotting all the other curves (since the basis
	# functions are identical, only the control point positions change)

	# Get the basis function lambda expressions
	R_pq = Basis.get_NURBS_basis_functions(cases_patch_info_list[0]["Control_Points_and_Weights"], 
		cases_patch_info_list[0]["P"], cases_patch_info_list[0]["Q"], cases_patch_info_list[0]["xiVector"], 
		cases_patch_info_list[0]["etaVector"])


	# Store the values of the basis functions at the required 
	# points on the knot domain. 
	numI = len(cases_patch_info_list[0]["Connectivity_Matrix"])
	numJ = len(cases_patch_info_list[0]["Connectivity_Matrix"][0])
	
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


	# Plot the patch isocurves
	for case_index in range(len(cases_patch_info_list)):

		# Get the patch information
		patch_info = cases_patch_info_list[case_index]

		# Get the case plotting information
		case_data_tuple = CONST_File_list[case_index]

		plot_patch_points(basis_values, patch_info["Control_Points_and_Weights"], case_data_tuple)


	plt.legend()
	plt.grid()

	plt.gca().set_xlim(CONST_PLOT_X_RANGE)
	plt.gca().set_ylim(CONST_PLOT_Y_RANGE)

	plt.show(block=True)


if __name__ == "__main__":
	main()

