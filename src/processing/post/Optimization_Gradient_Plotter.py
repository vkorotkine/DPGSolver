"""
Module: Optimization_Gradient_Plotter.py

Plot the gradient of the optimization at the first design step.
"""

import matplotlib.pyplot as plt
import os
import os.path
import numpy as np
import math

# Absolute path to the directory with all the optimization results
CONST_OPTIMIZATION_DIR_ABS_PATH = "/Users/manmeetbhabra/Documents/McGill/Research/DPGSolver/build_2D/output/paraview/euler/steady/NURBS_Airfoil/Constrained_Target_CL_cases/NACA0012_Target_CL_0.25_EqualDOF"

# The list of files and the label to associate with them when plotting them. Each tuple
# contains the name of the file first and the label second. An empty label will result 
# in no legend. The tuples are in the form:
# (file name, legend name, color, linestyle, Scatter Boolean)
CONST_File_list = [
	
	# Target CL
	#("ml2_P1_NURBS_Y_CM_eq_0.1255_vfc_eq_1.0_Objective_Gradient.txt", "P = 1, ml = 2, NURBS Metrics", "r", "-", ".", False),
	#("ml2_P1_NURBS_NIso_CM_eq_0.1255_vfc_eq_1.0_Objective_Gradient.txt", "P = 1, ml = 2, Standard (Isoparametric)", "c", "-", ".", False),
	#("ml2_P1_NURBS_NSuper_CM_eq_0.1255_vfc_eq_1.0_Objective_Gradient.txt", "P = 1, ml = 2, Standard (Superparametric)", "m", "-", ".", False),

	# ("ml2_P2_NURBS_Y_CM_eq_0.1255_vfc_eq_1.0_Objective_Gradient.txt", "P = 2, ml = 2, NURBS Metrics", "g", "-", ".", False),
	# ("ml2_P2_NURBS_NIso_CM_eq_0.1255_vfc_eq_1.0_Objective_Gradient.txt", "P = 2, ml = 2, Standard (Isoparametric)", "c", "-", ".", False),
	# ("ml2_P2_NURBS_NSuper_CM_eq_0.1255_vfc_eq_1.0_Objective_Gradient.txt", "P = 2, ml = 2, Standard (Superparametric)", "m", "-", ".", False),

	#("ml2_P3_NURBS_Y_CM_eq_0.1255_vfc_eq_1.0_Objective_Gradient.txt", "P = 3, ml = 2, NURBS Metrics", "b", "-", ".", False),
	#("ml2_P3_NURBS_NIso_CM_eq_0.1255_vfc_eq_1.0_Objective_Gradient.txt", "P = 3, ml = 2, Standard (Isoparametric)", "c", "-", ".", False),
	#("ml2_P3_NURBS_NSuper_CM_eq_0.1255_vfc_eq_1.0_Objective_Gradient.txt", "P = 3, ml = 2, Standard (Superparametric)", "m", "-", ".", False),

	# ("ml1_P2_NURBS_Y_CM_eq_0.1255_vfc_eq_1.0_Objective_Gradient.txt", "P = 2, ml = 1, NURBS Metrics", "r", "-", ".", False),
	# ("ml1_P2_NURBS_NIso_CM_eq_0.1255_vfc_eq_1.0_Objective_Gradient.txt", "P = 2, ml = 1, Standard (Isoparametric)", "c", "-", ".", False),
	# ("ml1_P2_NURBS_NSuper_CM_eq_0.1255_vfc_eq_1.0_Objective_Gradient.txt", "P = 2, ml = 1, Standard (Superparametric)", "m", "-", ".", False),

	# Target CL Equal DOF
	("ml4_P1_NURBS_Y_CM_eq_0.1255_vfc_eq_1.0_Objective_Gradient.txt", "P = 1, ml = 4, NURBS Metrics", "r", "-", ".", False),
	("ml3_P2_NURBS_Y_CM_eq_0.1255_vfc_eq_1.0_Objective_Gradient.txt", "P = 2, ml = 3, NURBS Metrics", "c", "-", ".", False),
	("ml2_P3_NURBS_Y_CM_eq_0.1255_vfc_eq_1.0_Objective_Gradient.txt", "P = 3, ml = 2, NURBS Metrics", "m", "-", ".", False),


	# Inverse Design
	#("ml2_P1_NURBS_Y_Objective_Gradient.txt", "P = 1, ml = 2, NURBS Metrics", "r", "-", ".", False),
	#("ml2_P1_NURBS_NIso_Objective_Gradient.txt", "P = 1, ml = 2, Standard (Isoparametric)", "c", "-", ".", False),
	#("ml2_P1_NURBS_NSuper_Objective_Gradient.txt", "P = 1, ml = 2, Standard (Superparametric)", "m", "-", ".", False),

	#("ml2_P2_NURBS_Y_Objective_Gradient.txt", "P = 2, ml = 2, NURBS Metrics", "r", "-", ".", False),
	#("ml2_P2_NURBS_NIso_Objective_Gradient.txt", "P = 2, ml = 2, Standard (Isoparametric)", "c", "-", ".", False),
	#("ml2_P2_NURBS_NSuper_Objective_Gradient.txt", "P = 2, ml = 2, Standard (Superparametric)", "m", "-", ".", False),

	#("ml2_P3_NURBS_Y_Objective_Gradient.txt", "P = 3, ml = 2, NURBS Metrics", "r", "-", ".", False),
	#("ml2_P3_NURBS_NIso_Objective_Gradient.txt", "P = 3, ml = 2, Standard (Isoparametric)", "c", "-", ".", False),
	#("ml2_P3_NURBS_NSuper_Objective_Gradient.txt", "P = 3, ml = 2, Standard (Superparametric)", "m", "-", ".", False),
	
	
]


def read_gradient_file(file_abs_path):

	"""
	Read the gradient file and store the values in a list of points.

	:param file_abs_path: Path (absolute) to the file to read
	"""

	# Will hold the convergence data 
	gradient_data = []

	with open(file_abs_path, "r") as fp:

		while True:

			line = fp.readline().rstrip('\n')

			if line == "":
				break

			line_elems = line.split()
			gradient_data.append(float(line_elems[0]))

	return gradient_data


def plot_data(file_gradient_data):

	"""
	Plot the gradient data
	"""

	plt.figure(1)

	for case_gradient_data in file_gradient_data:

		# Get the case information from CONST_File_list
		case_info_tuple = CONST_File_list[file_gradient_data.index(case_gradient_data)]

		# Get the values
		yVals = case_gradient_data
		xVals = [i for i in range(len(yVals))]

		plt.plot(xVals, yVals, c=case_info_tuple[2], linestyle=case_info_tuple[3], 
			label=case_info_tuple[1], marker=case_info_tuple[4])


	x_label = "Design Point"
	y_label = "Gradient of Objective"

	plt.xlabel(x_label, fontsize=14)
	plt.ylabel(y_label, fontsize=14)
	plt.grid()
	plt.legend()

	plt.show(block=True)


def main():
	
	"""
	The main function
	"""

	# Read each file and plot it

	file_gradient_data = []

	for data_tuple in CONST_File_list:

		file = data_tuple[0]
		file_abs_path = os.path.join(CONST_OPTIMIZATION_DIR_ABS_PATH, file)
		case_gradient_data = read_gradient_file(file_abs_path)

		file_gradient_data.append(case_gradient_data)

	plot_data(file_gradient_data)



if __name__ == "__main__":
	main()



