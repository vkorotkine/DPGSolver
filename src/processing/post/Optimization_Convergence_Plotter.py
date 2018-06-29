"""
Module: Optimization_Convergence_Plotter.py

Plot the optimization convergence results by reading the convergence files
from the output directories
"""

import matplotlib.pyplot as plt
import os
import os.path
import numpy as np
import math

# Absolute path to the directory with all the optimization results
CONST_OPTIMIZATION_DIR_ABS_PATH = "/Users/manmeetbhabra/Documents/McGill/Research/DPGSolver/build_debug_2D/output/paraview/euler/steady/NURBS_Airfoil"
append_path = ""
CONST_OPTIMIZATION_DIR_ABS_PATH = os.path.join(CONST_OPTIMIZATION_DIR_ABS_PATH, append_path)

# The list of files and the label to associate with them when plotting them. Each tuple
# contains the name of the file first and the label second. An empty label will result 
# in no legend. The tuples are in the form:
# (file name, legend name, color, linestyle, Scatter Boolean)
CONST_File_list = [
	
	# Target CL
	("NACA0012_TargetCL0.24_P2_16x10_NURBSMetricY_BFGSmaxnorm1E-2/optimization_convergence.txt", "P = 2, 16x10, NURBS Metrics", "r", "-", ".", False),
	("NACA0012_TargetCL0.24_P2_16x10_NURBSMetricN_BFGSmaxnorm1E-2/optimization_convergence.txt", "P = 2, 16x10, Standard", "c", "-", ".", False),
	("NACA0012_TargetCL0.24_P2Superparametric_16x10_NURBSMetricN_BFGSmaxnorm1E-2/optimization_convergence.txt", "P = 2 (Sup), 16x10, Standard", "m", "-", ".", False),

	#("NACA0012_TargetCL0.24_P3_20x10_NURBSMetricY_BFGSmaxnorm1E-2/convergence.txt", "P = 2, 16x10, NURBS Metrics", "r", "-", ".", False),
	#("NACA0012_TargetCL0.24_P3_20x10_NURBSMetricN_BFGSmaxnorm1E-2/convergence.txt", "P = 2, 16x10, Standard", "c", "-", ".", False),
	#("NACA0012_TargetCL0.24_P3Superparametric_20x10_NURBSMetricN_BFGSmaxnorm1E-2/convergence.txt", "P = 2 (Sup), 16x10, Standard", "m", "-", ".", False),

]


def read_convergence_file(file_abs_path):

	"""
	Read the convergence file and store the list of points. Read the header
	first to know what each column in the file corresponds to. Then place the elements
	in a dictionary such that that the key is the header string and value is the list
	holding the numbers associated to that column

	:param file_abs_path: Path (absolute) to the file to read
	"""

	# Will hold the convergence data 
	conv_data_dictionary = {}

	with open(file_abs_path, "r") as fp:

		# The header strings (they are separated by spaces)
		header_terms = fp.readline().rstrip("\n").split()

		# Create the dictionary to hold all the values in a list
		for header_val in header_terms:
			conv_data_dictionary[header_val] = []


		while True:

			line = fp.readline().rstrip('\n')

			if line == "":
				break

			line_elems = line.split()

			pt = []

			iter_num = int(line_elems[0])  # First element is the iteration always
			pt.append(iter_num)

			# Rest of the values are floats
			for i in range(1, len(line_elems)):
				pt.append(float(line_elems[i]))

			for header_index in range(len(header_terms)):
				conv_data_dictionary[header_terms[header_index]].append(pt[header_index])

	return conv_data_dictionary


def add_plot(file_convergence_data, x_value_key, y_value_key, semilogy=False):

	"""
	Add a plot to the current open figure. Create a plot of x_value_key vs 
	y_value_key.
	"""
	
	for case_data_dict in file_convergence_data:

		# Get the case information from CONST_File_list
		case_info_tuple = CONST_File_list[file_convergence_data.index(case_data_dict)]

		# Get the dictionary of values
		xVals = case_data_dict[x_value_key]
		yVals = case_data_dict[y_value_key]

		if semilogy:
			plt.semilogy(xVals, yVals, c=case_info_tuple[2], linestyle=case_info_tuple[3], 
				label=case_info_tuple[1], marker=case_info_tuple[4])
		else:
			plt.plot(xVals, yVals, c=case_info_tuple[2], linestyle=case_info_tuple[3], 
				label=case_info_tuple[1], marker=case_info_tuple[4])


	x_label = ""
	for s in x_value_key.split("_"):
		x_label += s + " "

	y_label = ""
	for s in y_value_key.split("_"):
		y_label += s + " "		

	plt.xlabel(x_label, fontsize=14)
	plt.ylabel(y_label, fontsize=14)
	plt.grid()
	plt.legend()


def plot_data(file_convergence_data):

	"""
	Plot the convergence data on the separate plots
	"""


	# log(Objective) Convergence vs. Iteration
	plt.figure(1)
	add_plot(file_convergence_data, "Design_Iteration", "Cost_Function", True)
	
	# log(Gradient) Convergence vs. Iteration
	plt.figure(2)
	add_plot(file_convergence_data, "Design_Iteration", "L2_Norm_Gradient_Cost_Function", True)

	# log(Objective) Convergence vs. CPU time
	plt.figure(3)
	add_plot(file_convergence_data, "CPU_Time(s)", "Cost_Function", True)
	
	# log(Gradient) Convergence vs. CPU time
	plt.figure(4)
	add_plot(file_convergence_data, "CPU_Time(s)", "L2_Norm_Gradient_Cost_Function", True)


	# Check if target CL exists and if so plot it
	CL_data_exists = True
	for case_data_dict in file_convergence_data:
		if "CL" not in case_data_dict:
			CL_data_exists = False

	if CL_data_exists:
		# CL vs. Iteration
		plt.figure(5)
		add_plot(file_convergence_data, "Design_Iteration", "CL")

		# CL vs CPU time
		plt.figure(6)
		add_plot(file_convergence_data, "CPU_Time(s)", "CL")

	plt.show(block=True)


def main():
	
	"""
	The main function
	"""

	# Read each file and plot it

	file_convergence_data = []

	for data_tuple in CONST_File_list:

		file = data_tuple[0]
		file_abs_path = os.path.join(CONST_OPTIMIZATION_DIR_ABS_PATH, file)
		conv_data = read_convergence_file(file_abs_path)

		file_convergence_data.append(conv_data)

	plot_data(file_convergence_data)



if __name__ == "__main__":
	main()



