"""
Module: Optimization_Gradient_Plotter.py

Plot the gradient of the optimization at the first design step.
"""

import matplotlib.pyplot as plt
import os
import os.path
import numpy as np
import math


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


def plot_gradient_data(file_list, file_gradient_data):

	"""
	Plot the gradient data
	"""

	for case_gradient_data in file_gradient_data:

		# Get the case information from file_list
		case_info_tuple = file_list[file_gradient_data.index(case_gradient_data)]

		# Get the values
		yVals = case_gradient_data
		xVals = [i for i in range(len(yVals))]

		plt.plot(xVals, yVals, c=case_info_tuple[4], linestyle=case_info_tuple[5], 
			label=case_info_tuple[3], marker=case_info_tuple[6])


	x_label = "Design Point"
	y_label = "Gradient of Objective"

	plt.xlabel(x_label, fontsize=14)
	plt.ylabel(y_label, fontsize=14)
	plt.grid()
	plt.legend()



