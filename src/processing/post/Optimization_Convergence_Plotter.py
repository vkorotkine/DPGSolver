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


def add_convergence_plot(file_list, file_convergence_data, x_value_key, y_value_key, semilogy=False):

	"""
	Add a plot to the current open figure. Create a plot of x_value_key vs 
	y_value_key.
	"""
	
	for case_data_dict in file_convergence_data:

		# Get the case information from file_list
		case_info_tuple = file_list[file_convergence_data.index(case_data_dict)]

		# Get the dictionary of values
		xVals = case_data_dict[x_value_key]
		yVals = case_data_dict[y_value_key]

		if semilogy:
			plt.semilogy(xVals, yVals, c=case_info_tuple[4], linestyle=case_info_tuple[5], 
				label=case_info_tuple[3], marker=case_info_tuple[6])
		else:
			plt.plot(xVals, yVals, c=case_info_tuple[4], linestyle=case_info_tuple[5], 
				label=case_info_tuple[3], marker=case_info_tuple[6])


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



