
"""
**************************************************************
					Test File Generator
**************************************************************

This python script will be used to generate .data files to be used 
for validating the B spline and NURBS implementation in the code

"""

import Basis
import matplotlib.pyplot as plt
import numpy


def plot_basis_1D(basis_function_list, xi_min, xi_max):

	"""
	Plot the basis functions on the domain xi_min to 
	xi_max

	:param basis_function_list: The list of basis functions (lambda expressions)
	:param xi_min: The lower limit of the domain to plot
	:param xi_max: The upper limit of the domain to plot

	:return : -
	"""

	xi_min += 1E-14
	xi_max -= 1E-14

	num_plot_pts = 50
	plot_pts = numpy.linspace(xi_min, xi_max, num_plot_pts)

	for N_b in basis_function_list:

		y_vals = []

		for pt in plot_pts:
			y_vals.append(N_b(pt))

		plt.plot(plot_pts, y_vals)


	plt.show(block=True)

def write_test_results(fp, P, basis_functions, xi_vector, xi_vals, w_vector=None):

	"""
	Write the block of data for a given test

	:param fp: The file pointer
	:param basis_functions: The list of lambda expressions for the basis functions
	:param xi_vals: The values on the parametric domain to evaluate the 
		1D basis functions at
	:param w_vector: Optional parameter for NURBS cases

	:return : -
	"""

	# Print knot information
	fp.write("knots %d\n" % (len(xi_vector)))
	for knot_val in xi_vector:
		fp.write("%.14e " % (knot_val))
	fp.write("\n\n")

	# The order of the basis functions
	fp.write("P %d\n" % P)
	fp.write("\n")

	# Print xi values
	fp.write("xi_vals %d\n" % len(xi_vals))
	for val in xi_vals:
		fp.write("%.14e " % val)
	fp.write("\n\n")

	# Print the weights
	if w_vector is not None:
		fp.write("weights %d\n" % (len(w_vector)))
		for w in w_vector:
			fp.write("%.14e " % (w))
		fp.write("\n\n")

	# Print basis function evaluated at the given values
	fp.write("Basis_vals %d \n" % (len(basis_functions)))

	for val in xi_vals:
		# loop through all the xi values to evaluate the bases at

		for N_b in basis_functions:
			# loop through the basis functions
			fp.write("%.14e " % (N_b(val)))

		fp.write("\n")
	
	fp.write("\n\n")


def output_BSpline_test_file():

	"""
	Output the varying tests for the B Spline implementation.
	Output the results for a 1D test case (any higher dimension 
	will just be a tensor product)

	Test cases will have varying P values:
	
	- Case (1) P = 2
		xi_vector = [-2,-2,-2,-1,0,1,2,2,2]
	- Case (2) P = 3
		xi_vector = [-2,-2,-2,-1,0,1,2,2,2]
	- Case (3) P = 4
		xi_vector = [-2,-2,-2,-1,0,1,2,2,2]
	
	:param -
	
	:return : -

	"""

	# The file to write the data to
	fp = open("BSpline_bases.data", "w")

	# ==========================
	#          Case 1
	# ==========================

	fp.write("--------------------------------------------------------------------------------\n")
	fp.write("B Spline Basis 1D - Case 1\n")
	fp.write("\n")

	# - In this case, we have n = m - p - 1 = 9 - 2 - 1 = 6 basis functions
	xi_vector = [-2,-2,-2,-1,0,1,2,2,2]
	xi_vals = [-2, -0.5, 0.25, 2]  # Location to evaluate the Basis functions at
	P_case_1 = 2

	BSpline_basis = Basis.get_BSpline_basis_functions_1D(P_case_1, xi_vector)

	write_test_results(fp, P_case_1, BSpline_basis, xi_vector, xi_vals)

	# ==========================
	#          Case 2
	# ==========================

	fp.write("--------------------------------------------------------------------------------\n")
	fp.write("B Spline Basis 1D - Case 2\n")
	fp.write("\n")

	# - In this case, we have n = m - p - 1 = 9 - 3 - 1 = 5 basis functions
	xi_vector = [-2,-2,-2,-1,0,1,2,2,2]
	xi_vals = [-2, -0.5, 0.25, 2]  # Location to evaluate the Basis functions at
	P_case_2 = 3

	BSpline_basis = Basis.get_BSpline_basis_functions_1D(P_case_2, xi_vector)

	write_test_results(fp, P_case_2, BSpline_basis, xi_vector, xi_vals)

	# ==========================
	#          Case 3
	# ==========================

	fp.write("--------------------------------------------------------------------------------\n")
	fp.write("B Spline Basis 1D - Case 3\n")
	fp.write("\n")

	# - In this case, we have n = m - p - 1 = 9 - 3 - 1 = 5 basis functions
	xi_vector = [-2,-2,-2,-1,0,1,2,2,2]
	xi_vals = [-2, -0.5, 0.25, 2]  # Location to evaluate the Basis functions at
	P_case_3 = 4

	BSpline_basis = Basis.get_BSpline_basis_functions_1D(P_case_3, xi_vector)

	write_test_results(fp, P_case_3, BSpline_basis, xi_vector, xi_vals)	

	#plot_basis_1D(BSpline_basis, xi_vector[0], xi_vector[-1])
	#return

	fp.close()


def output_NURBS_test_file():

	"""
	Output the varying tests for the NURBS implementation.
	Output the results for a 1D test case (any higher dimension 
	will just be a tensor product)

	Test cases will have varying P values:
	
	- Case (1) P = 2
		xi_vector = [-2,-2,-2,-1,0,1,2,2,2]  // knot vector
		w_vector = [0.5, 1.75, 0.8, 1.8, 1.75, 1.5]  // weight vector
	- Case (2) P = 3
		xi_vector = [-2,-2,-2,-1,0,1,2,2,2]
		w_vector = [0.5, 1.75, 0.8, 1.8, 1.75] 
	- Case (3) P = 4
		xi_vector = [-2,-2,-2,-1,0,1,2,2,2]
		w_vector = [0.5, 1.75, 0.8, 1.8]
	
	:param -
	
	:return : -

	"""

	# The file to write the data to
	fp = open("NURBS_bases.data", "w")

	# ==========================
	#          Case 1
	# ==========================

	fp.write("--------------------------------------------------------------------------------\n")
	fp.write("NURBS Basis 1D - Case 1\n")
	fp.write("\n")

	# - In this case, we have n = m - p - 1 = 9 - 2 - 1 = 6 basis functions
	xi_vector = [-2,-2,-2,-1,0,1,2,2,2]
	w_vector = [0.5, 1.75, 0.8, 1.8, 1.75, 1.5]
	xi_vals = [-2, -0.5, 0.25, 2]  # Location to evaluate the Basis functions at
	P_case_1 = 2

	NURBS_basis = Basis.get_NURBS_basis_functions_1D(P_case_1, xi_vector, w_vector)

	write_test_results(fp, P_case_1, NURBS_basis, xi_vector, xi_vals, w_vector)

	# ==========================
	#          Case 2
	# ==========================

	fp.write("--------------------------------------------------------------------------------\n")
	fp.write("NURBS Basis 1D - Case 2\n")
	fp.write("\n")

	xi_vector = [-2,-2,-2,-1,0,1,2,2,2]
	w_vector = [0.5, 1.75, 0.8, 1.8, 1.75]
	xi_vals = [-2, -0.5, 0.25, 2]  # Location to evaluate the Basis functions at
	P_case_2 = 3

	NURBS_basis = Basis.get_NURBS_basis_functions_1D(P_case_2, xi_vector, w_vector)

	write_test_results(fp, P_case_2, NURBS_basis, xi_vector, xi_vals, w_vector)

	# ==========================
	#          Case 3
	# ==========================

	fp.write("--------------------------------------------------------------------------------\n")
	fp.write("NURBS Basis 1D - Case 3\n")
	fp.write("\n")

	xi_vector = [-2,-2,-2,-1,0,1,2,2,2]
	w_vector = [0.5, 1.75, 0.8, 1.8]
	xi_vals = [-2, -0.5, 0.25, 2]  # Location to evaluate the Basis functions at
	P_case_3 = 4

	NURBS_basis = Basis.get_NURBS_basis_functions_1D(P_case_3, xi_vector, w_vector)

	write_test_results(fp, P_case_3, NURBS_basis, xi_vector, xi_vals, w_vector)	

	fp.close()


def main():
	output_BSpline_test_file()
	output_NURBS_test_file()


if __name__ == "__main__":
	main()

