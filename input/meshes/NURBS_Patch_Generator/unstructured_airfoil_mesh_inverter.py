"""Module: unstructured_airfoil_mesh_inverter.py

Module used to invert a NACA0012 airfoil mesh to get the mesh on the
parametric domain.
"""

import Basis
import NURBS_patch
import matplotlib.pyplot as plt
import numpy as np
import sys
import math


class Node(object):
	
	def __init__(self, pt_physical):
		
		"""Constructor for the node object
		
		Args:
		    pt (list): List of floats for the point
		"""
		
		self.pt_physical = pt_physical
		self.pt_parametric = None


class Element(object):

	def __init__(self, element_type, physical_tag, node_connectivity):

		self.element_type = element_type
		self.physical_tag = physical_tag
		self.node_connectivity = node_connectivity


	def __str__(self):

		return "%s %s %s" % (self.element_type, self.physical_tag, self.node_connectivity)



def read_mesh_file(mesh_file_abs_path):

	"""Read the mesh file and get the nodes and elements
	
	Args:
	    mesh_file_abs_path (str): Mesh file absolute path
	"""
	
	# Holds the list of nodes
	node_list = []

	# Dictionary with the list of element information
	# Dictionary keys are the different element types. Value
	# are lists of elements
	elements_dict = {}

	with open(mesh_file_abs_path, "r") as fp:

		while True:

			line = fp.readline()

			if line == "":
				break

			line = line.rstrip("\n")
			
			if line == "$Nodes":
				# Load the nodes information

				num_nodes = int(fp.readline().rstrip("\n"))

				for _ in range(num_nodes):
					l = fp.readline().rstrip("\n").split()
					
					pt_list = []
					for i in range(3):
						pt_list.append(float(l[i+1]))

					node_list.append(Node(pt_list))

			if line == "$Elements":
				# Load the element information

				num_elems = int(fp.readline().rstrip("\n"))

				for _ in range(num_elems):

					elem_ints = [int(x) for x in fp.readline().rstrip("\n").split()]

					elem_index = elem_ints[0]
					elem_type = elem_ints[1]
					elem_num_tags = elem_ints[2]

					tags = []
					for i in range(elem_num_tags):
						tags.append(elem_ints[3+i])

					node_connectivity = []
					for i in range(2+elem_num_tags+1, len(elem_ints)):
						node_connectivity.append(elem_ints[i])

					if elem_type in elements_dict:
						elements_dict[elem_type].append(Element(elem_type, tags[0], node_connectivity))

					else:
						elements_dict[elem_type] = [Element(elem_type, tags[0], node_connectivity)]


	return node_list, elements_dict


def read_patch_file(file_path):
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


def plot_mesh(node_list):

	x_vals = [n.pt[0] for n in node_list]
	y_vals = [n.pt[1] for n in node_list]

	plt.scatter(x_vals, y_vals, s=4)

	plt.show(block=True)


def get_Jac(xi_i, nurbs_patch_grad_xi, nurbs_patch_grad_eta):
		
	J = np.zeros((2,2))

	# Partial with respect to xi
	J[0,0] = nurbs_patch_grad_xi(xi_i[0,0], xi_i[1,0])[0]
	J[1,0] = nurbs_patch_grad_xi(xi_i[0,0], xi_i[1,0])[1]

	# Partial with respect to eta
	J[0,1] = nurbs_patch_grad_eta(xi_i[0,0], xi_i[1,0])[0]
	J[1,1] = nurbs_patch_grad_eta(xi_i[0,0], xi_i[1,0])[1]

	return J


def get_G(xi_i, x_p, nurbs_patch):

	G = np.zeros((2,1))

	G[0,0] = nurbs_patch(xi_i[0,0], xi_i[1,0])[0] - x_p[0,0]
	G[1,0] = nurbs_patch(xi_i[0,0], xi_i[1,0])[1] - x_p[1,0]

	return G


def invert_physical_point(guess_pt_grid, x_point, nurbs_patch, nurbs_patch_grad_xi, nurbs_patch_grad_eta):

	# Get the point to guess (should have the radially closest physical point)
	num_i_guess_pts = len(guess_pt_grid)
	num_j_guess_pts = len(guess_pt_grid[0])

	closest_pt = [guess_pt_grid[0][0], 
		(x_point[0] - guess_pt_grid[0][0]["physical"][0])**2 + (x_point[1] - guess_pt_grid[0][0]["physical"][1])**2]

	for i in range(num_i_guess_pts):
		for j in range(num_j_guess_pts):
			
			rad_distance = (x_point[0] - guess_pt_grid[i][j]["physical"][0])**2 + (x_point[1] - guess_pt_grid[i][j]["physical"][1])**2

			if rad_distance < closest_pt[1]:
				closest_pt[0], closest_pt[1] = guess_pt_grid[i][j], rad_distance

	# print "closest point : "
	# print closest_pt

	# Point on the physical domain to invert
	x_physical = np.zeros((2,1))
	x_physical[0,0], x_physical[1,0] = x_point[0], x_point[1]

	# The point to iteratively progress to the paramteric domain coordinate
	xi_i = np.zeros((2,1))
	xi_i[0,0], xi_i[1,0] = closest_pt[0]["parametric"][0], closest_pt[0]["parametric"][1]

	# print "xi_i : "
	# print xi_i

	# The Newton-Raphson Method:
	while True:

		J = get_Jac(xi_i, nurbs_patch_grad_xi, nurbs_patch_grad_eta)
		G = get_G(xi_i, x_physical, nurbs_patch)
			
		if math.isnan(G[0,0]) or math.isnan(G[1,0]):
			"Encountered Nan in inversion : pt_physical = %s" % (x_point)
			return None
			#raise ValueError("Encountered Nan in inversion : pt_physical = %s" % (x_point))

		delta_xi = np.linalg.solve(J, -G)

		if np.linalg.norm(delta_xi) < 1E-10:
			break

		xi_i = xi_i + delta_xi

		#print "xi_i : %s" % (xi_i) 

	return xi_i


def initialize_guess_point_grid(nurbs_patch):

	# Create a grid of points to be used to determine initial guesses
	num_guess_pts_xi, num_guess_pts_eta = 8, 8
	guess_domain_xi, guess_domain_eta = [-0.95, 0.95], [-0.95, 0.95]

	# Store dictionaries with the parametric domain and physical domain locations
	# for the given points (equally spaced)
	guess_pt_grid = []

	for i in range(num_guess_pts_xi):
		row = []
		for j in range(num_guess_pts_eta):
			row.append(None)
		guess_pt_grid.append(row)

	xi_guess_vals = np.linspace(guess_domain_xi[0], guess_domain_xi[1], num_guess_pts_xi)
	eta_guess_vals = np.linspace(guess_domain_eta[0], guess_domain_eta[1], num_guess_pts_eta)

	for i in range(num_guess_pts_xi):
		for j in range(num_guess_pts_eta):

			xi_val, eta_val = xi_guess_vals[i], eta_guess_vals[j]

			x_pts = nurbs_patch(xi_val, eta_val)

			guess_pt_grid[i][j] = {
				"parametric": (xi_val, eta_val),
				"physical": (x_pts[0], x_pts[1])
			}

	return guess_pt_grid


def main():

	node_list, elements_dict = read_mesh_file(CONST_INPUT_MSH_FILE)

	patch_information = read_patch_file(CONST_NURBS_PATCH_FILE)

	R_ij = Basis.get_NURBS_basis_functions(patch_information["Control_Points_and_Weights"], 
		patch_information["P"], patch_information["Q"], patch_information["xiVector"], 
		patch_information["etaVector"])

	R_ij_grad = Basis.get_grad_NURBS_basis_functions(patch_information["Control_Points_and_Weights"], 
		patch_information["P"], patch_information["Q"], patch_information["xiVector"], 
		patch_information["etaVector"])

	numI = len(patch_information["Control_Points_and_Weights"])
	numJ = len(patch_information["Control_Points_and_Weights"][0])
	
	R_ij_grad_xi, R_ij_grad_eta = [], []

	for i in range(numI):
		row1, row2 = [], []

		for j in range(numJ):
			row1.append(R_ij_grad[i][j][0])
			row2.append(R_ij_grad[i][j][1])

		R_ij_grad_xi.append(row1)
		R_ij_grad_eta.append(row2)


	nurbs_patch = lambda xi, eta: NURBS_patch.NURBS_patch(xi, eta, 
		patch_information["Control_Points_and_Weights"], R_ij)

	nurbs_patch_grad_xi = lambda xi, eta: NURBS_patch.NURBS_patch(xi, eta, 
		patch_information["Control_Points_and_Weights"], R_ij_grad_xi)

	nurbs_patch_grad_eta = lambda xi, eta: NURBS_patch.NURBS_patch(xi, eta, 
		patch_information["Control_Points_and_Weights"], R_ij_grad_eta)


	guess_pt_grid = initialize_guess_point_grid(nurbs_patch)

	# Test point to invert
	xi_point_true = [0.85, 0.85]
	x_point = nurbs_patch(xi_point_true[0], xi_point_true[1])

	# ===================================
	#        Test Patch Functions
	# ===================================

	# print "\nxi_point_true : %s" % (xi_point_true)
	# print "x_point : %s" % (x_point)

	# if False:
	# 	xi_vals = np.linspace(-1, 1, 100)
	# 	x_vals_1, y_vals_1 = [], []
	# 	x_vals_2, y_vals_2 = [], []

	# 	for i in range(100):

	# 		pt = nurbs_patch(xi_vals[i], -1)
	# 		x_vals_1.append(pt[0])
	# 		y_vals_1.append(pt[1])

	# 		pt = nurbs_patch(xi_vals[i], 1)
	# 		x_vals_2.append(pt[0])
	# 		y_vals_2.append(pt[1])

	# 	plt.plot(x_vals_1, y_vals_1)
	# 	plt.plot(x_vals_2, y_vals_2)

	# 	plt.grid()
	# 	plt.show(block=True)

	# print "\npartials analytical"
	# print "partial_xi : %s" % (nurbs_patch_grad_xi(xi_point_true[0], xi_point_true[1]))
	# print "partial_eta : %s" % (nurbs_patch_grad_eta(xi_point_true[0], xi_point_true[1]))

	# print "\npartials finite difference"
	# h = 1E-6
	# xvec_plus_xi_h = nurbs_patch(xi_point_true[0]+h, xi_point_true[1])
	# partial_xi_fd = [(1./h)*(xvec_plus_xi_h[0] - x_point[0]), (1./h)*(xvec_plus_xi_h[1] - x_point[1])]
	# print "partial_xi_fd : %s" % (partial_xi_fd)
	# xvec_plus_eta_h = nurbs_patch(xi_point_true[0], xi_point_true[1]+h)
	# partial_eta_fd = [(1./h)*(xvec_plus_eta_h[0] - x_point[0]), (1./h)*(xvec_plus_eta_h[1] - x_point[1])]
	# print "partial_eta_fd : %s" % (partial_eta_fd)


	# ===================================
	#        Invert Patch Point
	# ===================================

	xi_point = invert_physical_point(guess_pt_grid, x_point, nurbs_patch, nurbs_patch_grad_xi, nurbs_patch_grad_eta)

	print xi_point

	#plot_mesh(node_list)


if __name__ == "__main__":
	main()




