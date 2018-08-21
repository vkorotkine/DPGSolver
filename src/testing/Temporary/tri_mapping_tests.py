

import numpy
import matplotlib.pyplot as plt
import math


CONST_PHYS_TRIANGLE_POINTS = [(-1.0,-1.0), (-0.75,-1.0), (-1.0,-0.75)]

mapping_functions = [
	lambda xi, eta: (-0.5*xi - (math.sqrt(3)/6.)*eta + 1./3.),
	lambda xi, eta: (+0.5*xi - (math.sqrt(3)/6.)*eta + 1./3.),
	lambda xi, eta: (+(math.sqrt(3)/3.)*eta + 1./3.),
]

del_mapping_functions = [
	[lambda xi, eta: -0.5, lambda xi, eta: -math.sqrt(3.)/6.],
	[lambda xi, eta:  0.5, lambda xi, eta: -math.sqrt(3.)/6.],
	[lambda xi, eta:  0.0, lambda xi, eta:  math.sqrt(3.)/3.]
]

CONST_REF_TRIANGLE_POINTS = [(-1., -1./math.sqrt(3)), (0., 2./math.sqrt(3)), (1., -1./math.sqrt(3))]

"""
CONST_CUB_PTS_TEST = [
(-7.2527e-01, -4.1874e-01),
(7.2527e-01, -4.1874e-01),
(0, 8.3747e-01),
(0, -3.9011e-01),
(3.3785e-01, 1.9506e-01),
(-3.3785e-01, 1.9506e-01)
]
"""


CONST_CUB_PTS_TEST = [
(8.8730e-01, -3.8215e-01),
(5.0000e-01, 2.8868e-01),
(1.1270e-01, 9.5950e-01)
]


def map_points(xi_eta_pts):

	"""
	Map each point in the xi, eta domain onto the physical domain 
	and return the list of mapped points. 
	"""

	mapped_pts = []

	for xi_eta_pt in xi_eta_pts:

		xi, eta = xi_eta_pt[0], xi_eta_pt[1]

		mapped_pt = [0.0, 0.0]

		for i in range(3):
			mapped_pt[0] += CONST_PHYS_TRIANGLE_POINTS[i][0]*mapping_functions[i](xi,eta)
			mapped_pt[1] += CONST_PHYS_TRIANGLE_POINTS[i][1]*mapping_functions[i](xi,eta)

		mapped_pts.append(mapped_pt)

	return mapped_pts


def get_reference_element_triangle_domain_points():

	"""
	Get the points on the reference domain (xi, eta) in the reference triangle
	"""

	xi_eta_pts = []

	eta_vals = numpy.linspace(-1./math.sqrt(3), 2./math.sqrt(3), 15)
	n_pts_per_eta_level = 10

	for eta_val in eta_vals:
		h = 2./math.sqrt(3) - eta_val
		l = 1.0 * (h/math.sqrt(3))

		xi_vals_for_const_eta = numpy.linspace(-l, l, n_pts_per_eta_level)

		for xi_val in xi_vals_for_const_eta:
			xi_eta_pts.append((xi_val, eta_val))

	return xi_eta_pts


def get_mapping_metrics(xi_eta_pts):

	"""
	Get the metric terms from the mapping
	"""

	d_by_dxi = [] # Partial with respect to xi
	d_by_deta = [] # Partial with respect to eta

	for xi_eta_pt in xi_eta_pts:

		xi, eta = xi_eta_pt[0], xi_eta_pt[1]

		dpt_by_dxi = [0.0, 0.0]
		dpt_by_deta = [0.0, 0.0]

		for i in range(3):
			dpt_by_dxi[0] += CONST_PHYS_TRIANGLE_POINTS[i][0]*del_mapping_functions[i][0](xi,eta)
			dpt_by_dxi[1] += CONST_PHYS_TRIANGLE_POINTS[i][1]*del_mapping_functions[i][0](xi,eta)

			dpt_by_deta[0] += CONST_PHYS_TRIANGLE_POINTS[i][0]*del_mapping_functions[i][1](xi,eta)
			dpt_by_deta[1] += CONST_PHYS_TRIANGLE_POINTS[i][1]*del_mapping_functions[i][1](xi,eta)

		d_by_dxi.append(dpt_by_dxi)
		d_by_deta.append(dpt_by_deta)

	return d_by_dxi, d_by_deta


def main():

	xi_eta_pts = get_reference_element_triangle_domain_points()
	mapped_pts = map_points(xi_eta_pts)

	x_vals = [pt[0] for pt in mapped_pts]
	y_vals = [pt[1] for pt in mapped_pts]

	plt.scatter(x_vals, y_vals)

	plt.grid()
	plt.show(block=True)


def test_cubature():

	mapped_cub_pts = map_points(CONST_CUB_PTS_TEST)

	triangle_vertices = CONST_PHYS_TRIANGLE_POINTS[:]
	triangle_vertices.append(triangle_vertices[0][:])

	print "mapped_cub_pts: "
	for pt in mapped_cub_pts:
		print pt

	print " "

	d_by_dxi, d_by_deta = get_mapping_metrics(CONST_CUB_PTS_TEST)
	print "d_by_dxi: "
	for v in d_by_dxi:
		print v

	print " "
	
	print "d_by_deta: "
	for v in d_by_deta:
		print v

	
	# Plotting

	# Plot the triangle
	x_vals = [pt[0] for pt in triangle_vertices]
	y_vals = [pt[1] for pt in triangle_vertices]
	
	plt.plot(x_vals, y_vals)

	# Scatter the points
	x_vals = [pt[0] for pt in mapped_cub_pts]
	y_vals = [pt[1] for pt in mapped_cub_pts]

	plt.scatter(x_vals, y_vals)

	plt.grid()
	plt.show(block=True)


if __name__ == "__main__":
	#main()
	test_cubature()


