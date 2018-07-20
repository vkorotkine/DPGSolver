"""
Module: Pressure_Distribution_Plotter.py

Plot the pressure distribution on the surface
"""

import matplotlib.pyplot as plt
import os
import os.path
import numpy as np
import math

# Absolute path to the directory with all the optimization results
CONST_OPTIMIZATION_DIR_ABS_PATH = "/Users/manmeetbhabra/Documents/McGill/Research/DPGSolver/input/input_files/euler/steady/NURBS"

CONST_File_list = [
	("ml4_P4_NURBS_Y_Pressure_Distribution.txt", "ml4, P4", "r", "-", True),
	("ml3_P2_NURBS_Y_Pressure_Distribution.txt", "ml3, P2", "b", "-", True)
]

CONST_Airfoil_Boolean = False

def read_pressure_file(file_abs_path):

	"""
	Read the pressure file and store the sorted list of points. If this is an airfoil
	case, then also reverse the ordering of the pressure data over the bottom surface

	:param file_abs_path: Path (absolute) to the file to read
	"""

	# Will hold tuples of the range of s values (arc length) and the
	# list of geometry control points
	pressure_data = []

	with open(file_abs_path, "r") as fp:

		num_pts = int(fp.readline().rstrip('\n'))

		for i in range(num_pts):
			line = fp.readline().rstrip('\n')

			pt = [float(x) for x in line.split()]
			pressure_data.append(pt)



	# Arrange the curve data now by sorting it with respect to the 
	# arc length position
	pressure_data = sorted(pressure_data, key=lambda x: x[0])

	# If this is an airfoil case, then close the curve and invert all points with
	# s > 0.5 (bottom surface)

	if CONST_Airfoil_Boolean:

		# Close the curve
		first_point = [pressure_data[0][0], pressure_data[0][1]]
		pressure_data.append(first_point)

		# Loop over the points and if a value has an s value bigger than 0.5 (bottom surface)
		# then subtract it from 1 and set that as the new position (so points wrap around)
		for pt in pressure_data:
			if pt[0] > 0.5:
				pt[0] = 1.0 - pt[0]

		# Mutliply all arc lengths by 2 so that the values are now percentage chord (x/c) 
		# as the x axis
		for pt in pressure_data:
			pt[0] *= 2.0

		# Turn the data in Cp data
		for pt in pressure_data:
			pt[1] = (pt[1] - CONST_P_inf)/CONST_q_inf

	return pressure_data


def main():
	
	"""
	The main function
	"""

	# Read each file and plot it

	fig, ax = plt.subplots() # create a new figure with a default 111 subplot

	for data_tuple in CONST_File_list:

		file = data_tuple[0]
		file_abs_path = os.path.join(CONST_OPTIMIZATION_DIR_ABS_PATH, file)
		press_data = read_pressure_file(file_abs_path)

		# plot the data
		xVals = []
		yVals = []

		for pt in press_data:
			xVals.append(pt[0])
			yVals.append(pt[1])

		ax.plot(xVals, yVals, c=data_tuple[2], linestyle=data_tuple[3], label=data_tuple[1])
		if data_tuple[4]:
			ax.scatter(xVals, yVals, c=data_tuple[2], s=10)


	# If this is an airfoil Cp plot, flip the y axis
	if CONST_Airfoil_Boolean:
		plt.gca().invert_yaxis()

	# plt.title('Normalized Pressure Distribution')
	ax.set_xlabel('s', fontsize=14)
	
	if CONST_Airfoil_Boolean:
		ax.set_ylabel(r'$C_{p}$', fontsize=14)
	else:
		ax.set_ylabel(r'$p/p_{T}$', fontsize=14)

	ax.legend(loc=3, fontsize=14)
	ax.grid()
	ax.legend()

	plt.show(block=True)
	return

	# ===========================================
	# Plot Zoom

	from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
	axins = zoomed_inset_axes(ax, 3.5, loc=4) # zoom-factor: 2.5, location: upper-left

	#axins.xaxis.set_visible('False')
	#axins.yaxis.set_visible('False')

	for data_tuple in CONST_File_list:

		file = data_tuple[0]
		file_abs_path = os.path.join(CONST_OPTIMIZATION_DIR_ABS_PATH, file)
		press_data = read_pressure_file(file_abs_path)

		# plot the data
		xVals = []
		yVals = []

		for pt in press_data:
			xVals.append(pt[0])
			yVals.append(pt[1])

		axins.plot(xVals, yVals, c=data_tuple[2], linestyle=data_tuple[3], label=data_tuple[1])
		if data_tuple[4]:
			axins.scatter(xVals, yVals, c=data_tuple[2], s=10)

	x1, x2, y1, y2 = 0.45, 0.55, 0.973, 0.976 # specify the limits
	axins.set_xlim(x1, x2) # apply the x-limits
	axins.set_ylim(y1, y2) # apply the y-limits

	plt.yticks(visible=False)
	plt.xticks(visible=False)

	from mpl_toolkits.axes_grid1.inset_locator import mark_inset
	mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")

	# ===========================================


	plt.show(block=True)


if __name__ == "__main__":
	main()
