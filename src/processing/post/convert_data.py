import sys
import re
import numpy as np
import string

from numpy import linalg as LA
# Need to import matplot lib as plt for use of the plotting output

#np.set_printoptions(precision=3)
np.set_printoptions(linewidth=200)
np.set_printoptions(formatter={'float': lambda x: format(x, '6.2e')})

class output_class:
	"""Stores output related information."""

	def init_name(self,name):
		if (name == 'std'):
			self.n_entries = 1
			self.type   = 'std'
			self.format = 'table'
		else:
			print("Error: Unsupported (output_class)."); sys.exit()


class input_data:
	"""Stores input data related information."""

	def __init__(self):
		self.case  = list()
		self.n_vars = 0
		self.ml_max = 0
		self.p_max  = 0


def initialize_input(output,file_name_no_ext):
	"""Initializes input data based on the desired output."""

	data_i = [ input_data() for i in range(0,output.n_entries)]

	if (output.name == 'std'):
		data_i[0].case = file_name_no_ext
	else:
		print("Error. Unsupported: output.name=",output.name); EXIT

	return data_i

def skip_lines(N,f):
	for n in range(0,N):
		line = f.readline()

def assign_block(A,cases_run,f):
	"""Fill data block A based on Cases Run."""

	n_rows = len(cases_run)
	n_cols = len(cases_run[0])

	re_floats = '[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?'

	non_zero_rows = np.amax(cases_run,axis=1)

	for i in range(0,n_rows):
#		if (non_zero_rows[i] != 1):
#			continue

		line   = f.readline()
		vals_d = [float(s) for s in re.findall(re_floats, line)]
#		print("line: ",line)
#		print(vals_d)

		j2 = 0
		for j in range(0,n_cols):
			if (cases_run[i,j] == 1):
				A[i,j] = vals_d[j]
				j2 += 1


def read_data(data_i):
	"""Reads in data required by specific output."""

	re_int    = '\d+'
	re_ints   = '\\b\d+\\b'
	re_floats = '[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?'

	for n in range(0,len(data_i)):
		data = data_i[n]

		f_name = data.case + '.txt'

		var_names = list()

		with open(f_name) as f:
			for line in f:
				if ('n_vars' in line):
					n_vars = np.int_(re.search(re_int,line).group())
				if ('ml_max' in line):
					ml_max = np.int_(re.search(re_int,line).group())
				if ('p_max' in line):
					p_max  = np.int_(re.search(re_int,line).group())

					# Inialize np arrays once dimensions are known
					block_size2 = (ml_max+1,p_max+1)
					block_size3 = (ml_max+1,p_max+1,n_vars)

					cases_run   = np.int_(  np.zeros(block_size2))
					h           = np.float_(np.zeros(block_size2))
					l2_errors   = np.float_(np.zeros(block_size3))
					conv_orders = np.float_(np.zeros(block_size3))

				if ('Cases Run' in line):
					for i in range(0,ml_max+1):
						line = f.readline()
						cases_run[i,:] = np.int_([int(s) for s in re.findall(re_ints, line)])

				if ('Mesh Size' in line):
					assign_block(h,cases_run,f)

				if ('L2 Errors' in line):
					for k in range(0,n_vars):
						skip_lines(1,f)
						line = f.readline()
						var_names.append(line.replace('\n',''))

						assign_block(l2_errors[:,:,k],cases_run,f)

				if ('Convergence Orders' in line):
					for k in range(0,n_vars):
						skip_lines(2,f)
						assign_block(conv_orders[:,:,k],cases_run,f)


			data.n_vars = n_vars
			data.ml_max = ml_max
			data.p_max  = p_max

			data.var_names   = var_names
			data.cases_run   = cases_run
			data.h           = h
			data.l2_errors   = l2_errors
			data.conv_orders = conv_orders

def f_write(f,n_tabs,string):
	"""Write string to file with specified number of tabs and newline"""
	for n in range(0,n_tabs):
		f.write('\t')
	f.write(string + '\n')

def write_std_tables(f,data_i,output):
	"""Write output to std latex table format."""

	for n in range(0,len(data_i)):
		data = data_i[n]

		n_vars = data.n_vars
		p_max  = data.p_max
		ml_max = data.ml_max

		var_names   = data.var_names
		cases_run   = data.cases_run
		h          = data.h
		l2_errors   = data.l2_errors
		conv_orders = data.conv_orders

		f_write(f,0,'\\begin{table}[!ht]')
#		f_write(f,0,'\\begin{table}[!htbp]') # Include additional table placement options
		f_write(f,0,'\\begin{center}')
		string  = '\\caption{Errors and Convergence Orders - '
		string += "{\\color{red} " + output.MeshType + "}" + " Meshes "
		string += " {\\color{red} " + output.name + '}}'
		f_write(f,0,string)
		f_write(f,0,'\\resizebox{\\textwidth}{!}{')

		string = '\\begin{tabular}{| l | c | '
		for k in range(0,n_vars): string += 'c '
		string += '| '
		for k in range(0,n_vars): string += 'c '
		string += '| }'
		f_write(f,0,string)
		f_write(f,1,'\\hline')

		string  = '\\multicolumn{2}{|c|}{} & \\multicolumn{' + '{:d}'.format(n_vars) + '}{c|}{$L^2$ Error} & '
		string +=                           '\\multicolumn{' + '{:d}'.format(n_vars) + '}{c|}{Conv. Order} \\\\'
		f_write(f,1,string)
		f_write(f,1,'\\hline')

		string = 'Order ($k$) & Mesh Size ($h$) '
		for k in range(0,n_vars): string += '& ' + var_names[k].replace("L^2-","") + ' '
		for k in range(0,n_vars): string += '& ' + var_names[k].replace("L^2-","") + ' '
		string += '\\\\'
		f_write(f,1,string)
		f_write(f,1,'\\hline')

		non_zero_cols = np.amax(cases_run,axis=0)
		non_zero_rows = np.amax(cases_run,axis=1)
		for j in range(0,p_max+1):
			if (non_zero_cols[j] == 0):
				continue

			for i in range(0,ml_max+1):
				if (cases_run[i][j] == 0):
					continue
				if (non_zero_rows[i] == 0):
					continue

				n_tabs = 1
				string = ''
				if i == 0:
					string += str(j) + '\t'
					n_tabs = 0

				string += '& ' + '{:.2e}'.format(h[i,j]) + ' '
				for k in range(0,n_vars): string += '& ' + '{:.2e}'.format(l2_errors[i,j,k]) + ' '

				if i == 0:
					for k in range(0,n_vars): string += '& -    '
				else:
					for k in range(0,n_vars): string += '& ' + '{:2.2f}'.format(conv_orders[i,j,k]) + ' '
				string += '\\\\'
				f_write(f,n_tabs,string)

			f_write(f,1,'\\hline')

		f_write(f,0,'\\end{tabular}')
		f_write(f,0,'}')
		f_write(f,0,'\\end{center}')
		f_write(f,0,'\\end{table}')
		f_write(f,0,'')

def output_figure(data_i,output):
	"""Output figure in vector image format."""

	print(len(data_i))
	for n in range(0,len(data_i)):
		data = data_i[n]

		n_vars = data.n_vars
		p_max  = data.p_max
		ml_max = data.ml_max

		var_names   = data.var_names
		cases_run   = data.cases_run
		h          = data.h
		l2_errors   = data.l2_errors
		conv_orders = data.conv_orders

		non_zero_cols = np.amax(cases_run,axis=0)
		non_zero_rows = np.amax(cases_run,axis=1)

		ind_i = np.arange(0)
		for i in range(0,ml_max+1):
			if (non_zero_rows[i] == 0):
				ind_i = np.append(ind_i,[i])
		ind_i = np.delete(np.arange(ml_max+1),ind_i)

		for j in range(0,p_max+1):
			if (non_zero_cols[j] == 0):
				continue

			data_x = h[ind_i,j]
#			print(data_x)

			for k in range(0,n_vars):
				data_y = l2_errors[ind_i,j,k]

		# var_names go in legend
#				print(data_y)
#				plt.plot(data_x,data_y,'-o')

#			print(" ")


#	plt.xscale('log')
#	plt.yscale('log')
#	axes = plt.gca()
#	axes.set_xlim([1.05*np.amin(h),1.05*np.amax(h)])
#	plt.savefig('latex_figure.eps', format='eps', dpi=1200)
	sys.exit()


def output_data(output,data_i):

	if ('table' in output.format):
		f = open('latex_tables.txt', 'w')

	if ('std' in output.type):
		write_std_tables(f,data_i,output)
	else:
		print("Error: Unsupported (output.type)."); sys.exit()

	if ('table' in output.format):
		f.close()




def convert_data(output,file_name_no_ext):
	"""Convert data from c code input to appropriate output format"""

	# Initialize input data class
	data_i = initialize_input(output,file_name_no_ext)

	# Read in and store data
	read_data(data_i)

	# Output data (latex table or eps figure)
	output_data(output,data_i)



if __name__ == '__main__':
	"""
	Command line arguments:
	1. Name of the error data file without the file extension.
	2. Path to the directory containing error data files.
	"""

	if (len(sys.argv) != 3):
		print("\nIncorrect number of inputs. Should be:")
		print("\t1. full_path_to_error_file.")
		print("\t2. error_file_name.")
		print("\n\n")
		EXIT

	output_dir     = sys.argv[2]
	file_name_part = sys.argv[1]
	file_name = output_dir+'/'+file_name_part

	output = output_class()

	output.name = 'std'
	output.MeshType = "NONE"

	output.init_name(output.name)

	print("\nOutputting data in {} format.\n\n".format(output.format))
	convert_data(output,file_name)
