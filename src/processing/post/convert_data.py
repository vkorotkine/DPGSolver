#!/usr/bin/env python3

import sys
import re
import numpy as np
import string

from numpy import linalg as LA
# Need to import matplot lib as plt for use of the plotting output

from file_processing import skip_lines, f_write

#np.set_printoptions(precision=3)
np.set_printoptions(linewidth=200)
np.set_printoptions(formatter={'float': lambda x: format(x, '6.2e')})

class Output_Info:
	"""Stores output related information."""

	def __init__ (self,input_path,input_type,dim):
		self.root_path = input_path+"/" ###< Full path to root directory for the input files.
		self.type      = input_type

		assert int(dim) > 0 and int(dim) <= 3,"Invalid dimension: "+dim
		self.dim = int(dim)

		self.format = ""
		self.sub_index = [-1] ###< Index of the subset of the data to output if relevant.
		self.set_type_specific_data()
		self.set_var_names()

	def set_type_specific_data (self):
		if (self.type in ["std","ivs_adv","ivs_euler_sv","ivs_euler_sv_t","ivs_euler_gb",
		                  "adv_peterson","adv_man_1d"]):
			self.format = "table"
		else:
			assert 0,"Unsupported: "+str(self.type)

		if (self.type in ["ivs_adv","ivs_euler_sv","ivs_euler_gb","adv_peterson","adv_man_1d"]):
			self.sub_index = [0]
		elif (self.type in ["ivs_euler_sv_t"]):
			self.sub_index = [0,self.dim+2]

	def set_var_names (self):
		if (self.type == "std"):
			pass
		elif (self.type == "ivs_adv"):
			self.var_names = ["$L^2 (u)^{\\text{i}}$",
			                  "$L^2 (u)^{\\text{s}}_{p_c = 2(p+1)}$",
			                  "$L^2 (u)^{\\text{i}}_{\\bm{\\hat{n}_{ex}}}$",
			                  "$L^2 (u)^{\\text{s}}$",]
		elif (self.type == "ivs_euler_sv"):
			self.var_names = ["$L^2 (\\rho)^{\\text{i}}$",
			                  "$L^2 (\\rho)^{\\text{i}}_{\\bm{\\hat{n}_{ex}}}$",
			                  "$L^2 (\\rho)^{\\text{s}}$",]
		elif (self.type == "ivs_euler_sv_t"):
			self.var_names = ["$L^2 (\\rho)^{\\text{i}}$",
			                  "$L^2 (\\rho)^{\\text{i}}_{\\Gamma}$",
			                  "$L^2 (\\rho)^{\\text{s}}$",
			                  "$L^2 (s)^{\\text{i}}$",
			                  "$L^2 (s)^{\\text{i}}_{\\Gamma}$",
			                  "$L^2 (s)^{\\text{s}}$",]
		elif (self.type == "ivs_euler_gb"):
			self.var_names = ["$L^2 (s)^{\\text{i}}$",
			                  "$L^2 (s)^{\\text{s}}$",
			                  "$L^2 (s)^{\\text{s}}_{\\bm{\\hat{n}_{ex}}}$",]
		elif (self.type in ["adv_peterson","adv_man_1d"]):
			self.var_names = ["$L^2 (u)_{L^2}$",
			                  "$L^2 (u)_{\\text{OPG}}$",
			                  "$L^2 (u)_{\\text{DPG} - H_b^1}$",
			                  "$L^2 (u)_{\\text{DPG} - H_b^-}$",
			                  "$L^2 (u)_{\\text{DG}}$",]

class Input_Info:
	"""Stores input data related information."""

	def __init__ (self,input_path,input_type):
		self.root_path = input_path+"/" ###< Full path to root directory for the input files.
		self.type      = input_type     ###< Type of input files expected.

		self.rel_paths = [] ###< Array of relative paths from \ref Input_Info::root_path to the error files.
		self.set_type_specific_data()

		self.data = [] ###< Array of input data (1 entry per input file)

	def set_type_specific_data (self):
		if (self.type == "std"):
			self.rel_paths = ["",]
		elif (self.type == "ivs_adv"):
			self.rel_paths = ["advection_steady_vortex_dg_2d_ar5_iso/",
			                  "advection_steady_vortex_dg_2d_ar5_super_p_cub_p2/",
			                  "advection_steady_vortex_dg_2d_ar5_iso_exact_normals/",
			                  "advection_steady_vortex_dg_2d_ar5_super/",]
		elif (self.type == "ivs_euler_sv"):
			self.rel_paths = ["euler_supersonic_vortex_dg_2d_ar5_iso/",
			                  "euler_supersonic_vortex_dg_2d_ar5_iso_exact_normals/",
			                  "euler_supersonic_vortex_dg_2d_ar5_super/",]
		elif (self.type == "ivs_euler_sv_t"):
			self.rel_paths = ["euler_supersonic_vortex_dg_2d_ar1_iso_transonic_o/",
			                  "euler_supersonic_vortex_dg_2d_ar1_iso_transonic_o_boundary_only/",
			                  "euler_supersonic_vortex_dg_2d_ar1_super_transonic_o/",]
		elif (self.type == "ivs_euler_gb"):
			self.rel_paths = ["euler_gaussian_bump_dg_2d_ar5_iso/",
			                  "euler_gaussian_bump_dg_2d_ar5_super/",
			                  "euler_gaussian_bump_dg_2d_ar5_super_exact_normals/",]
		elif (self.type == "adv_peterson"):
			self.rel_paths = ["peterson_l2/",
			                  "peterson_opgc0/",
			                  "peterson_dpg_h1_upwind/",
			                  "peterson_dpg_adjoint/",
			                  "peterson_dg/",]
		elif (self.type == "adv_man_1d"):
			self.rel_paths = ["manufactured_1d_l2/",
			                  "manufactured_1d_opgc0/",
			                  "manufactured_1d_dpg_h1_upwind/",
			                  "manufactured_1d_dpg_adjoint/",
			                  "manufactured_1d_dg/",]
		else:
			assert 0,"Unsupported: "+str(self.type)

	def return_data_path (self,ind):
		""" Return the full path to the data file corresponding to the input "ind"ex. """
		return self.root_path+self.rel_paths[ind]

class Input_Data:
	""" Stores data from the input files. """

	def __init__ (self,input_path):
		self.input_path = input_path ###< Full path to where the input file is located.
		self.n_vars = 0
		self.ml_max = 0
		self.p_max  = 0

	def read_data (self):
		"""Reads in data required for the given input_type."""

		re_int    = '\d+'
		re_ints   = '\\b\d+\\b'
		re_floats = '[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?'

		f_name = self.input_path + "l2errs+convergence_standard.txt"

		var_names = list()

		with open(f_name) as f:
			for line in f:
				if ('n_vars' in line):
					n_vars = np.int_(re.search(re_int,line).group())
				if ('ml_max' in line):
					ml_max = np.int_(re.search(re_int,line).group())
				if ('p_max' in line):
					p_max  = np.int_(re.search(re_int,line).group())

					# Inialize np arrays now that dimensions are known
					block_size2 = (ml_max+1,p_max+1)
					block_size3 = (ml_max+1,p_max+1,n_vars)

					cases_run   = np.int_(  np.zeros(block_size2))
					h           = np.float_(np.zeros(block_size2))
					l2_errors   = np.float_(-float('inf')+np.zeros(block_size3))
					conv_orders = np.float_(-float('inf')+np.zeros(block_size3))

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

			self.n_vars = n_vars
			self.ml_max = ml_max
			self.p_max  = p_max

			self.var_names   = var_names
			self.cases_run   = cases_run
			self.h           = h
			self.l2_errors   = l2_errors
			self.conv_orders = conv_orders

def assign_block(A,cases_run,f):
	"""Fill data block A based on cases_run."""

	n_rows = len(cases_run)
	n_cols = len(cases_run[0])

	re_floats = '[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?'

	non_zero_rows = np.amax(cases_run,axis=1)

	for i in range(0,n_rows):
		line   = f.readline()
		vals_d = [float(s) for s in re.findall(re_floats, line)]

		j2 = 0
		for j in range(0,n_cols):
			if (cases_run[i,j] == 1):
				A[i,j] = vals_d[j]
				j2 += 1

def write_std_tables(f,data_i,output,type_out):
	"""Write output to std latex table format."""

	for n in range(0,len(data_i)):
		data = data_i[n]

		n_vars = data.n_vars
		p_max  = data.p_max
		ml_max = data.ml_max

		var_names   = data.var_names
		cases_run   = data.cases_run
		h           = data.h
		l2_errors   = data.l2_errors
		conv_orders = data.conv_orders

		f_write(f,0,'\\begin{table}[!ht]')
#		f_write(f,0,'\\begin{table}[!htbp]') # Include additional table placement options
		f_write(f,0,'\\begin{center}')
		string = ''
		if   (type_out == 'a'): string += '\\caption{Errors and Convergence Orders - '
		elif (type_out == 'e'): string += '\\caption{Errors - '
		elif (type_out == 'c'): string += '\\caption{Convergence Orders - '
		string += "{\\color{red} " + "TODO Mesh Type" + "}" + " Meshes "
		string += " {\\color{red} " + "TODO Data Type" + '}}'
		f_write(f,0,string)
		f_write(f,0,'\\resizebox{\\textwidth}{!}{')

		string = '\\begin{tabular}{| l | c | '
		if (type_out in ['a','e']):
			for k in range(0,n_vars): string += 'c '
			string += '| '
		if (type_out in ['a','c']):
			for k in range(0,n_vars): string += 'c '
			string += '| '
		string += '}'
		f_write(f,0,string)
		f_write(f,1,'\\hline')


		string  = '\\multicolumn{2}{|c|}{} & '
		if (type_out in ['a','e']):
			string += '\\multicolumn{' + '{:d}'.format(n_vars) + '}{c|}{$L^2$ Error} '
		if (type_out == 'a'):
			string += '& '
		if (type_out in ['a','c']):
			string += '\\multicolumn{' + '{:d}'.format(n_vars) + '}{c|}{Conv. Order} '
		string += '\\\\'
		f_write(f,1,string)
		f_write(f,1,'\\hline')

		string = 'Order ($p$) & Mesh Size ($h$) '
		if (type_out in ['a','e']):
			for k in range(n_vars): string += '& ' + var_names[k].replace("L^2-","") + ' '
		if (type_out in ['a','c']):
			for k in range(n_vars): string += '& ' + var_names[k].replace("L^2-","") + ' '
		string += '\\\\'
		f_write(f,1,string)
		f_write(f,1,'\\hline')

		non_zero_cols = np.amax(cases_run,axis=0)
		non_zero_rows = np.amax(cases_run,axis=1)
		for j in range(p_max+1):
			if (non_zero_cols[j] == 0):
				continue

			wrote_first_line = False
			for i in range(ml_max+1):
				if (cases_run[i][j] == 0):
					continue
				if (non_zero_rows[i] == 0):
					continue

				n_tabs = 1
				string = ''
				if (not wrote_first_line):
					string += str(j) + '\t'
					n_tabs = 0

				string += '& ' + '{:.2e}'.format(h[i,j]) + ' '
				if (type_out in ['a','e']):
					for k in range(0,n_vars):
						val_str = '{:.2e}'.format(l2_errors[i,j,k]) if not l2_errors[i,j,k] == -float('inf') else '-   '
						string += '& ' + val_str + ' '
				if (type_out in ['a','c']):
					for k in range(0,n_vars):
						val_str = '{:2.2f}'.format(conv_orders[i,j,k]) if conv_orders[i,j,k] > 0.0 else '-   '
						string += '& ' + val_str + ' '
				string += '\\\\'
				f_write(f,n_tabs,string)
				wrote_first_line = True

			f_write(f,1,'\\hline')

		f_write(f,0,'\\end{tabular}')
		f_write(f,0,'}')
		f_write(f,0,'\\end{center}')
		f_write(f,0,'\\end{table}')
		f_write(f,0,'')

def output_figure(data_i,output):
	"""Output figure in vector image format."""

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


def get_data_o (input_i,output_i):
	"""
	Return the pointer to existing data or a newly constructed subset of the input data depending on the output type.

	In the case of non-standard (std) outputs. It is assume that the various cases were run on the same meshes and for
	the same polynomial degrees.
	"""

	if (output_i.type == "std"):
		return input_i.data
	elif (output_i.type in ["ivs_adv","ivs_euler_sv","ivs_euler_sv_t","ivs_euler_gb","adv_peterson","adv_man_1d"]):
		data_i = input_i.data

		n_out  = len(data_i)
		v_out  = len(output_i.sub_index)
		ml_max = data_i[0].ml_max
		p_max  = data_i[0].p_max

		block_size3 = (ml_max+1,p_max+1,n_out*v_out)

		l2_errors   = np.zeros(block_size3)
		conv_orders = np.zeros(block_size3)

		ind_var = output_i.sub_index
		for v in range(v_out):
			for n in range(n_out):
				l2_errors[:,:,n+n_out*v]   = data_i[n].l2_errors[:,:,ind_var[v]]
				conv_orders[:,:,n+n_out*v] = data_i[n].conv_orders[:,:,ind_var[v]]

		data_o = [Input_Data(input_i.root_path)]
		data_o[0].n_vars = n_out*v_out
		data_o[0].p_max  = p_max
		data_o[0].ml_max = ml_max

		data_o[0].var_names   = output_i.var_names
		data_o[0].cases_run   = data_i[0].cases_run
		data_o[0].h           = data_i[0].h
		data_o[0].l2_errors   = l2_errors
		data_o[0].conv_orders = conv_orders

		return data_o
	else:
		assert 0,"Unsupported: "+str(output_i.format)

def output_data (output_i,input_i):
	data_o = get_data_o(input_i,output_i)
	if (output_i.format == "table"):
		f = open(output_i.root_path+"latex_tables.txt",'w')
		write_std_tables(f,data_o,output_i,'a')
		f.close()

		f = open(output_i.root_path+"latex_tables_conv_only.txt",'w')
		write_std_tables(f,data_o,output_i,'c')
		f.close()

		f = open(output_i.root_path+"latex_tables_err_only.txt",'w')
		write_std_tables(f,data_o,output_i,'e')
		f.close()
	else:
		assert 0,"Unsupported: "+str(output_i.format)

def convert_data (input_i,output_i):
	""" Read data from c code output and convert/output to appropriate format. """

	for i in range(len(input_i.rel_paths)):
		input_i.data.append(Input_Data(input_i.return_data_path(i)))
		input_i.data[i].read_data()

	# Output data (latex table or eps figure)
	output_data(output_i,input_i)

if __name__ == '__main__':
	"""
	This function serves to convert data generated by the code into output files which can be used in LaTeX. It is
	expected that a file called `l2errs+convergence_standard.txt` is present in each subfolder of the input data
	paths.

	Command line arguments:
	1. Path to the ROOT directory containing error data files.
	2. The error data file type to be processed. Options:
	- std:          generate data for a single data file.
	- ivs_adv:      'i'soparametric 'v's. 's'uperparametric linear advection case data.
	- ivs_euler_sv: 'i'soparametric 'v's. 's'uperparametric euler ('s'upersonic 'v'ortex)
	- ivs_euler_sv_t: 'i'soparametric 'v's. 's'uperparametric euler ('s'upersonic 'v'ortex) 't'ransonic
	- ivs_euler_gb: 'i'soparametric 'v's. 's'uperparametric euler ('g'aussian 'b'ump)
	- adv_peterson: Linear "Adv"ection "Peterson" test case.
	- adv_man_1d:   Linear "Adv"ection 1d "Man"ufacture test case.
	3. The dimension of the problem.
	"""

	assert len(sys.argv) == 4,"\nIncorrect number of inputs. Should be:\n"\
	                         +"\t1. full_path_to_error_file.\n"\
	                         +"\t2. error_file_type.\n"\
	                         +"\t3. dimension\n\n."

	input_path = sys.argv[1]
	data_type  = sys.argv[2]
	dimension  = sys.argv[3]
	input_i = Input_Info(input_path,data_type)

	output_i = Output_Info(input_path,data_type,dimension)

	print("\nOutputting data in {} format.\n\n".format(output_i.format))
	convert_data(input_i,output_i)
