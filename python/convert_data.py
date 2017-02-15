import sys
import re
import numpy as np
import matplotlib.pyplot as plt

from numpy import linalg as LA

#np.set_printoptions(precision=3)
np.set_printoptions(linewidth=200)
np.set_printoptions(formatter={'float': lambda x: format(x, '6.2e')})

class output_class:
	"""Stores output related information."""

	def init_name(self,name):
		if (name == 'std'):
			self.NEntries = 1
			self.type   = 'std'
			self.format = 'table'
		elif (name == 'Optimal Surface'):
			self.NEntries = 5
			self.type     = 'Poisson_H0'
			self.format   = 'table'
#			self.format   = 'figure'
		elif (name == 'Optimal Blending'):
			self.NEntries = 4
			self.type     = 'Poisson_H0'
			self.format   = 'table'
		elif (name == 'Suboptimal Surface'):
			self.NEntries = 4
			self.type     = 'Poisson_H0'
			self.format   = 'table'
		elif (name == 'Suboptimal Blending'):
			self.NEntries = 3
			self.type     = 'Poisson_H0'
			self.format   = 'table'
		elif (name == 'Optimal h'):
			self.NEntries = 5
			self.type     = 'Poisson_H0'
			self.format   = 'table'
		elif (name == 'Projected Conforming'):
			self.NEntries = 4
			self.type     = 'Poisson_H0'
			self.format   = 'table'
		elif (name == 'Projected NonConforming'):
			self.NEntries = 5
			self.type     = 'Poisson_H0'
			self.format   = 'table'
		elif (name == 'Optimal Euler'):
			self.NEntries = 4
			self.type     = 'Euler_s'
			self.format   = 'table'
		else:
			print("Error: Unsupported (output_class)."); sys.exit()


class input_data:
	"""Stores input data related information."""

	def __init__(self):
		self.case  = list()
		self.NVars = 0
		self.MLMax = 0
		self.PMax  = 0


def initialize_input(output):
	"""Initializes input data based on the desired output."""

	data_i = [ input_data() for i in range(0,output.NEntries)]

	if (output.name == 'std'):
		data_i[0].case = 'L2errs+Convergence'
	elif (output.name == 'Optimal Surface'):
		name_root = 'L2errs+Convergence_' + output.Geometry + '_' + output.MeshType
		data_i[0].case = name_root + '_L2s'
		data_i[1].case = name_root + '_L2_SB'
		data_i[2].case = name_root + '_SB_Normal'
		data_i[3].case = name_root + '_SB_ArcLength'
		data_i[4].case = name_root + '_SB_Radial'

		output.varNames = ['$L^2$ (s)','$L^2$ (SB)','SB$_{\\text{N}}$','SB$_{\\text{A}}$','SB$_{\\text{R}}$']
	elif (output.name == 'Optimal Blending'):
		name_root = 'L2errs+Convergence_' + output.Geometry + '_' + output.MeshType
		data_i[0].case = name_root + '_L2_SB'
		data_i[1].case = name_root + '_SB_Normal'
		data_i[2].case = name_root + '_Scott_Normal'
		data_i[3].case = name_root + '_Lenoir_Normal'

		output.varNames = ['$L^2$ (SB)','SB$_{\\text{N}}$','Scott$_{\\text{N}}$','Lenoir$_{\\text{N}}$']
	elif (output.name == 'Suboptimal Surface'):
		name_root = 'L2errs+Convergence_' + output.Geometry + '_' + output.MeshType
		data_i[0].case = name_root + '_SB_Radial'
		data_i[1].case = name_root + '_SB_Order_h_Radial'
		data_i[2].case = name_root + '_SB_Order_h_i3c-1'
		data_i[3].case = name_root + '_SB_Order_h_i2c-1'

		output.varNames = ['SB$_{\\text{R}}$','SB$_{\\text{H$_{\\text{R}}$}}$','SB$_{\\text{H$_{3,-1}$}}$','SB$_{\\text{H$_{2,-1}$}}$']
	elif (output.name == 'Suboptimal Blending'):
		name_root = 'L2errs+Convergence_' + output.Geometry + '_' + output.MeshType
		data_i[0].case = name_root + '_SB_Normal'
		data_i[1].case = name_root + '_Nielson_Normal'
		data_i[2].case = name_root + '_Nielson_Radial'

		output.varNames = ['SB$_{\\text{N}}$','Nielson$_{\\text{N}}$','Nielson$_{\\text{R}}$']
	elif (output.name == 'Optimal h'):
		name_root = ['L2errs+Convergence_' + output.Geometry + output.Ratio + '_' + output.MeshType, \
		             'L2errs+Convergence_' + output.Geometry + output.Ratio + '_' + 'CurvedQUAD']

		data_i[0].case = name_root[0] + '_L2_m=1-ns'
		data_i[1].case = name_root[0] + '_L2_m=1'
		data_i[2].case = name_root[0] + '_L2_m=k'
		data_i[3].case = name_root[0] + '_Normal_SB_m=k'
		data_i[4].case = name_root[1] + '_Normal_GH_m=k'

		output.varNames = ['$L^2_{m=1,\\text{ns}}$','$L^2_{m=1}$','$L^2_{m=k}$','SB$_{N,m=k}$','GH$_{N,m=k}$']
	elif (output.name == 'Projected Conforming'):
		name_root = ['L2errs+Convergence_' + output.Geometry + '3_ProjectedCurvedTRI', \
		             'L2errs+Convergence_' + output.Geometry + '3_ProjectedFullyCurvedTRI']

		name_end  = ['_L2_m=k', '_Normal_SB_m=k']

		data_i[0].case = name_root[0] + name_end[0]
		data_i[1].case = name_root[0] + name_end[1]
		data_i[2].case = name_root[1] + name_end[0]
		data_i[3].case = name_root[1] + name_end[1]

		output.varNames = ['$L^2_{m=k,\\text{ps}}$','SB$_{N,m=k,ps}$','$L^2_{m=k,\\text{pc}}$','SB$_{N,m=k,pc}$']
	elif (output.name == 'Projected NonConforming'):
		name_root = ['L2errs+Convergence_' + output.Geometry + '3_NonConforming_CurvedTRI', \
		             'L2errs+Convergence_' + output.Geometry + '3_NonConforming_ProjectedCurvedTRI', \
		             'L2errs+Convergence_' + output.Geometry + '3_NonConforming_ProjectedFullyCurvedTRI']

		name_end  = ['_L2_m=k', '_Normal_SB_m=k']

		data_i[0].case = name_root[0] + name_end[1]
		data_i[1].case = name_root[1] + name_end[0]
		data_i[2].case = name_root[1] + name_end[1]
		data_i[3].case = name_root[2] + name_end[0]
		data_i[4].case = name_root[2] + name_end[1]

		output.varNames = ['SB$_{N,m=k}$','$L^2_{m=k,\\text{ps}}$','SB$_{N,m=k,ps}$','$L^2_{m=k,\\text{pc}}$', \
		                   'SB$_{N,m=k,pc}$']
	elif (output.name == 'Optimal Euler'):
		name_root = ['L2errs+Convergence_' + output.Geometry + '2.25_' + output.MeshType, \
		             'L2errs+Convergence_' + output.Geometry + '1.00_' + output.MeshType]

		name_end  = ['original_Normal_SB_m=k','_Normal_SB_m=k'];

		data_i[0].case = name_root[0] + name_end[0]
		data_i[1].case = name_root[0] + name_end[1]
		data_i[2].case = name_root[1] + name_end[0]
		data_i[3].case = name_root[1] + name_end[1]

		output.varNames = ['SB$_{N,m=k,(2.25)}$','SB$_{N,m=k,(2.25r)}$','SB$_{N,m=k,(1.00)}$','SB$_{N,m=k,(1.00r)}$']
	else:
		print("Error: Unsupported (initialize_input)."); sys.exit()

	return data_i

def skip_lines(N,f):
	for n in range(0,N):
		line = f.readline()

def assign_block(A,CasesRun,f):
	"""Fill data block A based on CasesRun"""

	NRows = len(CasesRun)
	NCols = len(CasesRun[0])

	re_floats = '[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?'

	NonZeroRows = np.amax(CasesRun,axis=1)

	for i in range(0,NRows):
		if (NonZeroRows[i] != 1):
			continue

		line   = f.readline()
		vals_d = [float(s) for s in re.findall(re_floats, line)]

		j2 = 0
		for j in range(0,NCols):
			if (CasesRun[i,j] == 1):
				A[i,j] = vals_d[j2]
				j2 += 1


def read_data(data_i):
	"""Reads in data required by specific output."""

	re_int    = '\d+'
	re_ints   = '\\b\d+\\b'
	re_floats = '[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?'

	for n in range(0,len(data_i)):
		data = data_i[n]

		fName = data.case + '.txt'

		varNames = list()

		with open(fName) as f:
			for line in f:
				if ('NVars' in line):
					NVars = np.int_(re.search(re_int,line).group())
				if ('MLMax' in line):
					MLMax = np.int_(re.search(re_int,line).group())
				if ('PMax' in line):
					PMax  = np.int_(re.search(re_int,line).group())

					# Inialize np arrays once dimensions are known
					block_size2 = (MLMax+1,PMax+1)
					block_size3 = (MLMax+1,PMax+1,NVars)

					CasesRun   = np.int_(  np.zeros(block_size2))
					h          = np.float_(np.zeros(block_size2))
					L2Errors   = np.float_(np.zeros(block_size3))
					ConvOrders = np.float_(np.zeros(block_size3))

				if ('Cases Run' in line):
					for i in range(0,MLMax+1):
						line = f.readline()
						CasesRun[i,:] = np.int_([int(s) for s in re.findall(re_ints, line)])

				if ('Mesh Size' in line):
					assign_block(h,CasesRun,f)

				if ('L2Errors' in line):
					for k in range(0,NVars):
						skip_lines(1,f)
						line = f.readline()
						varNames.append(line.replace('\n',''))

						assign_block(L2Errors[:,:,k],CasesRun,f)

				if ('Convergence Orders' in line):
					for k in range(0,NVars):
						skip_lines(2,f)
						assign_block(ConvOrders[:,:,k],CasesRun,f)


			data.NVars = NVars
			data.MLMax = MLMax
			data.PMax  = PMax

			data.varNames   = varNames
			data.CasesRun   = CasesRun
			data.h          = h
			data.L2Errors   = L2Errors
			data.ConvOrders = ConvOrders

def f_write(f,Ntabs,string):
	"""Write string to file with specified number of tabs and newline"""
	for n in range(0,Ntabs):
		f.write('\t')
	f.write(string + '\n')

def write_std_tables(f,data_i,output):
	"""Write output to std latex table format."""

	for n in range(0,len(data_i)):
		data = data_i[n]

		NVars = data.NVars
		PMax  = data.PMax
		MLMax = data.MLMax

		varNames   = data.varNames
		CasesRun   = data.CasesRun
		h          = data.h
		L2Errors   = data.L2Errors
		ConvOrders = data.ConvOrders

		f_write(f,0,'\\begin{table}[!ht]')
#		f_write(f,0,'\\begin{table}[!htbp]') # Include additional table placement options
		f_write(f,0,'\\begin{center}')
		string = '\\caption{Errors and Convergence Orders - ' + output.MeshType + ' Meshes {\\color{red} '
		string += output.name + '}}'
		f_write(f,0,string)
		f_write(f,0,'\\resizebox{\\textwidth}{!}{')

		string = '\\begin{tabular}{| l | c | '
		for k in range(0,NVars): string += 'c '
		string += '| '
		for k in range(0,NVars): string += 'c '
		string += '| }'
		f_write(f,0,string)
		f_write(f,1,'\\hline')

		string  = '\\multicolumn{2}{|c|}{} & \\multicolumn{' + '{:d}'.format(NVars) + '}{c|}{$L^2$ Error} & '
		string +=                           '\\multicolumn{' + '{:d}'.format(NVars) + '}{c|}{Conv. Order} \\\\'
		f_write(f,1,string)
		f_write(f,1,'\\hline')

		string = 'Order ($k$) & Mesh Size ($h$) '
		for k in range(0,NVars): string += '& ' + varNames[k] + ' '
		for k in range(0,NVars): string += '& ' + varNames[k] + ' '
		string += '\\\\'
		f_write(f,1,string)
		f_write(f,1,'\\hline')

		NonZeroCols = np.amax(CasesRun,axis=0)
		NonZeroRows = np.amax(CasesRun,axis=1)
		for j in range(0,PMax+1):
			if (NonZeroCols[j] == 0):
				continue

			for i in range(0,MLMax+1):
				if (NonZeroRows[i] == 0):
					continue

				Ntabs = 1
				string = ''
				if i == 0:
					string += str(j) + '\t'
					Ntabs = 0

				string += '& ' + '{:.2e}'.format(h[i,j]) + ' '
				for k in range(0,NVars): string += '& ' + '{:.2e}'.format(L2Errors[i,j,k]) + ' '

				if i == 0:
					for k in range(0,NVars): string += '& -    '
				else:
					for k in range(0,NVars): string += '& ' + '{:2.2f}'.format(ConvOrders[i,j,k]) + ' '
				string += '\\\\'
				f_write(f,Ntabs,string)

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

		NVars = data.NVars
		PMax  = data.PMax
		MLMax = data.MLMax

		varNames   = data.varNames
		CasesRun   = data.CasesRun
		h          = data.h
		L2Errors   = data.L2Errors
		ConvOrders = data.ConvOrders

		NonZeroCols = np.amax(CasesRun,axis=0)
		NonZeroRows = np.amax(CasesRun,axis=1)

		ind_i = np.arange(0)
		for i in range(0,MLMax+1):
			if (NonZeroRows[i] == 0):
				ind_i = np.append(ind_i,[i])
		ind_i = np.delete(np.arange(MLMax+1),ind_i)

		for j in range(0,PMax+1):
			if (NonZeroCols[j] == 0):
				continue

			data_x = h[ind_i,j]
#			print(data_x)

			for k in range(0,NVars):
				data_y = L2Errors[ind_i,j,k]

		# varNames go in legend
#				print(data_y)
				plt.plot(data_x,data_y,'-o')

#			print(" ")


	plt.xscale('log')
	plt.yscale('log')
	axes = plt.gca()
	axes.set_xlim([1.05*np.amin(h),1.05*np.amax(h)])
	plt.savefig('latex_figure.eps', format='eps', dpi=1200)
	sys.exit()


def output_data(output,data_i):

	if ('table' in output.format):
		f = open('latex_tables.txt', 'w')

	if ('std' in output.type):
		write_std_tables(f,data_i,output)
	elif ('Poisson_H0' in output.type or 'Euler_s' in output.type):
		# Assemble data to print
		NOut  = len(data_i)
		MLMax = data_i[0].MLMax
		PMax  = data_i[0].PMax

		block_size3 = (MLMax+1,PMax+1,NOut)

		L2Errors   = np.zeros(block_size3)
		ConvOrders = np.zeros(block_size3)

		for n in range(0,NOut):
			data = data_i[n]

			L2Errors[:,:,n]   = data.L2Errors[:,:,0]
			ConvOrders[:,:,n] = data.ConvOrders[:,:,0]

		data_print = [ input_data() for i in range(0,1)]
		data_print[0].NVars = NOut
		data_print[0].PMax  = PMax
		data_print[0].MLMax = MLMax

		data_print[0].varNames   = output.varNames
		data_print[0].CasesRun   = data_i[1].CasesRun
		data_print[0].h          = data_i[1].h
		data_print[0].L2Errors   = L2Errors
		data_print[0].ConvOrders = ConvOrders

		data_print[0].MeshType   = output.MeshType
		data_print[0].case       = data_i[0].case

		if ('table' in output.format):
			write_std_tables(f,data_print,output)
		elif ('figure' in output.format):
			output_figure(data_print,output)
			print("Add support"); sys.exit()
#			plt.savefig('latex_figure.eps', format='eps', dpi=1200)
#			fig.savefig('latex_figure.svg', format='svg', dpi=1200)
		else:
			print("Error: Unsupported (output.format)."); sys.exit()
	else:
		print("Error: Unsupported (output.type)."); sys.exit()

	if ('table' in output.format):
		f.close()
	



def convert_data(output):
	"""Convert data from c code input to appropriate output format"""

	# Initialize input data class
	data_i = initialize_input(output)

	# Read in and store data
	read_data(data_i)

	# Output data (latex table or eps figure)
	output_data(output,data_i)



if __name__ == '__main__':
	output = output_class()

#	output.Geometry = 'dm1-Spherical_Section'
#	output.Geometry = 'Ellipsoidal_Section'
	output.Geometry = 'InviscidChannel_Joukowski'

	output.MeshType = 'CurvedTRI'
#	output.MeshType = 'CurvedQUAD'

#	output.name = 'std'
#	output.name = 'Optimal Surface'
#	output.name = 'Optimal Blending'
#	output.name = 'Suboptimal Surface'
#	output.name = 'Suboptimal Blending'
#	output.name = 'Optimal h'; output.Ratio = '1'
#	output.name = 'Projected Conforming'
#	output.name = 'Projected NonConforming'
	output.name = 'Optimal Euler';

	output.init_name(output.name)

	print("Outputting data in {} format".format(output.format))
	convert_data(output)
