"""
Generate Peterson meshes used for testing sharpness of suboptimal convergence of the DG method.

References:
	Peterson(1991)-A Note on the Convergence of the Discontinuous Galerkin Method for a Scalar Hyperbolic Equation
"""

import sys
import os
import math
import numpy as np
import re

sys.path.insert(0,"../../../src/processing/")
from support_functions import list_print

### Classes ###
class Mesh_c:
	"""
	Defines mesh related information.
	"""

	def __init__(self,sigma,ML,mesh_name_full):
		self.sigma          = sigma
		self.ML             = ML
		self.mesh_name_full = mesh_name_full

		h = math.pow(0.5,ML)
		m = round(math.pow(h,-sigma))

		self.NRows = int(math.pow(2,ML+1)+1)
		self.NCols = 2+m-1

		NRows = self.NRows
		self.NZCols     = np.zeros(NRows,dtype=np.int)
		self.NZCols_evn = np.zeros(NRows,dtype=np.int)
		self.NZCols_odd = np.zeros(NRows,dtype=np.int)

		self.node_index = np.zeros(NRows*NRows,dtype=np.int)

		self.coordinates = np.zeros(shape=(NRows*NRows,3))

		self.boundaries = []
		self.conn = []
		self.gmsh_elements = []

	def compute_nonzero_cols(self):
		ensure_no_intersection = 0 # Ensure that vertical lines do not intersect boundary elements

		Interval = (self.NRows-1)/(self.NCols-1)

		NZCols = self.NZCols
		for i in range(0,self.NCols):
			Index = round(i*Interval)
			if (ensure_no_intersection and (Index % 2 != 0)):
				if (Index < self.NRows/2):
					Index += 1
				else:
					Index -= 1
			NZCols[Index] = 1

		for i in range(0,self.NRows):
			if (NZCols[i] == 1):
				self.NZCols_evn[i] = 1
				self.NZCols_odd[i] = 1
			else:
				if (i % 2 == 0):
					self.NZCols_evn[i] = 1
				else:
					self.NZCols_odd[i] = 1

	def compute_global_node_list(self):
		node_index = self.node_index
		IndG = 0
		for row in range(0,self.NRows):
			if (row % 2 == 0):
				NZCols = self.NZCols_evn
			else:
				NZCols = self.NZCols_odd

			for col in range(0,self.NRows):
				IndT = row*self.NRows+col

				if (NZCols[col] == 0):
					node_index[IndT] = -1
				else:
					node_index[IndT] = IndG
					IndG += 1

	def compute_coordinates(self):
		coordinates = self.coordinates

		NRows = self.NRows

		IndG = 0
		for row in range(0,NRows):
			y = -1.0 + row/(NRows-1)*2.0
			for col in range(0,NRows):
				x = -1.0 + col/(NRows-1)*2.0

				coordinates[IndG][:] = [x,y,0.0]
				IndG += 1

	def compute_connectivity(self):
		def add_local_conn(self,neighbour_location,conn_E):
			"""
			Add local (E)lement (conn)ectivity array to global connectivity array. Swap the position of TRI nodes to
			avoid inverted elements.
			"""

			if (neighbour_location == 'above'):
				self.conn.append(conn_E)
			else:
				conn_E[1], conn_E[2] = conn_E[2], conn_E[1]
				self.conn.append(conn_E)

		def add_TRIs_conn(self,row,neighbour_location):
			"""
			Add TRI element connectivity entries to the global connectivity array for the current row.
			"""

			conn   = self.conn
			if (row % 2 == 0):
				NZCols = self.NZCols_evn
			else:
				NZCols = self.NZCols_odd

			base_c = row*self.NRows
			if (neighbour_location == 'above'):
				base_n = (row+1)*self.NRows
			elif (neighbour_location == 'below'):
				base_n = (row-1)*self.NRows
			else:
				EXIT

			iterator = iter(range(0,self.NRows-1))
			for col in iterator:
				if (NZCols[col] == 0):
					continue

				if (row % 2 != 0):
					if (col == 0):
						conn_E = [base_c+col,base_c+col+1,base_n+col]
						add_local_conn(self,neighbour_location,conn_E)
						continue
					elif (col == self.NRows-2):
						conn_E = [base_c+col,base_c+col+1,base_n+col+1]
						add_local_conn(self,neighbour_location,conn_E)
						continue

				if (NZCols[col+1] == 0):
					conn_E = [base_c+col,base_c+col+2,base_n+col+1]
					add_local_conn(self,neighbour_location,conn_E)
				else:
					conn_E = [base_c+col,base_c+col+1,base_n+col+1]
					add_local_conn(self,neighbour_location,conn_E)

					conn_E = [base_c+col+1,base_c+col+2,base_n+col+1]
					add_local_conn(self,neighbour_location,conn_E)
					next(iterator,None)

		def renumber_conn(self):
			node_index = self.node_index
			conn       = self.conn

			for item in conn:
				for i in range(0,len(item)):
					item[i] = node_index[item[i]]


		# Compute connectivity with indices assuming that all nodes are present
		for row in range(0,self.NRows):
			if (row != 0):
				add_TRIs_conn(self,row,'below')
			if (row != self.NRows-1):
				add_TRIs_conn(self,row,'above')

		# Update connectivity indices deleting those which are not present
#		list_print(self.conn,"conn",2)
		renumber_conn(self)
#		print(node_index)
#		list_print(self.conn,"conn",2)

	def compute_boundaries(self,input_dir):
		"""
		Compute the indices of the boundary elements and assign their gmsh tags.
		"""

		def add_boundaries(self,line_index,BC_index):
			boundaries = self.boundaries

			NRows = self.NRows
			if (math.floor(line_index / 1000) == 1):
				# Horizontal boundaries (some entries skipped over depending on NZCols)
				NZCols = self.NZCols_evn

				base_c = 0
				if (line_index == 1002):
					base_c = NRows*(NRows-1)

				for col in range(0,NRows-1):
					if (NZCols[col] == 0):
						continue

					if (NZCols[col+1] == 0):
						conn_E = [base_c+col,base_c+col+2]
					else:
						conn_E = [base_c+col,base_c+col+1]

					boundaries.append([2,10000+BC_index,line_index]+conn_E)

			elif (math.floor(line_index / 1000) == 2):
				# Vertical boundaries
				base_c = 0
				if (line_index == 2002):
					base_c = NRows-1

				NRows = self.NRows
				for row in range(0,NRows-1):
					conn_E = [base_c+row*NRows,base_c+(row+1)*NRows]

					boundaries.append([2,10000+BC_index,line_index]+conn_E)
			else:
				EXIT

		def renumber_boundaries(self):
			node_index = self.node_index
			boundaries = self.boundaries

			for item in boundaries:
				for i in range(len(item)-2,len(item)):
					item[i] = node_index[item[i]]

		def get_gmsh_number (var_name,input_dir):
			param_file_name = input_dir+"/parameters.geo"

			with open(param_file_name) as f:
				for line in f:
					if (var_name.upper() in line):
						return line.split()[2][:-1]

			print("\n\nDid not find a value for "+var_name.upper()+" in "+param_file_name+'\n')
			EXIT


		b_conditions = dict()
		b_c_members = ["bc_inflow","bc_outflow"]

		for s in b_c_members:
			b_conditions[s] = get_gmsh_number(s,input_dir)

		add_boundaries(self,1001,int(b_conditions["bc_inflow"]))
		add_boundaries(self,1002,int(b_conditions["bc_outflow"]))
		add_boundaries(self,2001,int(b_conditions["bc_outflow"]))
		add_boundaries(self,2002,int(b_conditions["bc_outflow"]))

		renumber_boundaries(self)

#		list_print(self.boundaries,"",2)

	def compute_gmsh_elements_array(self):
		boundaries = self.boundaries
		conn       = self.conn

		gmsh_elements = self.gmsh_elements

		IndE = 1

		# LINEs
		for item in boundaries:
			for i in range(len(item)-2,len(item)):
				item[i] += 1
			gmsh_elements.append([IndE,1]+item)
			IndE += 1

		# TRIs
		for item in conn:
			for i in range(0,len(item)):
				item[i] += 1
			gmsh_elements.append([IndE,2,2,9401,4001]+item)
			IndE += 1

#		list_print(gmsh_elements,"",2)

	def output(self):
		def f_write(f,string):
			f.write(string + '\n')

		mesh_name = self.mesh_name_full
		mesh_dir = re.search(r"/(([\w-]+/)*)",mesh_name).group(0)

		if not os.path.exists(mesh_dir):
			os.makedirs(mesh_dir)

		print("Generating mesh: ",mesh_name)
		f = open(mesh_name,'w')

		# Header
		f_write(f,'$MeshFormat')
		f_write(f,'2.2 0 8')
		f_write(f,'$EndMeshFormat')

		# Nodes
		coordinates = self.coordinates
		node_index = self.node_index

		f_write(f,'$Nodes')
		f_write(f,str(sum(x >= 0 for x in node_index)))
		for n in range(0,coordinates.shape[0]):
			if (node_index[n] < 0):
				continue
			Indn = node_index[n]
			f_write(f,str(Indn+1) + ' ' + ' '.join([str(x) for x in coordinates[n]]))
		f_write(f,'$EndNodes')

		# Elements
		gmsh_elements = self.gmsh_elements

		f_write(f,'$Elements')
		f_write(f,str(len(gmsh_elements)))
		for item in gmsh_elements:
			f_write(f,' '.join([str(x) for x in item]))
		f_write(f,'$EndElements')

		f.close()


### Functions ###

if __name__ == '__main__':
	"""
	Generate meshes and export in the gmsh format.

	Reference: Richter(2008)-On the Order of Convergence of the Discontinuous Galerkin Method for Hyperbolic Equations
	"""

	project_src_dir = sys.argv[1]
	input_dir       = project_src_dir+"/input/meshes"
	mesh_name_full  = sys.argv[2]

	NML   = 6
	NP    = 3

	sigma = 0.0
	for p in range(0,NP+1):
		mesh_name_full = re.sub("/p\d_","/p"+str(p)+"_",mesh_name_full)

		# The value of sigma is taken from (eq. (27), Richter(2008)).
		sigma = (2.0*p+1.0)/(2.0*p+2.0)

		for ML in range(0,NML):
			mesh_name_full = re.sub("_ml\d.msh","_ml"+str(ML)+".msh",mesh_name_full)

			Mesh = Mesh_c(sigma,ML,mesh_name_full)

			Mesh.compute_nonzero_cols()
			Mesh.compute_global_node_list()

			Mesh.compute_coordinates()

			Mesh.compute_boundaries(input_dir)
			Mesh.compute_connectivity()
			Mesh.compute_gmsh_elements_array()

			Mesh.output()
