import sys
import subprocess
import shlex
import re

class Gmsh_call:
	""" Container for gmsh call related information. """
	def __init__ (self):
		self.dim      = '' ###< The dimension.
		self.in_dir   = '' ###< The full path to the ROOT input directory.
		self.in_name  = '' ###< The name of the input file.
		self.out_name = '' ###< The name of the output file.
		self.args     = '' ###< The arguments to be passed to gmsh

	def set_input (self,project_src_dir,mesh_name):
		self.in_dir  = project_src_dir+'/input/meshes/'
		self.in_name = self.in_dir

		supported_geometries = ['n-cube',
		                        'n-cylinder_hollow_section',
		                       ]

		geometry = mesh_name.split('/')[0]
		if (geometry not in supported_geometries):
			print("The ",geometry," is not currently supported.\n")
			EXIT

		geo_name = re.search(r"(__)([\w]+)(__)",mesh_name).group(2)

		self.in_name += geometry+'/'+geo_name+".geo"

		self.dim = re.search(r"(^[\D]*)(\d)([\D]*$)",geo_name).group(2)

	def set_output (self,mesh_name_full):
		self.out_name = mesh_name_full

	def set_args (self,mesh_name,mesh_name_full):
		self.args  = ' ' + self.in_name
		self.args += " -" + self.dim
		self.args += set_gmsh_setnumbers(self.in_dir,mesh_name)
		self.args += ' -o ' + self.out_name

		out_dir = re.search(r"/(([\w-]+/)*)",mesh_name_full).group(0)

		subprocess.call(shlex.split('mkdir -p ' + out_dir))
		subprocess.call(shlex.split('gmsh' + self.args))

def set_gmsh_setnumbers (input_dir,mesh_name):
	""" Set the -setnumber inputs to be passed to gmsh. """
	gmsh_setnumbers = ''

	# Mesh level

	mesh_level = re.search(r"(^.*_ml)(\d+)(.*$)",mesh_name).group(2)
	gmsh_setnumbers += ' -setnumber MESH_LEVEL ' + mesh_level

	# Required parameters

	gmsh_setnumbers += ' -setnumber PDE_NAME '
	var_names = ['advection','poisson','euler','navierstokes']
	gmsh_setnumbers += get_gmsh_number(mesh_name,var_names,input_dir,0)

	gmsh_setnumbers += ' -setnumber PDE_SPEC '
	var_names = ['internal/supersonic_vortex',
	             'periodic/periodic_vortex',
	            ]
	gmsh_setnumbers += get_gmsh_number(mesh_name,var_names,input_dir,0)

	gmsh_setnumbers += ' -setnumber MESH_DOMAIN '
	var_names = ['straight','curved','parametric']
	gmsh_setnumbers += get_gmsh_number(mesh_name,var_names,input_dir,0)

	gmsh_setnumbers += ' -setnumber MESH_TYPE '
	var_names = ['line','tri','quad','tet','hex','wedge','pyr','mixed']
	gmsh_setnumbers += get_gmsh_number(mesh_name,var_names,input_dir,1)

	# Additional geom_spec parameters

	return gmsh_setnumbers


# Delete other function and change name of this function
def get_gmsh_number (mesh_name,var_names,input_dir,with_underscore):
	""" Get the number associated with the variable name as specified in the parameters.geo file. """

	param_file_name = input_dir+'/parameters.geo'

	for target in var_names:
		target_name = target
		if (with_underscore):
			target_name = '_'+target_name+'_'

		if (mesh_name.find(target_name) != -1):
			var_name = target.replace('/','_')
			break

	with open(param_file_name) as f:
		for line in f:
			if (var_name.upper() in line):
				return line.split()[2][:-1]

	print("\n\nDid not find a value for "+var_name.upper()+" in "+param_file_name+'\n')
	EXIT



if __name__ == '__main__':
	""" Generate the gmsh mesh passed as a command line argument using the appropriate .geo file. """


	if (len(sys.argv) > 3):
		print("Only a single mesh_name should be input at a time.\n")
		EXIT

	project_src_dir = sys.argv[1]
	mesh_name_full  = sys.argv[2]

	mesh_name = re.sub(".*meshes/","",mesh_name_full)

	parts = mesh_name.split('/')
	print(parts)

	gmsh_call = Gmsh_call()
	gmsh_call.set_input(project_src_dir,mesh_name)
	gmsh_call.set_output(mesh_name_full)
	gmsh_call.set_args(mesh_name,mesh_name_full)
