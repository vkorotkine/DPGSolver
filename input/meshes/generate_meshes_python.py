import sys
import subprocess
import shlex
import re

class Python_call:
	""" Container for python3 call related information. """
	def __init__ (self,project_src_dir,mesh_name_full):
		self.in_dir   = '' ###< The full path to the ROOT input directory.
		self.in_name  = '' ###< The name of the input file.
		self.out_name = '' ###< The name of the output file.

		### The arguments to be passed to python function
		self.args     = project_src_dir + ' ' + mesh_name_full

	def set_input (self,project_src_dir,mesh_name):
		self.in_dir  = project_src_dir+"/input/meshes/"
		self.in_name = self.in_dir

		supported_geometries = ["n-cube",
		                       ]

		geometry = mesh_name.split('/')[0]
		if (geometry not in supported_geometries):
			print("The ",geometry," is not currently supported.\n")
			EXIT

		geo_name = re.search(r"(__)([\w]+)(__)",mesh_name).group(2)

		self.in_name += geometry+'/'+geo_name+".py"

	def set_output (self,mesh_name_full):
		self.out_name = mesh_name_full

	def call_function (self,mesh_name_full):
		out_dir = re.search(r"/(([\w-]+/)*)",mesh_name_full).group(0)
		subprocess.call(shlex.split("mkdir -p " + out_dir))

		subprocess.call(shlex.split("python3 " + self.in_name + ' ' + self.args))




if __name__ == "__main__":
	""" Generate the python mesh passed as a command line argument using the appropriate .py file. """

	if (len(sys.argv) > 3):
		print("Only a single mesh_name should be input at a time.\n")
		EXIT

	project_src_dir = sys.argv[1]
	mesh_name_full  = sys.argv[2]
	mesh_name       = re.sub(".*meshes/","",mesh_name_full)

	python_call = Python_call(project_src_dir,mesh_name_full)

	python_call.set_input(project_src_dir,mesh_name)
	python_call.set_output(mesh_name_full)
	python_call.call_function(mesh_name_full)
