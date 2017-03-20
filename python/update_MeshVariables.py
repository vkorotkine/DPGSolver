import sys


### Classes ###
from meshfile_classes import Paths_class
from meshfile_classes import TestCase_class


### Functions ###
from support_functions import f_write

def write_MeshVariables(TestCases,Paths):
	"""Write MeshVariables file."""

	f = open(Paths.meshes+'MeshVariables', 'w')

	f_write(f,0,2,'# Variables and dependencies needed for DPG_ROOT/meshes/Makefile')

	for i in range(0,len(TestCases)):
		TestCase = TestCases[i]
		f_write(f,0,0,TestCase.VarName+' := ')
		if (TestCase.name.find('L2_proj_h') != -1):
			print("Not including "+TestCase.VarName+" mesh files as they are duplicates of those for L2_PROJ_P.")
			f_write(f,0,2,"")
		else:
			f_write(f,0,2,TestCase.MeshOutputs)

	f_write(f,0,2,'')

	for i in range(0,len(TestCases)):
		TestCase = TestCases[i]
		f_write(f,0,2,'$('+TestCase.VarName+') : '+TestCase.GeoDeps)


if __name__ == '__main__':
	"""Generate the MeshVariables file based on .ctrl files specified in the Makefile of the ROOT directory."""
	user = 'PZwan'

	Paths = Paths_class()
	Paths.set_paths(user)

	SubDirectories = sys.argv[1:]

	TestCases = []

	for i in range(0,len(SubDirectories)):
		TestCase = TestCase_class(SubDirectories[i])
		TestCase.set_paths(Paths)
		TestCase.add_MeshTypes(Paths,'all')
		
		TestCase.get_geo_dependencies()
		TestCase.get_mesh_outputs()

		TestCases.append(TestCase)

	# Output MeshVariables file
	write_MeshVariables(TestCases,Paths)
