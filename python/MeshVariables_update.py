import sys


### Classes ###
from meshfile_classes import Paths_class
from meshfile_classes import TestCase_class


### Functions ###
from support_functions import f_write

def write_MeshVariables(TestCases,Paths):
	"""Write MeshVariables file."""

	f = open(Paths.meshes+'MeshVariables', 'w')

	f_write(f,0,1,'# Variables and dependencies needed for DPG_ROOT/meshes/Makefile')
	f_write(f,0,2,'# Automatically generated (DO NOT MODIFY)')

	for i in range(0,len(TestCases)):
		TestCase = TestCases[i]
		f_write(f,0,2,TestCase.VarName+' := '+TestCase.MeshOutputs)

	f_write(f,0,2,'')

	for i in range(0,len(TestCases)):
		TestCase = TestCases[i]
		f_write(f,0,2,'$('+TestCase.VarName+') : '+TestCase.GeoDeps)

	f.close()


if __name__ == '__main__':
	"""Generate the MeshVariables file based on .ctrl files specified in the Makefile of the ROOT directory."""
	Paths = Paths_class()
	Paths.set_paths()

	SubDirectories = sys.argv[1:]

#	print("SD:",SubDirectories)
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
