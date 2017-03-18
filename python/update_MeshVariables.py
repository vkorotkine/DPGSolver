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
		if (TestCase.name.find('update_h') != -1):
			f_write(f,0,1,TestCase.VarName+' := '+TestCase.MeshOutputs)

	f_write(f,0,1,'')

	for i in range(0,len(TestCases)):
		TestCase = TestCases[i]
		if (TestCase.name.find('update_h') != -1):
			# ToBeDeleted: Remove the if condition below after control files are updated
#			if (TestCase.GeoDeps.find('n-Cube') == -1):
			f_write(f,0,1,'$('+TestCase.VarName+') : '+TestCase.GeoDeps)


if __name__ == '__main__':
	"""Generate the MeshVariables file based on .ctrl files specified in the Makefile of the ROOT directory."""
	user = 'PZwan'

	Paths = Paths_class()
	Paths.set_paths(user)

	SubDirectories = sys.argv[1:]
	print(SubDirectories)

	TestCases = []

	for i in range(0,len(SubDirectories)):
		TestCase = TestCase_class(SubDirectories[i])
		TestCase.set_paths(Paths)
		TestCase.add_MeshTypes(Paths,'all')
		
		if(TestCase.name.find('update_h') != -1):
			print(TestCase.name)
			TestCase.get_geo_dependencies()
			TestCase.get_mesh_outputs()

			print(TestCase.GeoDeps,'\n')
			print(TestCase.MeshOutputs,'\n')

		TestCases.append(TestCase)

	# Output MeshVariables file
	write_MeshVariables(TestCases,Paths)


