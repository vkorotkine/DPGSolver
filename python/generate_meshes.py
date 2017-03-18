import sys
import subprocess
import shlex

'''
Purpose:
	Generate .msh files necessary for the code using the gmsh .geo files.

Comments:
	Paths to necessary executables are set (in set_paths) according to the user and platform being used.
'''

### Classes ###
from meshfile_classes import Paths_class
from meshfile_classes import TestCase_class



### Functions ###
from support_functions import EXIT_TRACEBACK


def create_meshes(TestCase,Paths):
	for i in range(0,len(TestCase.MeshTypes)):
#	for i in range(0,1):
		MeshType = TestCase.MeshTypes[i]

		print(Paths.meshes)

		gmsh_args = ' ' + Paths.meshes + MeshType.InputName
		gmsh_args += ' -' + MeshType.dim
		gmsh_args = add_gmsh_setnumber(gmsh_args,MeshType,Paths)
		gmsh_args += ' -o ' + MeshType.OutputName
	
#		print(MeshType.OutputDir,'\n')
#		print(gmsh_args,'\n')
#		print(Paths.gmsh+gmsh_args,'\n')
#		EXIT_TRACEBACK()
		subprocess.call(shlex.split('mkdir -p ' + MeshType.OutputDir))
		subprocess.call(shlex.split(Paths.gmsh + gmsh_args))

		if (i == 1):
			print("Exiting\n\n")
			EXIT_TRACEBACK()

def add_gmsh_setnumber(gmsh_args,MeshType,Paths):
	""" Set numbers for gmsh command line arguments based on values in Parameters.geo. """
	# Need to treat case where one of these parameters is missing in the control file? ToBeDelted

	# MeshType
	gmsh_args += ' -setnumber MeshType '
	gmsh_args += get_gmsh_number(gmsh_args,MeshType.name,Paths)

	# PDEName
	gmsh_args += ' -setnumber PDEName '
	gmsh_args += get_gmsh_number(gmsh_args,MeshType.PDEName,Paths)

	# MeshCurving
	gmsh_args += ' -setnumber MeshCurving '
	gmsh_args += get_gmsh_number(gmsh_args,MeshType.MeshCurving,Paths)

	# MeshLevel
	gmsh_args += ' -setnumber MeshLevel ' + MeshType.MeshLevel

	return gmsh_args

def get_gmsh_number(gmsh_args,name,Paths):
	fName = Paths.meshes + 'Parameters.geo'

	Found = 0
	with open(fName) as f:
		for line in f:
			if (name.upper() in line):
				Found = 1
				return line.split()[2][:-1]

	print('Error: Did not find a value for '+name.upper()+' in '+fName+'\n')
	EXIT_BACKTRACE()


if __name__ == '__main__':
	""" Generate meshes based on the CaseList read from command line arguments. """

	user = 'PZwan'


	Paths = Paths_class()
	Paths.set_paths(user)

	CaseName = sys.argv[1]
	MeshName = sys.argv[2]

	if (len(sys.argv) > 3):
		print("Input arguments should be limitted to 'TestCase MeshFile'")
		EXIT_TRACEBACK()

	print('Generating '+MeshName+' for user '+user+'.\n')
	print(sys.argv)

	TestCase = TestCase_class(CaseName)

	TestCase.set_paths(Paths)
	TestCase.add_MeshTypes(Paths,MeshName)
	create_meshes(TestCase,Paths)

#	print(TestCase.name,TestCase.NMeshTypes)
#	for i in range(TestCase.NMeshTypes):
#		print(TestCase.MeshTypes[i].name,TestCase.MeshTypes[i].dim)


