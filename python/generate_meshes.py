import sys
import subprocess
import shlex
import re

'''
Purpose:
	Generate .msh files necessary for the code using the gmsh .geo files.

Comments:
	Paths to necessary executables are set in set_paths.
'''

### Classes ###
from meshfile_classes import Paths_class
from meshfile_classes import TestCase_class



### Functions ###
from support_functions import EXIT_TRACEBACK


def create_meshes(TestCase,Paths):
	for i in range(0,len(TestCase.MeshTypes)):
		MeshType = TestCase.MeshTypes[i]

		gmsh_args = ' ' + Paths.meshes + MeshType.InputName
		gmsh_args += ' -' + MeshType.dim
		gmsh_args = add_gmsh_setnumber(gmsh_args,MeshType,Paths)
		gmsh_args += ' -o ' + MeshType.OutputName

#		print("MON:",MeshType.OutputName)
		subprocess.call(shlex.split('mkdir -p ' + MeshType.OutputDir))
		subprocess.call(shlex.split(Paths.gmsh + gmsh_args))


def add_gmsh_setnumber(gmsh_args,MeshType,Paths):
	""" Set numbers for gmsh command line arguments based on values in Parameters.geo. """

	# Parameters which should always be present (using get_gmsh_number)

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


	# Other parameters

	# Extended domain
	gmsh_args += ' -setnumber Extended '
	if (MeshType.GeomSpecifier.find('Extended') == -1):
		gmsh_args += get_gmsh_number(gmsh_args,"Extension_Disabled",Paths)
	else:
		gmsh_args += get_gmsh_number(gmsh_args,"Extension_Enabled",Paths)

	gmsh_args += ' -setnumber Geom_AR '
	if (MeshType.GeomSpecifier.find('AR') != -1):
		AR_Num = re.search('\/AR_(.+?)\/',MeshType.GeomSpecifier)
		gmsh_args += get_gmsh_number(gmsh_args,"Geom_AR_"+AR_Num.group(1),Paths)
	else:
		gmsh_args += get_gmsh_number(gmsh_args,"Geom_AR_1",Paths)

	gmsh_args += ' -setnumber Geom_Adv '
	print('MeshType',MeshType.GeomSpecifier)
	if (MeshType.GeomSpecifier.find('YL') != -1):
		gmsh_args += get_gmsh_number(gmsh_args,"Geom_Adv_YL",Paths)
	elif (MeshType.GeomSpecifier.find('EllipticPipe') != -1):
		print('here?')
		gmsh_args += get_gmsh_number(gmsh_args,"Geom_Adv_Manufactured",Paths)
	else: # Default
		gmsh_args += get_gmsh_number(gmsh_args,"Geom_Adv_NONE",Paths)

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

	Paths = Paths_class()
	Paths.set_paths()

	CaseName = sys.argv[1]
	MeshName = sys.argv[2]

	if (len(sys.argv) > 3):
		print("Input arguments should be limitted to 'TestCase MeshFile'")
		EXIT_TRACEBACK()

	print('\n\n\nGenerating '+MeshName+'.\n\n')

	TestCase = TestCase_class(CaseName)

	TestCase.set_paths(Paths)
	TestCase.add_MeshTypes(Paths,MeshName)
	create_meshes(TestCase,Paths)
