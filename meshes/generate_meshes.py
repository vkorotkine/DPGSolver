import sys
import subprocess
from subprocess import call
import os

def EXIT_TRACEBACK():
	print("\nError: Unsupported")
	traceback.print_tb()


class MeshType_class:
	def __init__(self,name):
		self.name = name
		self.dim  = 0

	def set_dim(self):
		name = self.name
		if (name == 'TRI' or name == 'QUAD'):
			self.dim = 2
		elif (name == 'TET' or name == 'HEX' or name == 'WEDGE' or name == 'PYR'):
			self.dim = 3
		else:
			EXIT_TRACEBACK()

class TestCase_class:
	""" Class storing relevant data for each test case. """

	def __init__(self,name):
		self.name = name

		self.CurvedType = []
		self.NMeshTypes = 0
		self.MeshTypes  = []

	def add_MeshType(self,MeshType):
		self.NMeshTypes += 1

		TypeCurrent = MeshType_class(MeshType)
		TypeCurrent.set_dim()
		self.MeshTypes.append(TypeCurrent)


def add_MeshTypes(TestCase):

	MeshTypeList = ['TRI','QUAD','TET','HEX','WEDGE','PYR']

	if (TestCase.name == 'update_h'):
		for i in range(0,6):
			TestCase.add_MeshType(MeshTypeList[i])
	else:
		EXIT_TRACEBACK()

def create_meshes(TestCase):
	os.environ['gmsh'] = "/home/philip/Desktop/research/programs/gmsh/gmsh-2.12.0-Linux/bin/gmsh"
#	for i in range(0,TestCase.NMeshTypes):
	for i in range(0,1):
		call(["ls","-a"])
		print('\n')
		subprocess.check_call(['gmsh','-version'])
#		subprocess.check_call(['gmsh','-version'],
#		                      env=dict(os.environ, GMSH_VAR="visible in this subprocess"))
#		subprocess.check_call(['sqsub', '-np', sys.argv[1], '/path/to/executable'],
#		                      env=dict(os.environ, SQSUB_VAR="visible in this subprocess"))
#		call(['gmsh','TRI.geo','-'+str(TestCase.MeshTypes[i].dim]))



if __name__ == '__main__':
	print("Generating Meshes.")

	NCases   = 2
	CaseList = ['update_h','tmp']

	for i in range(0,NCases):
		TestCase = TestCase_class(CaseList[i])

		add_MeshTypes(TestCase)
		create_meshes(TestCase)

		print(TestCase.name,TestCase.NMeshTypes)
		for i in range(TestCase.NMeshTypes):
			print(TestCase.MeshTypes[i].name,TestCase.MeshTypes[i].dim)

	# Loop over all meshes for the current case
	# Create appropriate directory
	# Create mesh, rename and move to directory

