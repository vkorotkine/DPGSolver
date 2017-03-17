'''
	*** IMPORTANT ***   *** IMPORTANT ***   *** IMPORTANT ***   

	This script must be run using the Makefile such that the .pyc file is updated whenever this .py file is modified.
	Otherwise, changes made here will not affect how meshes are generated.

	*** IMPORTANT ***   *** IMPORTANT ***   *** IMPORTANT ***   
'''

import sys
import traceback
import subprocess
import shlex
import os
import platform

'''
Purpose:
	Generate .msh files necessary for the code using the gmsh .geo files.

Comments:
	Paths to necessary executables are set (in set_paths) according to the user and platform being used.
'''

### Special Functions ###
def EXIT_TRACEBACK():
	print('\nError: Unsupported')
	traceback.print_stack()
	sys.exit()


### Classes ###
class TestCase_class:
	""" Class storing relevant data for each test case. """

	def __init__(self,name):
		self.name = name

		self.MeshTypes  = []

		self.ctrlFilePath = ''

	def add_MeshTypes(self,Paths):
		NStd = 6
		MeshTypes   = ['TRI','QUAD','TET','HEX','WEDGE','PYR']

		if (self.name == 'update_h'):
			iRange = range(0,NStd)
		else:
			EXIT_TRACEBACK()

		for i in iRange:
			TypeCurrent = MeshType_class(MeshTypes[i])

			TypeCurrent.set_parameters(self,Paths)

			self.MeshTypes.append(TypeCurrent)

	def set_paths(self,Paths):
		if (self.name == 'update_h'):
			self.ctrlFilePath = Paths.cases+'test/'+self.name+'/Test_update_h_'

class MeshType_class:
	def __init__(self,name):
		self.name = name

		self.PDEName       = ''
		self.PDESpecifier  = ''
		self.Geometry      = ''
		self.GeomSpecifier = ''
		self.dim           = ''
		self.MeshLevel     = ''
		self.MeshCurving   = ''

		self.InputName  = ''
		self.OutputDir  = ''
		self.OutputName = ''

	def set_parameters(self,TestCase,Paths):
		fName = TestCase.ctrlFilePath+self.name+'.ctrl'

		with open(fName) as f:
			for line in f:
				if ('PDEName' in line):
					self.PDEName = line.split()[1]
				if ('PDESpecifier' in line):
					self.PDESpecifier = line.split()[1]
				if ('Geometry' in line):
					self.Geometry = line.split()[1]
				if ('GeomSpecifier' in line):
					self.GeomSpecifier = line.split()[1]
				if ('MeshCurving' in line):
					self.MeshCurving = line.split()[1]
				if ('Dimension' in line):
					self.dim = line.split()[1]
				if ('MeshLevel' in line):
					self.MeshLevel = line.split()[1]

		self.InputName  = self.Geometry + '/' + self.Geometry + self.dim + 'D.geo'
		self.OutputDir  = self.Geometry + '/' + self.PDEName + '/' \
		                + self.PDESpecifier + '/' + self.GeomSpecifier + '/'
		self.OutputDir  = Paths.meshes + self.OutputDir
		self.OutputName = self.OutputDir + self.Geometry + self.dim + 'D_'
		
		if (self.MeshCurving.find('Straight') != 0):
			self.OutputName += self.MeshCurving
		self.OutputName += self.name + self.MeshLevel + 'x.msh'


#		print('\nfName, OutputName:')
#		print(fName)
#		print(self.OutputDir,self.OutputName)
#		EXIT_TRACEBACK()

class Paths_class:
	def __init__(self):
		self.gmsh     = ''
		self.DPG_ROOT = ''
		self.cases    = ''
		self.meshes   = ''

	def set_paths(self,user):
		""" Set paths based on the current user for:
		    	gmsh,
		"""
		if (user == 'PZwan'):
			if any('Darwin' in string for string in platform.uname()):
				self.gmsh     = '/Applications/Gmsh.app/Contents/MacOS/gmsh'
				self.DPG_ROOT = '/Users/philip/Desktop/DPGSolver/'
			else:
				# Ubuntu (home)
				self.gmsh = '/home/philip/Desktop/research/programs/gmsh/gmsh-2.12.0-Linux/bin/gmsh'
				EXIT_TRACEBACK()
			self.meshes = self.DPG_ROOT+'meshes/'
		else:
			print('Add an option for yourself as a user.')
			EXIT_TRACEBACK()

		self.cases = self.DPG_ROOT+'cases/'


### Functions ###


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

	CaseList = sys.argv[1:]

	print('Generating Meshes for user '+user+str(CaseList)+'.\n')

	for i in range(0,len(CaseList)):
		TestCase = TestCase_class(CaseList[i])

		TestCase.set_paths(Paths)
		TestCase.add_MeshTypes(Paths)
		create_meshes(TestCase,Paths)

		print(TestCase.name,TestCase.NMeshTypes)
		for i in range(TestCase.NMeshTypes):
			print(TestCase.MeshTypes[i].name,TestCase.MeshTypes[i].dim)
