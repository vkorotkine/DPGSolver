import platform

from support_functions import EXIT_TRACEBACK

class Paths_class:
	def __init__(self):
		self.gmsh          = ''
		self.DPG_ROOT      = ''
		self.meshes        = ''
		self.cases         = ''
		self.control_files = ''

	def set_paths(self,user):
		""" Set paths based on the current user.  """
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

		self.cases         = self.DPG_ROOT+'cases/'
		self.control_files = self.cases+'control_files/'

class CtrlFile_class:
	""" Class storing relevant data for a control file group."""

	def __init__(self,name):
		self.name    = name
		self.VarName = ''

		self.MeshTypes  = []
		self.Path = ''
		self.GeoDeps = ''

	def add_MeshTypes(self,Paths):
		NStd = 6
		MeshTypes   = ['TRI','QUAD','TET','HEX','WEDGE','PYR']

		if (self.name.find('update_h') != -1):
			self.VarName = 'UPDATE_H'
			iRange = range(0,NStd)
		else:
			return # ToBeDeleted
			EXIT_TRACEBACK()

		for i in iRange:
			TypeCurrent = MeshType_class(MeshTypes[i])

			TypeCurrent.set_parameters(self,Paths)

			self.MeshTypes.append(TypeCurrent)

	def set_paths(self,Paths):
		if (self.name.find('update_h') != -1):
			self.Path = Paths.DPG_ROOT+self.name+'/Test_update_h_'

	def get_geo_dependencies(self):
		MeshTypes = self.MeshTypes
		for i in range(0,len(MeshTypes)):
			if (self.GeoDeps.find(MeshTypes[i].InputName) == -1):
				self.GeoDeps += MeshTypes[i].InputName + ' '


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

	def set_parameters(self,CtrlFile,Paths):
		fName = CtrlFile.Path+self.name+'.ctrl'

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
