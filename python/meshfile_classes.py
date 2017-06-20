import platform
import os

from support_functions import EXIT_TRACEBACK

class Paths_class:
	def __init__(self):
		self.gmsh          = ''
		self.DPG_ROOT      = ''
		self.meshes        = ''
		self.cases         = ''
		self.control_files = ''

	def set_paths(self):
		""" Set paths.  """
		self.DPG_ROOT = os.path.realpath(__file__)
		for i in range(0,2):
			self.DPG_ROOT = os.path.dirname(self.DPG_ROOT)
		self.DPG_ROOT += '/'

		self.gmsh          = 'gmsh' # Be sure to add the gmsh executable location to your $PATH
		self.meshes        = self.DPG_ROOT+'meshes/'
		self.cases         = self.DPG_ROOT+'cases/'
		self.control_files = self.cases+'control_files/'


def find_MeshType(N,MeshTypes,MeshTypesPrefix,MeshName,MeshCurving):
	Found = 0
	for i in range(0,N):
		if (MeshName.find(MeshCurving[i]+MeshTypes[i]) != -1):
			# Eliminate occurence of finding Curved in ToBeCurved ...
			if (MeshName.find('ToBeCurved') != -1 and MeshCurving[i].find('ToBeCurved') == -1):
				continue

			Found = 1
			MeshType    = [MeshTypes[i]]
			MeshCurving = [MeshCurving[i]]
			if (len(MeshTypesPrefix) == 1):
				MeshTypesPrefix = [MeshTypesPrefix[0]]
			else:
				MeshTypesPrefix = [MeshTypesPrefix[i]]

			break;

	if (Found == 0):
		print("Did not find the MeshType.\n")
		print(MeshName,"\n",MeshCurving,"\n",MeshTypes,"\n")
		EXIT_TRACEBACK()

	return [MeshType, MeshTypesPrefix, MeshCurving]

class TestCase_class:
	""" Class storing relevant data for a control file group."""

	def __init__(self,name):
		self.name    = name
		self.VarName = ''

		self.MeshTypes  = []
		self.Path = ''

		self.GeoDeps     = ''
		self.MeshOutputs = ''

	def add_MeshTypes(self,Paths,MeshName):
		MeshCurving     = []
		MeshTypes       = []
		MeshTypesPrefix = []
		if (self.name.find('update_h') != -1 or
		    self.name.find('L2_proj')  != -1):
			if (self.name.find('update_h') != -1):
				self.VarName = 'UPDATE_H'
			elif (self.name.find('L2_proj_p') != -1):
				self.VarName = 'L2_PROJ_P'
			elif (self.name.find('L2_proj_h') != -1):
				self.VarName = 'L2_PROJ_H'

			NTotal = 6
			MeshTypesPrefix = ['' for i in range(NTotal)]
			MeshCurving     = ['' for i in range(NTotal)]
			MeshTypes       = ['TRI','QUAD','TET','HEX','WEDGE','PYR']
			if   (MeshName.find('all')    != -1):
				iRange = range(0,NTotal)
			else:
				iRange = range(0,1)
				[MeshTypes,MeshTypesPrefix,MeshCurving] = find_MeshType(NTotal,MeshTypes,MeshTypesPrefix,MeshName,MeshCurving)
		elif (self.name.find('linearization') != -1):
			self.VarName = 'LINEARIZATION'

			NTotal = 3
			MeshTypesPrefix = ['' for i in range(NTotal)]
			MeshCurving     = ['ToBeCurved' for i in range(NTotal)]
			MeshTypes       = ['MIXED2D','MIXED3D_TP','MIXED3D_HW']
			if   (MeshName.find('all')    != -1):
				iRange = range(0,NTotal)
			else:
				iRange = range(0,1)
				[MeshTypes,MeshTypesPrefix,MeshCurving] = find_MeshType(NTotal,MeshTypes,MeshTypesPrefix,MeshName,MeshCurving)
		elif (self.name.find('Advection') != -1):
			if (self.name.lower().find('test') != -1):
				self.VarName = 'ADVECTION_TEST'

				MeshCurving.extend(('' for i in range (0,2)))
				MeshTypes.extend(('TRI','QUAD'))
				MeshTypesPrefix.extend(('Default_n-Cube_' for i in range(0,2)))

				NTotal = len(MeshTypes)
				if (MeshName.find('all') != -1):
					iRange = range(0,NTotal)
				else:
					iRange = range(0,1)
					[MeshTypes,MeshTypesPrefix,MeshCurving] = find_MeshType(NTotal,MeshTypes,MeshTypesPrefix,MeshName,MeshCurving)
			else:
				EXIT_TRACEBACK()
		elif (self.name.find('Poisson') != -1):
			if (self.name.lower().find('test') != -1):
				self.VarName = 'POISSON_TEST'

				MeshCurving.extend(('Curved'     for i in range (0,4)))
				MeshTypes.extend(('MIXED2D','TRI','QUAD','MIXED2D'))
				MeshTypesPrefix.extend(('n-Ball_HollowSection_'      for i in range(0,1)))
				MeshTypesPrefix.extend(('n-Ellipsoid_HollowSection_' for i in range(0,3)))

				NTotal = len(MeshTypes)
				if (MeshName.find('all') != -1):
					iRange = range(0,NTotal)
				else:
					iRange = range(0,1)
					[MeshTypes,MeshTypesPrefix,MeshCurving] = find_MeshType(NTotal,MeshTypes,MeshTypesPrefix,MeshName,MeshCurving)
			else:
				EXIT_TRACEBACK()
		elif (self.name.find('Euler') != -1):
			if (self.name.lower().find('test') != -1):
				self.VarName = 'EULER_TEST'

				MeshCurving.extend((''           for i in range (0,2)))
				MeshCurving.extend(('Curved'     for i in range (0,1)))
				MeshCurving.extend(('ToBeCurved' for i in range (0,6)))

				MeshTypes.extend(('TRI','QUAD'))
				MeshTypes.append(('MIXED2D'))
				MeshTypes.extend(('MIXED2D','MIXED3D_TP','MIXED3D_HW','TET','HEX','WEDGE'))

				MeshTypesPrefix.extend(('PeriodicVortex_'   for i in range(0,2)))
				MeshTypesPrefix.extend(('SupersonicVortex_' for i in range(0,7)))

				NTotal = len(MeshTypes)
				if   (MeshName.find('all')    != -1):
					iRange = range(0,NTotal)
				else:
					iRange = range(0,1)
					[MeshTypes,MeshTypesPrefix,MeshCurving] = find_MeshType(NTotal,MeshTypes,MeshTypesPrefix,MeshName,MeshCurving)
			else:
				print(self.name)
				EXIT_TRACEBACK()
		elif (self.name.find('NavierStokes') != -1):
			if (self.name.lower().find('test') != -1):
				self.VarName = 'NAVIERSTOKES_TEST'

				MeshCurving.extend((''           for i in range (0,1)))
				MeshCurving.extend(('ToBeCurved' for i in range (0,3)))

				MeshTypes.append(('QUAD'))
				MeshTypes.extend(('TRI','QUAD','MIXED2D'))

				MeshTypesPrefix.extend(('PlaneCouette_'  for i in range(0,1)))
				MeshTypesPrefix.extend(('TaylorCouette_' for i in range(0,3)))

				NTotal = len(MeshTypes)
				if (MeshName.find('all') != -1):
					iRange = range(0,NTotal)
				else:
					iRange = range(0,1)
					[MeshTypes,MeshTypesPrefix,MeshCurving] = find_MeshType(NTotal,MeshTypes,MeshTypesPrefix,MeshName,MeshCurving)
			else:
				print(self.name)
				EXIT_TRACEBACK()
		else:
			EXIT_TRACEBACK()

		for i in iRange:
			TypeCurrent = MeshType_class(MeshTypes[i],MeshTypesPrefix[i]+MeshCurving[i])

			TypeCurrent.set_parameters(self,Paths)

			self.MeshTypes.append(TypeCurrent)

	def set_paths(self,Paths):
		if (self.name.find('update_h') != -1):
			self.name = 'Test_update_h_'
			self.Path = Paths.control_files+'test/update_h/'
		elif (self.name.find('L2_proj_p') != -1):
			self.name = 'Test_L2_proj_p_'
			self.Path = Paths.control_files+'test/L2_proj_p/'
		elif (self.name.find('L2_proj_h') != -1):
			self.name = 'Test_L2_proj_h_'
			self.Path = Paths.control_files+'test/L2_proj_h/'
		elif (self.name.find('linearization') != -1):
			self.name = 'Test_linearization_'
			self.Path = Paths.control_files+'test/linearization/'
		elif (self.name.find('Advection') != -1):
			if (self.name.find('test') != -1):
				self.name = 'Test_Advection_'
				self.Path = Paths.control_files+'test/Advection/'
			else:
				print("name: Advection")
				EXIT_TRACEBACK()
		elif (self.name.find('Poisson') != -1):
			if (self.name.find('test') != -1):
				self.name = 'Test_Poisson_'
				self.Path = Paths.control_files+'test/Poisson/'
			else:
				self.name = 'Poisson'
				print("name:",self.name)
				EXIT_TRACEBACK()
		elif (self.name.find('Euler') != -1):
			if (self.name.find('test') != -1):
				self.name = 'Test_Euler_'
				self.Path = Paths.control_files+'test/Euler/'
			else:
				print("name:",self.name)
				EXIT_TRACEBACK()
		elif (self.name.find('NavierStokes') != -1):
			if (self.name.find('test') != -1):
				self.name = 'Test_NavierStokes_'
				self.Path = Paths.control_files+'test/NavierStokes/'
			else:
				print("name:",self.name)
				EXIT_TRACEBACK()
		else:
			print("name:",self.name)
			EXIT_TRACEBACK()

	def get_geo_dependencies(self):
		MeshTypes = self.MeshTypes
		for i in range(0,len(MeshTypes)):
			if (self.GeoDeps.find(MeshTypes[i].InputName) == -1):
				self.GeoDeps += MeshTypes[i].InputName + ' '

	def get_mesh_outputs(self):
		MeshTypes = self.MeshTypes
		for i in range(0,len(MeshTypes)):
			if (self.MeshOutputs.find(MeshTypes[i].InputName) == -1):
				self.MeshOutputs += MeshTypes[i].OutputName_from_meshesROOT + ' '


class MeshType_class:
	""" Both GeomSpecifier and PDESpecifier are included in the output path in order to allow for meshes having the same
	    geometry but different boundary conditions to be stored. This may result in redundant mesh storage.
	"""
	def __init__(self,name,prefix):
		self.name   = name
		self.prefix = prefix

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
		self.OutputName_from_meshesROOT = ''

	def set_parameters(self,TestCase,Paths):
		fName = TestCase.Path+TestCase.name+self.prefix+self.name+'.ctrl'

		with open(fName) as f:
			for line in f:
				if ('PDEName' in line):
					self.PDEName = line.split()[1]
				if ('PDESpecifier' in line):
					self.PDESpecifier = line.split()[1]
				if ('Geometry' in line):
					self.Geometry = line.split()[1]
					self.InputName = line.split()[2]
				if ('GeomSpecifier' in line):
					self.GeomSpecifier = line.split()[1]
				if ('MeshCurving' in line):
					self.MeshCurving = line.split()[1]
				if ('Dimension' in line):
					self.dim = line.split()[1]
				if ('MeshLevel' in line):
					self.MeshLevel = line.split()[1]

		self.OutputDir  = self.Geometry + '/' + self.PDEName + '/'
		if (self.PDESpecifier.find("NONE") == -1):
			self.OutputDir += self.PDESpecifier + '/'
		if (self.GeomSpecifier.find("NONE") == -1):
			self.OutputDir += self.GeomSpecifier + '/'
		self.OutputName = self.OutputDir + self.Geometry + self.dim + 'D_'
		self.OutputDir  = Paths.meshes + self.OutputDir

		if (self.MeshCurving.find('Straight') != 0):
			self.OutputName += self.MeshCurving
		self.OutputName += self.name + self.MeshLevel + 'x.msh'

		self.OutputName_from_meshesROOT = self.OutputName
		self.OutputName = Paths.meshes + self.OutputName
