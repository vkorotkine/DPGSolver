DPG_ROOT    = '../../'
PYTHON_ROOT = DPG_ROOT + 'python/'
MESHES_ROOT = DPG_ROOT + 'meshes/'

import sys

sys.path.insert(0,PYTHON_ROOT)
from support_functions import f_write
from support_functions import EXIT_TRACEBACK

class TestCase_c:
	"""
	Defines test case related information.
	"""

	def __init__(self,name):
		self.name = name

		self.VarName     = ''
		self.MeshOutputs = ''

	def set_variables(self):
		if (self.name == 'Peterson'):
			self.VarName = 'ADVECTION_TEST_PETERSON'

			MESHES_SPECIFIC = 'n-Cube/Advection/Steady/Peterson/YL/'

			MeshOutputs = ''

			NML = 6 # Should match NML used in generate_peterson_meshes.py
			for ML in range(0,NML):
				MeshOutputs += MESHES_SPECIFIC + 'n-Cube2D_TRI' + str(ML) + 'x.msh '

			self.MeshOutputs = MeshOutputs

		else:
			EXIT_TRACEBACK()

def write_MeshVariables(TestCases):
	"""
	Write MeshVariables_python file.
	"""

	f = open(MESHES_ROOT + 'MeshVariables_python','w')

	f_write(f,0,2,'# Variables and dependencies needed for DPG_ROOT/meshes/Makefile for python generated meshes')

	for i in range(0,len(TestCases)):
		TestCase = TestCases[i]
		f_write(f,0,2,TestCase.VarName+' := '+TestCase.MeshOutputs)

	f.close()


if __name__ == '__main__':
	"""
	Generate MeshVariables_python file corresponding to meshes which are generated using python (as opposed to gmsh).
	This currently includes meshes for the following cases:
		Advection_Peterson
	"""

	TestCases = []
	TestCases.append(TestCase_c('Peterson'))

	for TestCase in TestCases:
		TestCase.set_variables()

	write_MeshVariables(TestCases)
