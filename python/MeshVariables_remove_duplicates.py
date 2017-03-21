import re

### Classes ###
from meshfile_classes import Paths_class

### Functions ###
from support_functions import EXIT_TRACEBACK
from support_functions import f_write

if __name__ == '__main__':
	"""Remove duplicate msh file entries to avoid warnings in the Makefile."""
	user = 'PZwan'

	Paths = Paths_class()
	Paths.set_paths(user)

	fName = Paths.meshes+'MeshVariables'
	with open(fName, 'r') as f:
		data = f.read().split('\n')

	# Re-write to file removing duplicate .msh entries
	f = open(Paths.meshes+'MeshVariables','w')

	MeshList = set()
	for line in range(len(data)):
		MeshesPresent = re.search('(.+?).msh',data[line])
		if (not MeshesPresent):
			f_write(f,0,1,data[line])
		else:
			data_line = data[line].split()

			f_write(f,0,0,' '.join(data_line[0:2])+' ')

			MeshesCurrent = data_line[2:]
			MeshesOutput  = []
			for i in range(len(MeshesCurrent)):
				if (not MeshesCurrent[i] in MeshList):
					MeshesOutput.append(MeshesCurrent[i])

				MeshList.add(MeshesCurrent[i])

			f_write(f,0,1,' '.join(MeshesOutput))
