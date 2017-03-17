import sys


### Classes ###
from meshfile_classes import Paths_class
from meshfile_classes import CtrlFile_class


### Functions ###
if __name__ == '__main__':
	""" Generate the MeshVariables file based on .ctrl files specified in the Makefile of the ROOT directory."""
	user = 'PZwan'

	Paths = Paths_class()
	Paths.set_paths(user)

	SubDirectories = sys.argv[1:]
	print(SubDirectories)

	for i in range(0,len(SubDirectories)):
		CtrlFile = CtrlFile_class(SubDirectories[i])
		CtrlFile.set_paths(Paths)
		CtrlFile.add_MeshTypes(Paths)
		
		if(CtrlFile.name.find('update_h') != -1):
			print(CtrlFile.name)
			CtrlFile.get_geo_dependencies()
			print(CtrlFile.GeoDeps)
