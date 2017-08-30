import sys
sys.path.insert(0,'../')

### Classes ###
#from meshfile_classes import Paths_class
#from meshfile_classes import TestCase_class

### Functions ###
from support_functions import f_write


class Paths:
	""" Container for paths. """
	def __init__ (self,project_bin_dir,project_src_dir):
		self.project_bin_dir = project_bin_dir ### ${PROJECT_BINARY_DIR} from CMake
		self.project_src_dir = project_src_dir ### ${PROJECT_SOURCE_DIR} from CMake

		### The root directory in which to generate the meshes.
		self.mesh_root = project_bin_dir+'/meshes'

		### The root directory in which to find the control files.
		self.ctrl_root = project_src_dir+'/input/cases/control_files'

class Mesh_Info:
	""" Container for mesh related information found in the control file. """
	def __init__ (self,ctrl_file,paths):
		self.ctrl_file_name = paths.ctrl_root+'/'+ctrl_file ###< Control file name (including full path).
		self.mesh_root      = paths.mesh_root+'/'           ###< Defined in \ref Paths.

		self.info = dict() ###< Dictionary which will hold information.

		### Information members to be found.
		self.info_members = ['pde_name',  ###< The pde name.
		                     'pde_spec',  ###< The pde specifier.
		                     'geom_name', ###< The geometry name.
		                     'geom_spec', ###< The geometry specifier.

		                     'dimension', ###< The dimension of the mesh.

		                     'mesh_generator', """< The name of the file used to generate the mesh (including relative
		                                            path) """
		                     'mesh_format',    ###< The format of the mesh
		                     'mesh_domain',    ###< The domain type.
		                     'mesh_type',      ###< The element types present in the mesh
		                     'mesh_level',     ###< The level of refinement of the mesh
		                    ]

	def read_data (self):
		with open(self.ctrl_file_name) as ctrl_file:
			for line in ctrl_file:
				line_split = line.split()
				if (line_split and line_split[0] in self.info_members):
					self.info[line_split[0]] = line_split[1]

	def assemble_mesh_name (self):
		mesh_file_name  = self.mesh_root
		mesh_file_name += self.info['geom_name']+'/'
		mesh_file_name += self.info['pde_name']+'/'
		if (self.info['pde_spec'] != "NONE"):
			mesh_file_name += self.info['pde_spec']+'/'
		if (self.info['geom_spec'] != "NONE"):
			mesh_file_name += self.info['geom_spec']+'/'

		supported_mesh_domains = ['straight','curved','parametric']
		if (self.info['mesh_domain'] not in supported_mesh_domains):
			EXIT_ERROR

		mesh_file_name += self.info['mesh_domain']+'_'
		mesh_file_name += self.info['dimension']+'d_'
		mesh_file_name += self.info['mesh_type']+'_'
		mesh_file_name += "ml"+self.info['mesh_level']
		mesh_file_name += ".msh"

		return mesh_file_name


if __name__ == '__main__':
	""" Generate the mesh dependency file based on the list of input .ctrl files.  """

	paths = Paths(sys.argv[1],sys.argv[2])

	mesh_depend_name = sys.argv[3]
	mesh_make_name   = sys.argv[4]
	ctrl_files       = sys.argv[5:]

	dependencies = dict()
	for ctrl_file in ctrl_files:
		mesh_info = Mesh_Info(ctrl_file,paths)
		mesh_info.read_data()

		mesh_file_name      = mesh_info.assemble_mesh_name()
		mesh_generator_name = paths.project_src_dir+'/input/meshes/'+mesh_info.info['mesh_generator']

		if (mesh_file_name not in dependencies.keys()):
			dependencies[mesh_file_name] = set()
		dependencies[mesh_file_name].add(mesh_generator_name)

	# Write to file
	meshes = set()
	depend_file = open(mesh_depend_name,'w')

	f_write(depend_file,0,2,"# Mesh file dependencies")
	for key, value in dependencies.items():
		write_item = key+" : "+(' '.join(str(s) for s in value))

		f_write(depend_file,0,1,write_item)

		meshes.add(key)

	f_write(depend_file,0,1,'')
	if ("Depend_Geo" in mesh_make_name):
		meshes_name = "MESHES_GEO"
	elif ("Depend_Python" in mesh_make_name):
		meshes_name = "MESHES_PYTHON"

	f_write(depend_file,0,1,meshes_name+" := "+' '.join(str(s) for s in meshes))

	depend_file.close()
