// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
///	\file

#include "mesh_readers.h"

#include <stdlib.h>
#include <string.h>

#include "file_processing.h"
#include "allocators.h"
#include "containers.h"

#include "constants_gmsh.h"
#include "Macros.h"

static struct Mesh_Data* mesh_reader_gmsh (const char*const mesh_name_full, const unsigned int d);

struct Mesh_Data* mesh_reader (const char*const mesh_name_full, const unsigned int d)
{
	if (strstr(mesh_name_full,".msh"))
		return mesh_reader_gmsh(mesh_name_full,d);

	EXIT_UNSUPPORTED;
}


// Static functions ************************************************************************************************* //

// Gmsh ************************************************************************************************************* //

static void             fill_nodes (double*const node_row, const char* line, const unsigned int d);
static struct Matrix_d* read_nodes (FILE* mesh_file, const unsigned int d);

/// \brief Holds data relating to elements in the gmsh file.
struct Element_Data {
	unsigned int*  elem_types,
	            ** elem_tags;

	struct Vector_ui** node_nums;
};

/// Expected number of tags for elements in the gmsh file.
#define GMSH_N_TAGS 2

/// \brief Constructor.
static struct Element_Data* constructor_Element_Data (const unsigned int n_elems)
{
	struct Element_Data* elem_data = malloc(1 * sizeof *elem_data); // returned

	elem_data->elem_types = mallocator(UINT_T,1,n_elems);             // keep
	elem_data->elem_tags  = mallocator(UINT_T,2,n_elems,GMSH_N_TAGS); // keep
	elem_data->node_nums  = malloc(n_elems * sizeof *(elem_data->node_nums)); // keep

	return elem_data;
}

/** \brief Get the number of nodes specifying the geometry for the element of the given type.
 *
 *	The convention for the element type numbering is that of gmsh.
 */
static unsigned int get_n_nodes (const unsigned int elem_type)
{
	switch (elem_type) {
		case POINT: return 1; break;
		case LINE:  return 2; break;
		case TRI:   return 3; break;
		case QUAD:  return 4; break;
		case TET:   return 4; break;
		case HEX:   return 8; break;
		case WEDGE: return 6; break;
		case PYR:   return 5; break;
		default:
			EXIT_UNSUPPORTED;
			break;
	}

// Element classes
}

static void fill_elements (const size_t row, struct Element_Data*const elem_data, const char* line)
{
//	char* endptr = NULL;
	unsigned int n_tags;

	discard_line_values(&line,1);

	read_line_values_ui(&line,1,&elem_data->elem_types[row]);

	read_line_values_ui(&line,1,&n_tags);
	if (n_tags != GMSH_N_TAGS)
		EXIT_UNSUPPORTED;

	for (size_t n = 0; n < n_tags; n++)
		read_line_values_ui(&line,1,&elem_data->elem_tags[n][row]);

	unsigned int n_nodes = get_n_nodes(elem_data->elem_types[row]);
	elem_data->node_nums[row] = constructor_empty_Vector_ui(n_nodes); // keep

	read_line_values_ui(&line,n_nodes,elem_data->node_nums[row]->data);

print_Vector_ui(elem_data->node_nums[row]);
}

/// \brief Read the element data from the mesh file.
static struct Element_Data* read_elements(FILE* mesh_file)
{
	char line[STRLEN_MAX];
	char* endptr = NULL;

	fgets(line,sizeof(line),mesh_file);
	size_t n_elems = strtol(line,&endptr,10);

	struct Element_Data* elem_data = constructor_Element_Data(n_elems);

	size_t row = 0;
	while (fgets(line,sizeof(line),mesh_file)) {
		if (strstr(line,"$EndElements"))
			break;

		fill_elements(row,elem_data,line);

		if (row++ == n_elems)
			EXIT_UNSUPPORTED;
	}
	return elem_data;
}

/** /brief Read data from a mesh in gmsh format.
 *
 *	When dealing with **periodic** meshes, adherence to the numbering convention specified below is required for the
 *	correspondence of periodic entities to be identified.
 *
 *	Numbering (Geometry):
 *		- 0001-1000: Points (all).
 *		- 1001-2000: Lines x (constant yz).
 *		- 2001-3000: Lines y (constant xz).
 *		- 3001-4000: Lines z (constant xy).
 *		- 4001-5000: Surfaces xy (constant z).
 *		- 5001-6000: Surfaces xz (constant y).
 *		- 6001-7000: Surfaces yz (constant x).
 *		- 7001-8000: Volumes xyz.
 *
 *	Numbering (Physical):
 *		- \*0051 : Periodic (-x)
 *		- \*0052 : Periodic (+x)
 *		- \*0053 : Periodic (-y)
 *		- \*0054 : Periodic (+y)
 *		- \*0055 : Periodic (-z)
 *		- \*0056 : Periodic (+z)
 *
 *		\* : 1 (Straight), 2 (Curved). This is the same treatment as for all other BCs.
 *
 */
static struct Mesh_Data* mesh_reader_gmsh (const char*const mesh_name_full, const unsigned int d)
{
	struct Matrix_d* nodes = NULL;
	struct Element_Data* element_data = NULL;

	FILE* mesh_file = fopen_checked(mesh_name_full); // closed

	char line[STRLEN_MAX];
	while (fgets(line,sizeof(line),mesh_file)) {
		if (strstr(line,"$Nodes"))
			nodes = read_nodes(mesh_file,d); // tbd

		if (strstr(line,"$Elements"))
			element_data = read_elements(mesh_file); // tbd
//printf("%s",line);
	}

	fclose(mesh_file);

	struct Mesh_Data* mesh_data = malloc(1 * sizeof *mesh_data);

// Make casting functions for setting const lvalues.
	*(struct Matrix_d**)&mesh_data->nodes = nodes;
DO_NOTHING_P(element_data);
EXIT_UNSUPPORTED;

	return mesh_data;
}

/** \brief Fill one row of the nodes \ref Matrix_d.
 *	Note that the first entry of the line is the node index and is discarded.
 */
static void fill_nodes (double*const node_row, const char* line, const unsigned int d)
{
	discard_line_values(&line,1);

	char* endptr = NULL;
	for (unsigned int dim = 0; dim < d; dim++) {
		node_row[dim] = strtod(line,&endptr);
		line = endptr;
	}
}

/// \brief Read the nodes (xyz coordinates) from the mesh file.
static struct Matrix_d* read_nodes (FILE* mesh_file, const unsigned int d)
{
	char line[STRLEN_MAX];
	char* endptr = NULL;

	fgets(line,sizeof(line),mesh_file);
	size_t n_nodes = strtol(line,&endptr,10);
	struct Matrix_d* nodes = constructor_empty_Matrix_d('R',n_nodes,d);

	size_t row = 0;
	while (fgets(line,sizeof(line),mesh_file)) {
		if (strstr(line,"$EndNodes"))
			break;

		fill_nodes(get_row_Matrix_d(row,nodes),line,d);

		if (row++ == n_nodes)
			EXIT_UNSUPPORTED;
	}
	return nodes;
}


///\{ \name Definitions for the gmsh geometry numbering convention.
#define GMSH_XLINE_MIN  1001
#define GMSH_YLINE_MIN  2001
#define GMSH_ZLINE_MIN  3001
#define GMSH_XYFACE_MIN 4001
#define GMSH_XZFACE_MIN 5001
#define GMSH_YZFACE_MIN 6001
#define GMSH_XYZVOL_MIN 7001
///\}

