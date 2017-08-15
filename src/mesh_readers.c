// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
///	\file

#include "mesh_readers.h"
#include "Multiarray.h"
#include "Matrix.h"
#include "Vector.h"

#include <stdlib.h>
#include <string.h>

#include "Macros.h"
#include "file_processing.h"
#include "allocators.h"
#include "const_cast.h"

#include "constants_gmsh.h"

static struct Mesh_Data* mesh_reader_gmsh (const char*const mesh_name_full, const unsigned int d);

struct Mesh_Data* mesh_reader (const char*const mesh_name_full, const unsigned int d)
{
	if (strstr(mesh_name_full,".msh"))
		return mesh_reader_gmsh(mesh_name_full,d);

	EXIT_UNSUPPORTED;
}

void destructor_Mesh_Data (struct Mesh_Data* mesh_data)
{
	destructor_Matrix_d((struct Matrix_d*)mesh_data->nodes);

	destructor_Vector_ui((struct Vector_ui*)mesh_data->elem_types);
	destructor_Matrix_ui((struct Matrix_ui*)mesh_data->elem_tags);
	destructor_Multiarray_Vector_ui((struct Multiarray_Vector_ui*)mesh_data->node_nums);

	destructor_Matrix_ui((struct Matrix_ui*)mesh_data->periodic_corr);

	free(mesh_data);
}


// Static functions ************************************************************************************************* //

// Gmsh ************************************************************************************************************* //

static struct Matrix_d*     read_nodes    (FILE* mesh_file, const unsigned int d);
static struct Element_Data* read_elements (FILE* mesh_file);
static struct Matrix_ui*    read_periodic (FILE* mesh_file, const unsigned int d);

static void                 fill_nodes           (double*const node_row, char* line, const unsigned int d);
static unsigned int         get_n_nodes          (const unsigned int elem_type);
static void                 fill_elements        (const size_t row, struct Element_Data*const elem_data, char* line);
static void                 reorder_nodes        (const unsigned int elem_type, struct Vector_ui* node_nums);
static void                 skip_periodic_entity (FILE* file, char**const line, const size_t line_size);

/// \brief Holds data relating to elements in the gmsh file.
struct Element_Data {
	size_t n_elems;

	struct Vector_ui* elem_types;
	struct Matrix_ui* elem_tags;

	struct Multiarray_Vector_ui* node_nums;
};

/// Expected number of tags for elements in the gmsh file.
#define GMSH_N_TAGS 2

/// \brief Constructor.
static struct Element_Data* constructor_Element_Data (const unsigned int n_elems)
{
	struct Element_Data* elem_data = malloc(1 * sizeof *elem_data); // returned

	elem_data->n_elems = n_elems;

	elem_data->elem_types = constructor_empty_Vector_ui(n_elems);                 // keep
	elem_data->elem_tags  = constructor_empty_Matrix_ui('R',n_elems,GMSH_N_TAGS); // keep
	elem_data->node_nums  = constructor_empty_Multiarray_Vector_ui(1,n_elems);    // keep

	return elem_data;
}

static struct Mesh_Data* constructor_Mesh_Data
	(struct Matrix_d* nodes, struct Element_Data* elem_data, struct Matrix_ui* periodic_corr)
{
	struct Mesh_Data* mesh_data = malloc(1 * sizeof *mesh_data);

	const_constructor_move_Matrix_d(&mesh_data->nodes,nodes);

	const_constructor_move_Vector_ui(&mesh_data->elem_types,elem_data->elem_types);
	const_constructor_move_Matrix_ui(&mesh_data->elem_tags,elem_data->elem_tags);
	const_constructor_move_Multiarray_Vector_ui(&mesh_data->node_nums,elem_data->node_nums);

	const_constructor_move_Matrix_ui(&mesh_data->periodic_corr,periodic_corr);

	free(elem_data);

	return mesh_data;
}

/** /brief Read data from a mesh in gmsh format.
 *
 *	\todo Move these comments to the appropriate function.
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
	struct Matrix_d*     nodes = NULL;
	struct Element_Data* elem_data = NULL;
	struct Matrix_ui*    periodic_corr = NULL;

	FILE* mesh_file = fopen_checked(mesh_name_full); // closed

	char line[STRLEN_MAX];
	while (fgets(line,sizeof(line),mesh_file)) {
		if (strstr(line,"$Nodes"))
			nodes = read_nodes(mesh_file,d); // keep

		if (strstr(line,"$Elements"))
			elem_data = read_elements(mesh_file); // keep

		if (strstr(line,"Periodic"))
			periodic_corr = read_periodic(mesh_file,d); // keep
	}
	fclose(mesh_file);

	struct Mesh_Data* mesh_data = constructor_Mesh_Data(nodes,elem_data,periodic_corr); // keep

	return mesh_data;
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

/// \brief Read the element data from the mesh file.
static struct Element_Data* read_elements (FILE* mesh_file)
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

/** \brief Reads periodic entity data.
 *
 *	The use of `line_ptr` in the code below is required for passing the address of the `char[]` line to a `char**` as
 *	explained in [this SO answer][multidim_arrays].
 *
 *	[multidim_arrays]: https://stackoverflow.com/questions/1584100/converting-multidimensional-arrays-to-pointers-in-c
 */
static struct Matrix_ui* read_periodic (FILE* mesh_file, const unsigned int d)
{
	char line[STRLEN_MAX];
	char* endptr = NULL;

	fgets(line,sizeof(line),mesh_file);
	size_t n_periodic_all = strtol(line,&endptr,10);

	// Skip over lower dimensional periodic entities if present
	size_t n_periodic_low = 0;
	while (fgets(line,sizeof(line),mesh_file)) {
		if (strstr(line,"$EndPeriodic"))
			EXIT_UNSUPPORTED;

		size_t dim_entity = strtol(line,&endptr,10);

		if (dim_entity < d-1) {
			char* line_ptr[1] = {line};
			skip_periodic_entity(mesh_file,line_ptr,sizeof(line));
			++n_periodic_low;
		} else {
			break;
		}
	}

	// Store periodic entity correspondence d-1 dimensional periodic entities.
	size_t n_periodic = n_periodic_all-n_periodic_low;

	if (n_periodic == 0)
		EXIT_UNSUPPORTED;

	struct Matrix_ui* periodic_corr = constructor_empty_Matrix_ui('R',n_periodic,2);

	size_t row = 0;
	do {
		if (strstr(line,"$EndPeriodic"))
			break;
		char* line_ptr[1] = {line};
		discard_line_values(line_ptr,1);
		read_line_values_ui(line_ptr,periodic_corr->extents[1],get_row_Matrix_ui(row,periodic_corr),false);
		skip_periodic_entity(mesh_file,line_ptr,sizeof(line));

		if (row++ == n_periodic)
			EXIT_UNSUPPORTED;
	} while (fgets(line,sizeof(line),mesh_file));

	return periodic_corr;
}

/** \brief Fill one row of the nodes \ref Matrix_d.
 *	Note that the first entry of the line is the node index and is discarded.
 */
static void fill_nodes (double*const node_row, char* line, const unsigned int d)
{
	discard_line_values(&line,1);

	char* endptr = NULL;
	for (unsigned int dim = 0; dim < d; dim++) {
		node_row[dim] = strtod(line,&endptr);
		line = endptr;
	}
}

/// \brief Fill one row of the members of \ref elem_data.
static void fill_elements (const size_t row, struct Element_Data*const elem_data, char* line)
{
	unsigned int n_tags;

	discard_line_values(&line,1);

	read_line_values_ui(&line,1,&elem_data->elem_types->data[row],false);

	read_line_values_ui(&line,1,&n_tags,false);
	if (n_tags != elem_data->elem_tags->extents[1])
		EXIT_UNSUPPORTED;

	read_line_values_ui(&line,n_tags,get_row_Matrix_ui(row,elem_data->elem_tags),false);

	unsigned int n_nodes = get_n_nodes(elem_data->elem_types->data[row]);
	reserve_Vector_ui(elem_data->node_nums->data[row],n_nodes);

	read_line_values_ui(&line,n_nodes,elem_data->node_nums->data[row]->data,true);
	reorder_nodes(elem_data->elem_types->data[row],elem_data->node_nums->data[row]);
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
}

/// \brief Reorder the nodes such that they correspond to the ordering convention of this code.
static void reorder_nodes (const unsigned int elem_type, struct Vector_ui* node_nums)
{
	const unsigned int n_nodes_max = 8;
	size_t gmsh_ordering[n_nodes_max];

	switch (elem_type) {
		case POINT: case LINE: case TRI: case TET: case WEDGE:
			// Do nothing.
			return;
			break;
		case QUAD: {
			const size_t gmsh_ordering_l[] = {0,1,3,2};
			memcpy(gmsh_ordering,gmsh_ordering_l,sizeof(gmsh_ordering_l));
			break;
		} case HEX: {
			const size_t gmsh_ordering_l[] = {0,1,3,2,4,5,7,6};
			memcpy(gmsh_ordering,gmsh_ordering_l,sizeof(gmsh_ordering_l));
			break;
		} case PYR: {
			const size_t gmsh_ordering_l[] = {0,1,3,2,4};
			memcpy(gmsh_ordering,gmsh_ordering_l,sizeof(gmsh_ordering_l));
			break;
		} default:
			EXIT_UNSUPPORTED;
			break;
	}

	reorder_Vector_ui(node_nums,gmsh_ordering);
}

/// \brief Skips the current periodic entity
static void skip_periodic_entity (FILE* file, char**const line, const size_t line_size)
{
	char* endptr = NULL;

	fgets(*line,line_size,file);
	size_t n_skip = strtol(*line,&endptr,10);

	skip_lines(file,line,line_size,n_skip);
}
