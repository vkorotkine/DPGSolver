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
#include "Mesh.h"
#include "file_processing.h"
#include "allocators.h"
#include "const_cast.h"

#include "constants_core.h"
#include "constants_elements.h"

// Static function declarations ************************************************************************************* //

/// Expected number of tags for elements in the gmsh file.
#define GMSH_N_TAGS 2

/// /brief Read data from a mesh in gmsh format.
static struct Mesh_Data* mesh_reader_gmsh
	(const char*const mesh_name_full, ///< The name of the mesh including the full path.
	 const unsigned int d             ///< The dimension.
	);

// Interface functions ********************************************************************************************** //

struct Mesh_Data* mesh_reader (const char*const mesh_name_full, const unsigned int d)
{
	if (strstr(mesh_name_full,".msh"))
		return mesh_reader_gmsh(mesh_name_full,d);

	EXIT_UNSUPPORTED;
}

void destructor_Mesh_Data (struct Mesh_Data* mesh_data)
{
	destructor_Vector_ui((struct Vector_ui*)mesh_data->elem_per_dim);
	destructor_Matrix_d((struct Matrix_d*)mesh_data->nodes);

	destructor_Vector_ui((struct Vector_ui*)mesh_data->elem_types);
	destructor_Matrix_ui((struct Matrix_ui*)mesh_data->elem_tags);
	destructor_Multiarray_Vector_ui((struct Multiarray_Vector_ui*)mesh_data->node_nums);

	if (mesh_data->periodic_corr)
		destructor_Matrix_ui((struct Matrix_ui*)mesh_data->periodic_corr);

	free(mesh_data);
}

// Static functions ************************************************************************************************* //

// Gmsh ************************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Holds data relating to elements in the gmsh file.
struct Element_Data {
	size_t n_elems;

	struct Vector_ui* elem_types; ///< Defined in \ref Mesh_Data.
	struct Matrix_ui* elem_tags;  ///< Defined in \ref Mesh_Data.

	struct Multiarray_Vector_ui* node_nums; ///< Defined in \ref Mesh_Data.
};

/** \brief Read the nodes (xyz coordinates) from the mesh file.
 *	\return A \ref Matrix_d\* containing the \ref Mesh_Data::nodes.
 */
static struct Matrix_d* read_nodes
	(FILE* mesh_file,     ///< The mesh file.
	 const unsigned int d ///< The dimension.
	);

/** \brief Read the element data from the mesh file.
 *	\return A \ref Element_Data\* container holding the element related mesh data.
 */
static struct Element_Data* read_elements
	(FILE* mesh_file ///< The mesh file.
	);

/** \brief Reads periodic entity data.
 *	\return A \ref Matrix_ui\* containing the \ref Mesh_Data::periodic_corr.
 *
 *	The use of `line_ptr` in the code below is required for passing the address of the `char[]` line to a `char**` as
 *	explained in [this SO answer][multidim_arrays].
 *
 *	[multidim_arrays]: https://stackoverflow.com/questions/1584100/converting-multidimensional-arrays-to-pointers-in-c
 */
static struct Matrix_ui* read_periodic
	(FILE* mesh_file,     ///< The mesh file.
	 const unsigned int d ///< The dimension.
	);

/** \brief Constructor for the \ref Mesh_Data.
 *	\return Standard. */
static struct Mesh_Data* constructor_Mesh_Data
	(struct Matrix_d* nodes,         ///< Defined in \ref Mesh_Data.
	 struct Element_Data* elem_data, ///< \ref Element_Data.
	 struct Matrix_ui* periodic_corr ///< Defined in \ref Mesh_Data.
	);

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

// Level 1 ********************************************************************************************************** //

/** \brief Fill one row of the nodes \ref Matrix_d.
 *	Note that the first entry of the line is the node index and is discarded.
 */
static void fill_nodes
	(double*const node_row, ///< The current row.
	 char* line,            ///< The current line of the file.
	 const unsigned int d   ///< The dimension.
	);

/** \brief Constructor for \ref Element_Data\*.
 *	\return Standard.
 */
static struct Element_Data* constructor_Element_Data
	(const unsigned int n_elems ///< The number of elements.
	);

/// \brief Fill one row of the members of \ref elem_data.
static void fill_elements
	(const size_t row,                    ///< The current row.
	 struct Element_Data*const elem_data, ///< The \ref Element_Data.
	 char* line                           ///< The current line of the file.
	);

/// \brief Skips the current periodic entity
static void skip_periodic_entity
	(FILE* file,            ///< The file.
	 char**const line,      ///< The pointer to the current line.
	 const size_t line_size ///< The size of the line array.
	);

/** \brief See return.
 *	\return The number of elements of each dimension.
 */
static struct Vector_ui* count_elements_per_dim
	(const struct const_Vector_ui*const elem_types ///< Defined in \ref Conn_info.
	);

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

static struct Mesh_Data* constructor_Mesh_Data
	(struct Matrix_d* nodes, struct Element_Data* elem_data, struct Matrix_ui* periodic_corr)
{
	struct Mesh_Data* mesh_data = malloc(sizeof *mesh_data); // returned

	const_constructor_move_Matrix_d(&mesh_data->nodes,nodes);

	const_constructor_move_Vector_ui(&mesh_data->elem_types,elem_data->elem_types);
	const_constructor_move_Matrix_ui(&mesh_data->elem_tags,elem_data->elem_tags);
	const_constructor_move_Multiarray_Vector_ui(&mesh_data->node_nums,elem_data->node_nums);

	if (periodic_corr)
		const_constructor_move_Matrix_ui(&mesh_data->periodic_corr,periodic_corr);
	else
		*(struct const_Matrix_ui**)&mesh_data->periodic_corr = NULL;

	struct Vector_ui* elem_per_dim = count_elements_per_dim(mesh_data->elem_types); // keep
	const_constructor_move_Vector_ui(&mesh_data->elem_per_dim,elem_per_dim);

	const size_t d = nodes->extents[1];
	const unsigned int ind_v = get_first_volume_index(mesh_data->elem_per_dim,d);

	const_cast_ui(&mesh_data->d,d);
	const_cast_ui(&mesh_data->ind_v,ind_v);


	free(elem_data);

	return mesh_data;
}

// Level 2 ********************************************************************************************************** //

/** \brief Get the number of nodes specifying the geometry for the element of the given type.
 *	\return See brief.
 *
 *	The convention for the element type numbering is that of gmsh.
 */
static unsigned int get_n_nodes
	(const unsigned int elem_type ///< The element type.
	);

/// \brief Reorder the nodes such that they correspond to the ordering convention of this code.
static void reorder_nodes
	(const unsigned int elem_type, ///< Defined in \ref Mesh_Data.
	 struct Vector_ui* node_nums   ///< Defined in \ref Mesh_Data.
	);

static void fill_nodes (double*const node_row, char* line, const unsigned int d)
{
	discard_line_values(&line,1);

	char* endptr = NULL;
	for (unsigned int dim = 0; dim < d; dim++) {
		node_row[dim] = strtod(line,&endptr);
		line = endptr;
	}
}

static struct Element_Data* constructor_Element_Data (const unsigned int n_elems)
{
	struct Element_Data* elem_data = malloc(1 * sizeof *elem_data); // returned

	elem_data->n_elems = n_elems;

	elem_data->elem_types = constructor_empty_Vector_ui(n_elems);                 // keep
	elem_data->elem_tags  = constructor_empty_Matrix_ui('R',n_elems,GMSH_N_TAGS); // keep
	elem_data->node_nums  = constructor_empty_Multiarray_Vector_ui(1,n_elems);    // keep

	return elem_data;
}

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

static void skip_periodic_entity (FILE* file, char**const line, const size_t line_size)
{
	char* endptr = NULL;

	fgets(*line,line_size,file);
	size_t n_skip = strtol(*line,&endptr,10);

	skip_lines(file,line,line_size,n_skip);
}

static struct Vector_ui* count_elements_per_dim (const struct const_Vector_ui*const elem_types)
{
	struct Vector_ui* count = constructor_empty_Vector_ui(DMAX+1); // returned
	set_to_zero_Vector_ui(count);
	for (size_t i = 0; i < elem_types->extents[0]; i++) {
		const unsigned int elem_type = elem_types->data[i];

		switch (elem_type) {
		case POINT:
			count->data[0]++;
			break;
		case LINE:
			count->data[1]++;
			break;
		case TRI: case QUAD:
			count->data[2]++;
			break;
		case TET: case HEX: case WEDGE: case PYR:
			count->data[3]++;
			break;
		default:
			EXIT_UNSUPPORTED;
			break;
		}
	}

	return count;
}

// Level 3 ********************************************************************************************************** //

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

static void reorder_nodes (const unsigned int elem_type, struct Vector_ui* node_nums)
{
	const unsigned int n_nodes_max = 8;
	unsigned int gmsh_ordering[n_nodes_max];

	switch (elem_type) {
		case POINT: case LINE: case TRI: case TET: case WEDGE:
			// Do nothing.
			return;
			break;
		case QUAD: {
			const unsigned int gmsh_ordering_l[] = {0,1,3,2};
			memcpy(gmsh_ordering,gmsh_ordering_l,sizeof(gmsh_ordering_l));
			break;
		} case HEX: {
			const unsigned int gmsh_ordering_l[] = {0,1,3,2,4,5,7,6};
			memcpy(gmsh_ordering,gmsh_ordering_l,sizeof(gmsh_ordering_l));
			break;
		} case PYR: {
			const unsigned int gmsh_ordering_l[] = {0,1,3,2,4};
			memcpy(gmsh_ordering,gmsh_ordering_l,sizeof(gmsh_ordering_l));
			break;
		} default:
			EXIT_UNSUPPORTED;
			break;
	}
	reorder_Vector_ui(node_nums,gmsh_ordering);
}
