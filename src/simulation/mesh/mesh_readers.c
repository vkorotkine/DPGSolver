/* {{{
This file is part of DPGSolver.

DPGSolver is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or any later version.

DPGSolver is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
Public License for more details.

You should have received a copy of the GNU General Public License along with DPGSolver.  If not, see
<http://www.gnu.org/licenses/>.
}}} */
///	\file

#include "mesh_readers.h"

#include <stdlib.h>
#include <string.h>

#include "macros.h"
#include "definitions_core.h"
#include "definitions_alloc.h"
#include "definitions_elements.h"

#include "multiarray.h"
#include "matrix.h"
#include "vector.h"

#include "mesh.h"
#include "file_processing.h"
#include "const_cast.h"

// Static function declarations ************************************************************************************* //

#define OUTPUT_MESH_DATA false ///< Flag for whether data should be output.
#define GMSH_N_TAGS 2          ///< Expected number of tags for elements in the gmsh file.

/// \brief Holds data relating to elements in the gmsh file.
struct Element_Data {
	ptrdiff_t n_elems; ///< The number of physical elements.

	struct Vector_i* elem_types; ///< Defined in \ref Mesh_Data.
	struct Matrix_i* elem_tags;  ///< Defined in \ref Mesh_Data.

	struct Multiarray_Vector_i* node_nums; ///< Defined in \ref Mesh_Data.
};

/// \brief Container for locally computed \ref Mesh_Connectivity members.
struct Mesh_Data_l {
	struct Matrix_d*     nodes;         ///< Defined in \ref Mesh_Data.
	struct Element_Data* elem_data;     ///< \ref Element_Data.
	struct Matrix_i*     periodic_corr; ///< Define in \ref Mesh_Data.
};

/// /brief Read data from a mesh in gmsh format.
static void mesh_reader_gmsh
	(const char*const mesh_name_full,     ///< The name of the mesh including the full path.
	 const int d,                         ///< The dimension.
	 struct Mesh_Data_l*const mesh_data_l ///< \ref Mesh_Data_l.
	);

/** \brief See return.
 *	\return The number of elements of each dimension.
 */
static struct Vector_i* count_elements_per_dim
	(const struct const_Vector_i*const elem_types ///< Defined in \ref Conn_info.
	);

// Interface functions ********************************************************************************************** //

struct Mesh_Data* constructor_Mesh_Data (const char*const mesh_name_full, const int d)
{
	struct Mesh_Data_l mesh_data_l;
	if (strstr(mesh_name_full,".msh"))
		mesh_reader_gmsh(mesh_name_full,d,&mesh_data_l);
	else
		EXIT_UNSUPPORTED;

	struct Mesh_Data* mesh_data = calloc(1,sizeof *mesh_data); // returned

	const_constructor_move_Matrix_d(&mesh_data->nodes,mesh_data_l.nodes);

	const_constructor_move_Vector_i(&mesh_data->elem_types,mesh_data_l.elem_data->elem_types);
	const_constructor_move_Matrix_i(&mesh_data->elem_tags,mesh_data_l.elem_data->elem_tags);
	const_constructor_move_Multiarray_Vector_i(&mesh_data->node_nums,mesh_data_l.elem_data->node_nums);

	if (mesh_data_l.periodic_corr)
		const_constructor_move_Matrix_i(&mesh_data->periodic_corr,mesh_data_l.periodic_corr);
	else
		*(struct const_Matrix_i**)&mesh_data->periodic_corr = NULL;

	struct Vector_i* elem_per_dim = count_elements_per_dim(mesh_data->elem_types); // keep
	const_constructor_move_Vector_i(&mesh_data->elem_per_dim,elem_per_dim);

	const ptrdiff_t ind_v = get_first_volume_index(mesh_data->elem_per_dim,d);

	const_cast_i(&mesh_data->d,d);
	const_cast_ptrdiff(&mesh_data->ind_v,ind_v);

	if (OUTPUT_MESH_DATA) {
		print_const_Vector_i(mesh_data->elem_per_dim);
		print_const_Matrix_d(mesh_data->nodes);
		print_const_Vector_i(mesh_data->elem_types);
		print_const_Matrix_i(mesh_data->elem_tags);
		print_const_Multiarray_Vector_i(mesh_data->node_nums);
		if (mesh_data->periodic_corr)
			print_const_Matrix_i(mesh_data->periodic_corr);
		EXIT_ERROR("Disable output to continue.\n");
	}

	free(mesh_data_l.elem_data);

	return mesh_data;
}

void destructor_Mesh_Data (struct Mesh_Data* mesh_data)
{
	destructor_const_Vector_i(mesh_data->elem_per_dim);
	destructor_const_Matrix_d(mesh_data->nodes);

	destructor_const_Vector_i(mesh_data->elem_types);
	destructor_const_Matrix_i(mesh_data->elem_tags);
	destructor_const_Multiarray_Vector_i(mesh_data->node_nums);

	if (mesh_data->periodic_corr)
		destructor_const_Matrix_i(mesh_data->periodic_corr);

	free(mesh_data);
}

// Static functions ************************************************************************************************* //

// Gmsh ************************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Read the nodes (xyz coordinates) from the mesh file.
 *	\return A \ref Matrix_T\* containing the \ref Mesh_Data::nodes.
 */
static struct Matrix_d* read_nodes
	(FILE* mesh_file, ///< The mesh file.
	 const int d      ///< The dimension.
	);

/** \brief Read the element data from the mesh file.
 *	\return A \ref Element_Data\* container holding the element related mesh data.
 */
static struct Element_Data* read_elements
	(FILE* mesh_file ///< The mesh file.
	);

/** \brief Reads periodic entity data.
 *	\return A \ref Matrix_T\* containing the \ref Mesh_Data::periodic_corr.
 *
 *	The use of `line_ptr` in the code below is required for passing the address of the `char[]` line to a `char**` as
 *	explained in [this SO answer][multidim_arrays].
 *
 *	[multidim_arrays]: https://stackoverflow.com/questions/1584100/converting-multidimensional-arrays-to-pointers-in-c
 */
static struct Matrix_i* read_periodic
	(FILE* mesh_file, ///< The mesh file.
	 const int d      ///< The dimension.
	);

static void mesh_reader_gmsh (const char*const mesh_name_full, const int d, struct Mesh_Data_l*const mesh_data_l)
{
	mesh_data_l->nodes         = NULL;
	mesh_data_l->elem_data     = NULL;
	mesh_data_l->periodic_corr = NULL;

	FILE* mesh_file = fopen_checked(mesh_name_full); // closed

	char line[STRLEN_MAX];
	while (fgets(line,sizeof(line),mesh_file)) {
		if (strstr(line,"$Nodes"))
			mesh_data_l->nodes = read_nodes(mesh_file,d); // keep

		if (strstr(line,"$Elements"))
			mesh_data_l->elem_data = read_elements(mesh_file); // keep

		if (strstr(line,"Periodic"))
			mesh_data_l->periodic_corr = read_periodic(mesh_file,d); // keep
	}
	fclose(mesh_file);
}

static struct Vector_i* count_elements_per_dim (const struct const_Vector_i*const elem_types)
{
	struct Vector_i* count = constructor_empty_Vector_i(DMAX+1); // returned
	set_to_zero_Vector_i(count);
	for (ptrdiff_t i = 0; i < elem_types->ext_0; i++) {
		const int elem_type = elem_types->data[i];

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
			EXIT_ERROR("Unsupported: %d\n",elem_type);
			break;
		}
	}

	return count;
}

// Level 1 ********************************************************************************************************** //

/** \brief Fill one row of the nodes \ref Matrix_T.
 *	Note that the first entry of the line is the node index and is discarded.
 */
static void fill_nodes
	(double*const node_row, ///< The current row.
	 char* line,            ///< The current line of the file.
	 const int d            ///< The dimension.
	);

/** \brief Constructor for \ref Element_Data\*.
 *	\return Standard.
 */
static struct Element_Data* constructor_Element_Data
	(const ptrdiff_t n_elems ///< The number of elements.
	);

/// \brief Fill one row of the members of \ref Element_Data.
static void fill_elements
	(const ptrdiff_t row,                 ///< The current row.
	 struct Element_Data*const elem_data, ///< The \ref Element_Data.
	 char* line                           ///< The current line of the file.
	);

/// \brief Skips the current periodic entity
static void skip_periodic_entity
	(FILE* file,         ///< The file.
	 char**const line,   ///< The pointer to the current line.
	 const int line_size ///< The size of the line array.
	);

static struct Matrix_d* read_nodes (FILE* mesh_file, const int d)
{
	char line[STRLEN_MAX];
	char* endptr = NULL;

	if (fgets(line,sizeof(line),mesh_file) != NULL) {};
	ptrdiff_t n_nodes = strtol(line,&endptr,10);
	struct Matrix_d* nodes = constructor_empty_Matrix_d('R',n_nodes,d);

	ptrdiff_t row = 0;
	while (fgets(line,sizeof(line),mesh_file) != NULL) {
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

	if (fgets(line,sizeof(line),mesh_file) != NULL) {};
	ptrdiff_t n_elems = strtol(line,&endptr,10);

	struct Element_Data* elem_data = constructor_Element_Data(n_elems);

	ptrdiff_t row = 0;
	while (fgets(line,sizeof(line),mesh_file) != NULL) {
		if (strstr(line,"$EndElements"))
			break;

		fill_elements(row,elem_data,line);

		if (row++ == n_elems)
			EXIT_UNSUPPORTED;
	}
	return elem_data;
}

static struct Matrix_i* read_periodic (FILE* mesh_file, const int d)
{
	char line[STRLEN_MAX];
	char* endptr = NULL;

	if (fgets(line,sizeof(line),mesh_file) != NULL) {};
	ptrdiff_t n_periodic_all = strtol(line,&endptr,10);

	// Skip over lower dimensional periodic entities if present
	ptrdiff_t n_periodic_low = 0;
	while (fgets(line,sizeof(line),mesh_file) != NULL) {
		if (strstr(line,"$EndPeriodic"))
			EXIT_UNSUPPORTED;

		ptrdiff_t dim_entity = strtol(line,&endptr,10);

		if (dim_entity < d-1) {
			char* line_ptr[1] = {line};
			skip_periodic_entity(mesh_file,line_ptr,sizeof(line));
			++n_periodic_low;
		} else {
			break;
		}
	}

	// Store periodic entity correspondence d-1 dimensional periodic entities.
	ptrdiff_t n_periodic = n_periodic_all-n_periodic_low;

	if (n_periodic == 0)
		EXIT_UNSUPPORTED;

	struct Matrix_i* periodic_corr = constructor_empty_Matrix_i('R',n_periodic,2);

	ptrdiff_t row = 0;
	do {
		if (strstr(line,"$EndPeriodic"))
			break;
		char* line_ptr[1] = {line};
		discard_line_values(line_ptr,1);
		read_line_values_i(line_ptr,periodic_corr->ext_1,get_row_Matrix_i(row,periodic_corr),false);
		skip_periodic_entity(mesh_file,line_ptr,sizeof(line));

		if (row++ == n_periodic)
			EXIT_UNSUPPORTED;
	} while (fgets(line,sizeof(line),mesh_file) != NULL);

	return periodic_corr;
}

// Level 2 ********************************************************************************************************** //

/** \brief Get the number of nodes specifying the geometry for the element of the given type.
 *	\return See brief.
 *
 *	The convention for the element type numbering is that of gmsh.
 */
static int get_n_nodes
	(const int elem_type ///< The element type.
	);

/// \brief Reorder the nodes such that they correspond to the ordering convention of this code.
static void reorder_nodes
	(const int elem_type, ///< Defined in \ref Mesh_Data.
	 struct Vector_i* node_nums   ///< Defined in \ref Mesh_Data.
	);

static void fill_nodes (double*const node_row, char* line, const int d)
{
	discard_line_values(&line,1);

	char* endptr = NULL;
	for (int dim = 0; dim < d; dim++) {
		node_row[dim] = strtod(line,&endptr);
		line = endptr;
	}
}

static struct Element_Data* constructor_Element_Data (const ptrdiff_t n_elems)
{
	struct Element_Data* elem_data = calloc(1,sizeof *elem_data); // returned

	elem_data->n_elems = n_elems;

	elem_data->elem_types = constructor_empty_Vector_i(n_elems);                    // keep
	elem_data->elem_tags  = constructor_empty_Matrix_i('R',n_elems,GMSH_N_TAGS);    // keep
	elem_data->node_nums  = constructor_empty_Multiarray_Vector_i(true,1,&n_elems); // keep

	return elem_data;
}

static void fill_elements (const ptrdiff_t row, struct Element_Data*const elem_data, char* line)
{
	ptrdiff_t n_tags;

	discard_line_values(&line,1);

	read_line_values_i(&line,1,&elem_data->elem_types->data[row],false);

	read_line_values_l(&line,1,&n_tags,false);
	if (n_tags != elem_data->elem_tags->ext_1)
		EXIT_UNSUPPORTED;

	read_line_values_i(&line,n_tags,get_row_Matrix_i(row,elem_data->elem_tags),false);

	int n_nodes = get_n_nodes(elem_data->elem_types->data[row]);
	resize_Vector_i(elem_data->node_nums->data[row],n_nodes);

	read_line_values_i(&line,n_nodes,elem_data->node_nums->data[row]->data,true);
	reorder_nodes(elem_data->elem_types->data[row],elem_data->node_nums->data[row]);
}

static void skip_periodic_entity (FILE* file, char**const line, const int line_size)
{
	char* endptr = NULL;

	if (fgets(*line,line_size,file) != NULL) {};
	ptrdiff_t n_skip = strtol(*line,&endptr,10);

	skip_lines_ptr(file,line,line_size,(int)n_skip);
}

// Level 3 ********************************************************************************************************** //

static int get_n_nodes (const int elem_type)
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

static void reorder_nodes (const int elem_type, struct Vector_i* node_nums)
{
	int gmsh_ordering[NVEMAX];

	switch (elem_type) {
		case POINT: case LINE: case TRI: case TET: case WEDGE:
			// Do nothing.
			return;
			break;
#if DIM > 1
		case QUAD: {
			const int gmsh_ordering_l[] = {0,1,3,2};
			memcpy(gmsh_ordering,gmsh_ordering_l,sizeof(gmsh_ordering_l));
			break;
		}
#endif
#if DIM > 2
		case HEX: {
			const int gmsh_ordering_l[] = {0,1,3,2,4,5,7,6};
			memcpy(gmsh_ordering,gmsh_ordering_l,sizeof(gmsh_ordering_l));
			break;
		} case PYR: {
			const int gmsh_ordering_l[] = {0,1,3,2,4};
			memcpy(gmsh_ordering,gmsh_ordering_l,sizeof(gmsh_ordering_l));
			break;
		}
#endif
		default:
			EXIT_UNSUPPORTED;
			break;
	}
	reorder_Vector_i(node_nums,gmsh_ordering);
}
