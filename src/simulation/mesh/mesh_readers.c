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
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "mpi.h"

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
#define PARTITION_OWN 3      ///< Index at which the partitions the element belongs to starts in the Gmsh tags
#define MAX_PHYS_ID 1 ///< Maximum number of physical IDs associated with an entity.

/// \brief Holds data relating to elements in the gmsh file.
struct Element_Data {
	ptrdiff_t n_elems; ///< The number of physical elements.

	ptrdiff_t ind_g; ///< Defined in \ref Mesh_Data
	struct Vector_i* elem_global_id;  ///< Defined in \ref Mesh_Data
	struct Vector_i* elem_owner_part; ///< Defined in \ref Mesh_Data
	struct Multiarray_Vector_i* elem_ghost_part; ///< Defined in \ref Mesh_Data

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
static void mesh_reader_gmsh_parallel
	(const char*const mesh_name_full,     ///< The name of the mesh including the full path.
	 const int d,                         ///< The dimension.
	 struct Mesh_Data_l*const mesh_data_l ///< \ref Mesh_Data_l.
	);

/// /brief Read data from a mesh in gmsh format.
static void mesh_reader_gmsh
	(const char*const mesh_name_full,     ///< The name of the mesh including the full path.
	 const int d,                         ///< The dimension.
	 struct Mesh_Data_l*const mesh_data_l ///< \ref Mesh_Data_l.
	);

static void mesh_reader_gmsh4
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
	printf("start of mesh reader");
	if(0==1) mesh_reader_gmsh(mesh_name_full,d,&mesh_data_l);
	if(0==1) mesh_reader_gmsh_parallel(mesh_name_full,d,&mesh_data_l);
	if(0==1) mesh_reader_gmsh4(mesh_name_full,d,&mesh_data_l);
	if (strstr(mesh_name_full,".msh"))
		//mesh_reader_gmsh(mesh_name_full,d,&mesh_data_l);
		//mesh_reader_gmsh_parallel(mesh_name_full,d,&mesh_data_l);
		mesh_reader_gmsh4(mesh_name_full,d,&mesh_data_l);
	else
		EXIT_ERROR("Unsupported (mesh_name: %s)\n",mesh_name_full);
	printf("End of mesh reader");

	struct Mesh_Data* mesh_data = calloc(1,sizeof *mesh_data); // returned

	const_constructor_move_Vector_i(&mesh_data->elem_global_id,mesh_data_l.elem_data->elem_global_id);
	const_constructor_move_Vector_i(&mesh_data->elem_owner_part,mesh_data_l.elem_data->elem_owner_part);
	const_constructor_move_Multiarray_Vector_i(&mesh_data->elem_ghost_part,mesh_data_l.elem_data->elem_ghost_part);

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

	destructor_const_Vector_i(mesh_data->elem_global_id);
	destructor_const_Vector_i(mesh_data->elem_owner_part);
	destructor_const_Multiarray_Vector_i(mesh_data->elem_ghost_part);

	destructor_const_Matrix_d(mesh_data->nodes);

	destructor_const_Vector_i(mesh_data->elem_types);
	destructor_const_Matrix_i(mesh_data->elem_tags);
	destructor_const_Multiarray_Vector_i(mesh_data->node_nums);

	if (mesh_data->periodic_corr)
		destructor_const_Matrix_i(mesh_data->periodic_corr);

	free(mesh_data);
}

void reorder_nodes_gmsh (const int elem_type, struct Vector_i* node_nums)
{
	int gmsh_ordering[NVEMAX];

	switch (elem_type) {
		case POINT: case LINE: case TRI: case TET: case WEDGE:
			return; // Do nothing.
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

int get_n_nodes (const int elem_type)
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
 *  When reading a value, the pointer moves down the char* array.
 *  line_ptr[0] = line is used to reset the point to the beginning of the line.
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
	(const ptrdiff_t n_elems, ///< The number of elements.
	 const ptrdiff_t n_tags   ///< Maximum number of tags
	);

/// \brief Fill one row of the members of \ref Element_Data.
static bool fill_elements
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

	fgets_checked(line,sizeof(line),mesh_file);
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

	fgets_checked(line,sizeof(line),mesh_file);
	ptrdiff_t n_elems = strtol(line,&endptr,10);

	struct Element_Data* elem_data = constructor_Element_Data(n_elems,2);

	ptrdiff_t row = 0;
	bool filled;
	while (fgets(line,sizeof(line),mesh_file) != NULL) {
		if (strstr(line,"$EndElements"))
			break;

		filled = fill_elements(row,elem_data,line);
		if (filled) row++;
		if (row > n_elems) EXIT_UNSUPPORTED;
	}
	return elem_data;
}

static struct Matrix_i* read_periodic (FILE* mesh_file, const int d)
{
	char line[STRLEN_MAX];
	char* endptr = NULL;

	fgets_checked(line,sizeof(line),mesh_file);
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

static void fill_nodes (double*const node_row, char* line, const int d)
{
	discard_line_values(&line,1);

	char* endptr = NULL;
	for (int dim = 0; dim < d; dim++) {
		node_row[dim] = strtod(line,&endptr);
		line = endptr;
	}
}

static struct Element_Data* constructor_Element_Data (const ptrdiff_t n_elems, const ptrdiff_t n_tags)
{
	struct Element_Data* elem_data = calloc(1,sizeof *elem_data); // returned

	elem_data->n_elems = n_elems;

	elem_data->elem_global_id = constructor_empty_Vector_i(n_elems); // keep
	elem_data->elem_owner_part = constructor_empty_Vector_i(n_elems); // keep
	elem_data->elem_ghost_part = constructor_empty_Multiarray_Vector_i(true,1,&n_elems);

	elem_data->elem_types = constructor_empty_Vector_i(n_elems);                    // keep
	elem_data->elem_tags  = constructor_empty_Matrix_i('R',n_elems,n_tags);// keep
	elem_data->node_nums  = constructor_empty_Multiarray_Vector_i(true,1,&n_elems); // keep

	return elem_data;
}

static bool fill_elements (const ptrdiff_t row, struct Element_Data*const elem_data, char* line)
{
	int elem_type;
	ptrdiff_t n_tags;
	int tags[GMSH_N_TAGS_MAX];

	discard_line_values(&line,1); 

	read_line_values_i(&line,1,&elem_type,false); // Element type

	read_line_values_l(&line,1,&n_tags,false); // Number of tags
//	if (n_tags != GMSH_N_TAGS)
//		EXIT_UNSUPPORTED;
//	if (n_tags != elem_data->elem_tags->ext_1)
//		EXIT_UNSUPPORTED;

	read_line_values_i(&line,n_tags,tags,false); // Tags

//	int mpi_rank;
//    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
//	bool skip_element = true;
//	for (ptrdiff_t i_tag=PARTITION_OWN; i_tag<n_tags; ++i_tag) {
//		ptrdiff_t i_part = tags[i_tag]; // i_part indexing is 1-based
//		if ((i_part-1)==mpi_rank || (-i_part-1)==mpi_rank) skip_element = false;
//	}
//	if (skip_element) return false;

	elem_data->elem_types->data[row] = elem_type;
	int* elem_tags = get_row_Matrix_i(row,elem_data->elem_tags);
	for (ptrdiff_t i_tag=0; i_tag<n_tags; ++i_tag) {
		elem_tags[i_tag] = tags[i_tag];
	}

	int n_nodes = get_n_nodes(elem_data->elem_types->data[row]);
	resize_Vector_i(elem_data->node_nums->data[row],n_nodes);

	read_line_values_i(&line,n_nodes,elem_data->node_nums->data[row]->data,true);
	reorder_nodes_gmsh(elem_data->elem_types->data[row],elem_data->node_nums->data[row]);

	return true;
}

static void skip_periodic_entity (FILE* file, char**const line, const int line_size)
{
	char* endptr = NULL;

	fgets_checked(*line,sizeof(*line),file);
	ptrdiff_t n_skip = strtol(*line,&endptr,10);

	skip_lines_ptr(file,line,line_size,(int)n_skip);
}

static void mesh_reader_gmsh_parallel (const char*const mesh_name_full, const int d, struct Mesh_Data_l*const mesh_data_l)
{
	int mpi_size, mpi_rank;
	MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
	MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);

	mesh_data_l->nodes         = NULL;
	mesh_data_l->elem_data     = NULL;
	mesh_data_l->periodic_corr = NULL;

	printf("Mesh name: %s\n", mesh_name_full);
	FILE* mesh_file = fopen_checked(mesh_name_full); // closed

	char line[STRLEN_MAX];
	char* line_ptr[1] = {line}; // Need to pass point to next functions
	char* endptr = NULL;

	// Go through file and count the number of partitions, elements and tags for allocation purposes
	int max_partitions = 0;
	int max_n_tags = 0;
	ptrdiff_t n_tags;
	ptrdiff_t n_elems_total = 0;
	ptrdiff_t n_elems_proc  = 0;
	ptrdiff_t n_elems_dom   = 0;
	ptrdiff_t n_elems_ghost = 0;
	int tags[GMSH_N_TAGS_MAX];
	// Check maximum partition number
	// Count number of elements within the partition
	while (fgets(line,sizeof(line),mesh_file)) {
		// Find Element list
		if (strstr(line,"$Elements")) {
			// Get n_elems_total
			fgets_checked(line,sizeof(line),mesh_file);
			n_elems_total = strtol(line,&endptr,10);

			while (fgets(line,sizeof(line),mesh_file)) {
				if (strstr(line,"$EndElements")) break;
				// Discard element ID and element type
				line_ptr[0] = line;
				discard_line_values(line_ptr,2);
				read_line_values_l(line_ptr,1,&n_tags,false);
				assert(n_tags>1 && n_tags<=GMSH_N_TAGS_MAX);
				if (n_tags==2) {
					max_partitions = 1;
					max_n_tags     = 2;
					n_elems_dom    = n_elems_total;
					n_elems_ghost  = 0;
					break;
				}
				read_line_values_i(line_ptr,n_tags,tags,false);
				// Tag 0: physical entity
				// Tag 1: elementary entity
				// Tag 2: number of partitions
				// Tag 3: partition it belongs to -> PARTITION_OWN
				// Tag 4-n_tag: partitions in which it is a ghost cell
				for (int i_tag=PARTITION_OWN; i_tag<n_tags; ++i_tag) {
					int i_part = tags[i_tag]; // i_part indexing is 1-based
					if (i_part>max_partitions) max_partitions = i_part;
					if (n_tags>max_n_tags) max_n_tags = (int) n_tags;
					if ((i_part-1)==mpi_rank) n_elems_dom++;
					if ((-i_part-1)==mpi_rank) n_elems_ghost++;
				}
			}
			break;
		}
	}
	rewind(mesh_file);
	ptrdiff_t n_elems_dom_sum = 0;
	printf("MPI_RANK:%d, n_elems_dom: %td \n", mpi_rank, n_elems_dom);
	MPI_Allreduce( &n_elems_dom, &n_elems_dom_sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  	if (n_elems_total!=n_elems_dom_sum) 
		EXIT_ERROR("Non matching number of total elements (%td) versus sum of domain elements (%td)", 
		n_elems_total, n_elems_dom_sum);

  	bool equal_part_proc = (max_partitions==mpi_size);
  	printf("nproc = %d, nparts = %d \n",mpi_size,max_partitions);
  	if (!equal_part_proc) {
  		EXIT_UNSUPPORTED;

  		//char new_mesh_name_full[STRLEN_MAX];
  		//sprintf(new_mesh_name_full, "%s_partitioned%d.msh", mesh_name_full, mpi_size);
  		//if(mpi_rank==1) {
  		//	printf("Non-matching number of partitions (%d) "
  		//		   "and number of processors (%d). \n", max_partitions, mpi_size);
  		//	char command[STRLEN_MAX];
  		//	//sprintf(command, "gmsh -string \"Merge /home/ddong/Codes/DPGSolver/build_debug_2D/meshes/%s; PartitionMesh %d;\" -o %s -0",
  		//	sprintf(command, "gmsh -bgm %s -part %d -o %s -0",
  		//			mesh_name_full, mpi_size, new_mesh_name_full);
  		//	printf("Partitioning using ParMetis. Command: %s \n", command);
  		//	system(command);
  		//}
  		//MPI_Barrier(MPI_COMM_WORLD);
  		//EXIT_UNSUPPORTED;
  		//mesh_reader_gmsh_parallel (new_mesh_name_full, d, mesh_data_l);
  		//return;
  	}

	n_elems_proc = n_elems_dom+n_elems_ghost;
	struct Element_Data* elem_data = constructor_Element_Data(n_elems_proc,max_n_tags); // keep

	ptrdiff_t row = 0;
	while (fgets(line,sizeof(line),mesh_file) != NULL) {
		if (strstr(line,"$Elements")) {
			while (fgets(line,sizeof(line),mesh_file) != NULL) {
				if (strstr(line,"$EndElements"))
					break;

				ptrdiff_t n_tags;
				line_ptr[0] = line;
				discard_line_values(line_ptr,1); // Element number
				read_line_values_i(line_ptr,1,&elem_data->elem_types->data[row],false); // Element type
				read_line_values_l(line_ptr,1,&n_tags,false); // Number of tags
				read_line_values_i(line_ptr,n_tags,get_row_Matrix_i(row,elem_data->elem_tags),false); // Tags

				bool skip_element = true;
				for (ptrdiff_t i_tag=PARTITION_OWN; i_tag<n_tags; ++i_tag) {
					ptrdiff_t i_part = tags[i_tag]; // i_part indexing is 1-based
					if ((i_part-1)==mpi_rank || (-i_part-1)==mpi_rank) skip_element = false;
				}
				if (skip_element) continue;

				int n_nodes = get_n_nodes(elem_data->elem_types->data[row]);
				resize_Vector_i(elem_data->node_nums->data[row],n_nodes);

				read_line_values_i(line_ptr,n_nodes,elem_data->node_nums->data[row]->data,true);
				reorder_nodes_gmsh(elem_data->elem_types->data[row],elem_data->node_nums->data[row]);

				if (row++ == n_elems_proc)
					EXIT_UNSUPPORTED;
			}
			printf("DONE ELEMENTS READ");
		}
		if (strstr(line,"$Nodes")) {
			printf("STARTING NODES READ\n");
			mesh_data_l->nodes = read_nodes(mesh_file,d); // keep
			printf("DONE NODES READ\n");
		}

		if (strstr(line,"Periodic")) {
			printf("STARTING PERIODIC READ\n");
			mesh_data_l->periodic_corr = read_periodic(mesh_file,d); // keep
			printf("DONE PERIODIC READ\n");
		}
	}
	mesh_data_l->elem_data = elem_data;

    //UNUSED(n_elems);
    UNUSED(d);
	fclose(mesh_file);
	printf("DONE GMSH READ");
}

static int get_mesh_format(FILE* mesh_file)
{
	// Ensure that the mesh is in the MSH4 format and is not binary

	int gmsh_mesh_version;
	int is_binary;

	char line[STRLEN_MAX];
	char* line_ptr[1] = {line}; // Need to pass pointer of pointer to next functions

	fgets_checked(line,sizeof(line),mesh_file);
	read_line_values_i(line_ptr,1,&gmsh_mesh_version,false);
	read_line_values_i(line_ptr,1,&is_binary,false);

	//int MESH_FORMAT = 4;
	//if(gmsh_mesh_version!=MESH_FORMAT) EXIT_ERROR("Incorrect Gmsh format %d. Should be MSH4.\n\n", gmsh_mesh_version);
	int BINARY_FORMAT = 0;
	if(is_binary!=BINARY_FORMAT) EXIT_ADD_SUPPORT; // Might want to support binary

	fgets_checked(line,sizeof(line),mesh_file);
	if (!strstr(line,"$EndMeshFormat")) {
		EXIT_ERROR("Should only have data here if file is binary");
	}
	else {
		return gmsh_mesh_version;
	}
}


struct Elementary_Entity {
	struct Vector_i* tags;      ///< The list of elementary entity tags.
	struct Vector_i* num_owners; ///< Number of partitions the entity belongs to.
	struct Matrix_i* owners;   ///< The list of the associated partitions. (0-based index)
	struct Matrix_i* phys_id;        ///< The list of physical ID. (BC or Connectivity)
};
static struct Elementary_Entity* constructor_Elementary_Entities (const ptrdiff_t n_entities, const int n_partitions)
{
	struct Elementary_Entity* elementary_entities = calloc(1,sizeof *elementary_entities); // returned

	elementary_entities->tags = constructor_empty_Vector_i(n_entities);
	elementary_entities->num_owners = constructor_empty_Vector_i(n_entities);
	elementary_entities->owners = constructor_empty_Matrix_i('R',n_entities, n_partitions);
	elementary_entities->phys_id = constructor_empty_Matrix_i('R',n_entities, MAX_PHYS_ID);

	return elementary_entities;
}
void destructor_Elementary_Entity (struct Elementary_Entity* elementary_entities)
{
	destructor_Vector_i(elementary_entities->tags);
	destructor_Vector_i(elementary_entities->num_owners);
	destructor_Matrix_i(elementary_entities->owners);
	destructor_Matrix_i(elementary_entities->phys_id);

	free(elementary_entities);
}
static struct Elementary_Entity* read_elementary_entities (FILE* mesh_file)
{
	struct Elementary_Entity* elementary_entities = NULL;

	// Read and construct the (Partitioned)Entities
	// Ensures that the number of partitions correspond to the number of processors
	int mpi_size; MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
	int mpi_rank; MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
	bool is_partitioned = (mpi_size > 1);
	char line[STRLEN_MAX];
	char* line_ptr[1] = {line};
	char* endptr = NULL;

	int num_partitions = 1;
	// PartitionedEntities header
	if (is_partitioned) {
		// Get the number of partitions and ensure it is equal to the number of processors
		fgets_checked(line,sizeof(line),mesh_file);
		num_partitions = (int) strtol(line,&endptr,10);
		if (num_partitions != mpi_size) 
			EXIT_ERROR("Number of partitions (%d) must equal number of processors (%d). \n", num_partitions, mpi_size);

		// Discard the rest of the header
		int num_ghost_partitions;
		fgets_checked(line,sizeof(line),mesh_file);
		num_ghost_partitions = (int) strtol(line,&endptr,10);
		skip_lines(mesh_file, num_ghost_partitions);
	}
	// Read number of points, curves, surfaces, and volumes
	// Sum them to obtain total number of entities
	fgets_checked(line,sizeof(line),mesh_file);
	struct Vector_i* num_entities_per_dim = constructor_empty_Vector_i(DMAX+1); //destructed
	line_ptr[0] = line;
	read_line_values_i(line_ptr,DMAX+1,num_entities_per_dim->data,false);
	int n_entities = sum_Vector_i(num_entities_per_dim);
	destructor_Vector_i(num_entities_per_dim);

	// Allocate entities
	elementary_entities = constructor_Elementary_Entities(n_entities, num_partitions); // returned
	int n_phys_temp;
	// Read entities
	for (ptrdiff_t i_ent = 0; i_ent < n_entities; i_ent++) {
		fgets_checked(line,sizeof(line),mesh_file);
		line_ptr[0] = line;
		// tag(int) ( parentDim(int) parentTag(int) numPartitions(int) partitionTag[...](int) )
		//		boxMinX(double) boxMinY(double) boxMinZ(double) boxMaxX(double) boxMaxY(double) boxMaxZ(double)
		//		numPhyiscals(unsigned long) physicalTag[...](int)
		read_line_values_i(line_ptr,1,&elementary_entities->tags->data[i_ent],false);

		if (is_partitioned) {
			discard_line_values(line_ptr,2); // parentDim(int) parentTag(int)
			read_line_values_i(line_ptr,1,&elementary_entities->num_owners->data[i_ent],false); // numPartitions(int)
			read_line_values_i(line_ptr,elementary_entities->num_owners->data[i_ent],
				get_row_Matrix_i(i_ent,elementary_entities->owners),true); // partitionTag[...](int)
		}
		else {
			elementary_entities->num_owners->data[i_ent] = 1;
			assert(mpi_rank == 0);
			set_row_Matrix_i(i_ent, elementary_entities->owners, &mpi_rank);
		}
		discard_line_values(line_ptr,6); // boxMinX/Y/Z, boxMaxX/Y/Z
		
		read_line_values_i(line_ptr,1,&n_phys_temp,false);
		if (n_phys_temp > 1) EXIT_ERROR("Entity with multiple physical entities (n = %d). "
										 "Only support 1 for now. \n", n_phys_temp);
		read_line_values_i(line_ptr,n_phys_temp,&elementary_entities->phys_id->data[i_ent],false);
		// ignore BREPVert/Curve/Surfaces
	}
	// Make sure that the next line correspond to the end of the entities section
	char entities_end_key[STRLEN_MAX];
	(is_partitioned) ? sprintf(entities_end_key,"$EndPartitionedEntities") : sprintf(entities_end_key,"$EndEntities");
	fgets_checked(line,sizeof(line),mesh_file);
	if (!strstr(line,entities_end_key)) EXIT_ERROR("Started reading $(Partitioned)Entities, but did not reach %s.\n", entities_end_key);

	if (elementary_entities == NULL) EXIT_ERROR("$Entities (mpi_size==1) or $PartitionedEntities (mpi_size>1) not found");

	return elementary_entities;
}

static int get_entity_partition(const ptrdiff_t entity_tag, const struct Elementary_Entity*const entities) {
	// Returns the partition of entity_tag
	// Note: Used for ghost cell. Assume the entity is only part of 1 partition
	//       An entity of dimension (d-1) can be shared by multiple partitions.
	for (ptrdiff_t i_ent = 0; i_ent < entities->tags->ext_0; i_ent++) {
		if (entities->tags->data[i_ent] == entity_tag) {
			if(entities->num_owners->data[i_ent] != 1) 
				EXIT_ERROR("get_entity_partition only makes sense when the entity is associated to only 1 partition");
			return get_val_Matrix_i(i_ent, 0, entities->owners);
		}
	}
	EXIT_ERROR("Element entity tag (%td) does not correspond to any $(Partitioned)Entity", entity_tag);
	return -1;
}

static int get_entity_phys_id(const ptrdiff_t entity_tag, const struct Elementary_Entity*const entities) {
	// Returns the phys_id of entity_tag
	for (ptrdiff_t i_ent = 0; i_ent < entities->tags->ext_0; i_ent++) {
		if (entities->tags->data[i_ent] == entity_tag) {
			return entities->phys_id->data[i_ent];
		}
	}
	EXIT_ERROR("Did not find element entity tag (%td) when fetching phys_id.\n", entity_tag);
	return -1;
}
static struct Vector_i* get_elements_id (FILE* mesh_file, struct Elementary_Entity*const entities, const int mpi_rank)
{
	// Constructs the elements_id vector and get the list of elements that should be read by (mpi_rank)
	// Using push_back with the current implementation result in O(n^2) operations since reallocation 
	// copies the entire vector every time.
	// For large grids, might need to change the vector_T resize or push_back operation.

	// Get the entity tag and the number of elements associated with the entity
	// If the entity is owned by the partition/processor
	//		Add the number of elements to the count
	//		Add the element tags to the vector of element ID.
	// NOTE: An entity may belong to multiple processors. 
	//       Therefore, the total number of elements is not the sum of elements per processor.
	ptrdiff_t n_elems_partition = 0;

	char line[STRLEN_MAX];
	char* line_ptr[1] = {line}; // Need to pass point to next functions
	ptrdiff_t n_ent_blocks;
	ptrdiff_t n_elems_total;
	ptrdiff_t entity_tag;
	ptrdiff_t n_elems_entity;
	// Currently, only Vector_i are available. Therefore, we can't read a vector of ptrdiff_t.
	// If there are more than 2 billion elements, this will break.
	int element_id;


	// Used to check that the correct total number of elements have been read.
	ptrdiff_t n_elems_double_count = 0;

	line_ptr[0] = ptr_fgets_checked(line,sizeof(line),mesh_file);
	read_line_values_l(line_ptr,1,&n_ent_blocks,false);
	read_line_values_l(line_ptr,1,&n_elems_total,false);
	struct Vector_i* elements_id = constructor_default_Vector_i(); // returned
	ptrdiff_t i_elem = 0;
	for (ptrdiff_t i_ent_block = 0; i_ent_block < n_ent_blocks; i_ent_block++) {
		line_ptr[0] = ptr_fgets_checked(line,sizeof(line),mesh_file);
		read_line_values_l(line_ptr,1,&entity_tag,false);
		discard_line_values(line_ptr, 2); // elem_dim, elem_type
		read_line_values_l(line_ptr,1,&n_elems_entity,false);
		bool found_entity = false;
		bool own_entity   = false;
		// Check if partition owns the entity
		for (ptrdiff_t i_ent = 0; i_ent < entities->tags->ext_0; i_ent++) {
			if (entities->tags->data[i_ent] == entity_tag) {
				found_entity = true;
				n_elems_double_count += (entities->num_owners->data[i_ent]-1) * n_elems_entity;
				for (ptrdiff_t i_owner = 0; i_owner < entities->num_owners->data[i_ent]; i_owner++) {
					if(get_val_Matrix_i(i_ent, i_owner, entities->owners) == mpi_rank) {
						own_entity = true;
						break;
					}
				}
				break;
			}
		}
		if (!found_entity) EXIT_ERROR("Element entity tag (%td) does not correspond to any $(Partitioned)Entity", entity_tag);

		if (own_entity) {
			n_elems_partition += n_elems_entity;
			for (ptrdiff_t i_elem_ent = 0; i_elem_ent < n_elems_entity; i_elem_ent++) {
				line_ptr[0] = ptr_fgets_checked(line,sizeof(line),mesh_file);
				read_line_values_i(line_ptr,1,&element_id,true);
				push_back_Vector_i(elements_id,element_id,false,false);
				i_elem++;
			}
		}
		else {
			skip_lines(mesh_file, (int) n_elems_entity);
		}
	}
	assert(n_elems_partition == i_elem);

	// Sanity check
	ptrdiff_t n_elems_partition_sum = 0;
  	MPI_Allreduce(&n_elems_partition, &n_elems_partition_sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	ptrdiff_t n_elems_sum = n_elems_partition_sum-n_elems_double_count;
	if(n_elems_sum != n_elems_total) 
  		EXIT_ERROR("Non matching number of total elements (%td) versus sum of domain elements (%td)", 
  		n_elems_total, n_elems_sum);

	resize_Vector_i(elements_id, n_elems_partition);
	return elements_id;
}

static struct Vector_i* get_ghost_elements_id (FILE* mesh_file, const int mpi_rank)
{
	// Constructs the ghost_elements_id vector and get the list of ghost elements that should be read by (mpi_rank)
	// Using push_back since we do not expect a massive amount of ghost cells
	ptrdiff_t n_ghost_elements_total;
	ptrdiff_t n_ghost_elements_partition = 0;

	char line[STRLEN_MAX];
	char* line_ptr[1] = {line}; // Need to pass point to next functions
	line_ptr[0] = ptr_fgets_checked(line,sizeof(line),mesh_file);
	read_line_values_l(line_ptr,1,&n_ghost_elements_total,false);

	struct Vector_i* ghost_elements_id = constructor_default_Vector_i(); // returned
	// If there are more than 2 billion elements. Might need ptrdiff_t Vector.
	int element_id, n_ghost_owner, actual_element_owner, ghost_element_owner;
	for (int i_ghost = 0; i_ghost < n_ghost_elements_total; i_ghost++) {
		line_ptr[0] = ptr_fgets_checked(line,sizeof(line),mesh_file);
		// Need to follow up with GMSH team to make sure that 
		// actual_element_owner and ghost_element_owner are indeed in that order
		read_line_values_i(line_ptr,1,&element_id,true);
		read_line_values_i(line_ptr,1,&actual_element_owner,true); // 0-based
		read_line_values_i(line_ptr,1,&n_ghost_owner,false);
		for (int i_part = 0; i_part < n_ghost_owner; i_part++) {
			read_line_values_i(line_ptr,1,&ghost_element_owner,true); // 0-based
			if (ghost_element_owner == mpi_rank) {
				n_ghost_elements_partition++;
				push_back_Vector_i(ghost_elements_id,element_id,false,false);
				break;
			}
		}
	}
	return ghost_elements_id;
}

static void read_elements_partition 
	(FILE* mesh_file, const struct Element_Data*const elem_data, const struct Vector_i*const elements_id, 
	 const struct Vector_i*const ghost_elements_id, const struct Elementary_Entity*const entities, int mpi_rank)
{ 
	char line[STRLEN_MAX];
	char* line_ptr[1] = {line};
	line_ptr[0] = ptr_fgets_checked(line,sizeof(line),mesh_file);

	ptrdiff_t n_ent_blocks, n_elems_entity;
	int entity_tag, elem_id, elem_type, owner;
	read_line_values_l(line_ptr,1,&n_ent_blocks,false);
	//ptrdiff_t n_elems_total;
	//read_line_values_l(line_ptr,1,&n_elems_total,false);

	ptrdiff_t i_elem = 0;
	ptrdiff_t i_ghost = 0;
	bool own_elem, own_ghost;
	for (ptrdiff_t i_ent_block = 0; i_ent_block < n_ent_blocks; i_ent_block++) {

		line_ptr[0] = ptr_fgets_checked(line,sizeof(line),mesh_file);

		read_line_values_i(line_ptr,1,&entity_tag,false);
		discard_line_values(line_ptr,1); // dimEntity
		read_line_values_i(line_ptr,1,&elem_type,false);
		read_line_values_l(line_ptr,1,&n_elems_entity,false);

		for (ptrdiff_t i_elem_ent = 0; i_elem_ent < n_elems_entity; i_elem_ent++) {
			line_ptr[0] = ptr_fgets_checked(line,sizeof(line),mesh_file);
			read_line_values_i(line_ptr,1,&elem_id,true);

			own_elem = false; own_ghost = false;
			own_elem = find_val_Vector_i((struct const_Vector_i*)elements_id, elem_id, true);
			if (ghost_elements_id!=NULL && !own_elem) {
				own_ghost = find_val_Vector_i((struct const_Vector_i*)ghost_elements_id, entity_tag, true);
			}
			if (own_elem) {
				owner = mpi_rank;
			}
			else {
				owner = get_entity_partition(entity_tag, entities);
			}
			if (own_elem || own_ghost) {
				ptrdiff_t* index = (own_elem) ? &i_elem : &i_ghost;

				elem_data->elem_global_id->data[*index]  = elem_id;
				elem_data->elem_types->data[*index]      = elem_type;
				elem_data->elem_owner_part->data[*index] = owner;
				int* tags_row = get_row_Matrix_i(*index, elem_data->elem_tags);
				tags_row[0] = get_entity_phys_id(entity_tag, entities);

				int n_nodes = get_n_nodes(elem_data->elem_types->data[*index]);
				resize_Vector_i(elem_data->node_nums->data[*index],n_nodes);

				read_line_values_i(line_ptr,n_nodes,elem_data->node_nums->data[*index]->data,true);
				reorder_nodes_gmsh(elem_data->elem_types->data[*index],elem_data->node_nums->data[*index]);
							elem_data->elem_global_id->data[i_elem] = entity_tag;

				(*index)++;
			}


		}
	}
	return;
}

static struct Matrix_d* read_nodes_partition (FILE* mesh_file, const struct Vector_i*const node_nums_unique, const int dim)
{
	struct Matrix_d* nodes = constructor_empty_Matrix_d('R',node_nums_unique->ext_0,dim); // returned

	char line[STRLEN_MAX];
	char* line_ptr[1] = {line};

	ptrdiff_t n_ent_blocks, n_nodes_entity;

	struct const_Vector_i* const_node_nums_unique = (struct const_Vector_i*) node_nums_unique; // find_val_Vector_i uses const_Vector_i

	line_ptr[0] = ptr_fgets_checked(line,sizeof(line),mesh_file);
	read_line_values_l(line_ptr,1,&n_ent_blocks,false);

	int i_node = 0;
	int node_id;
	for (ptrdiff_t i_ent_block = 0; i_ent_block < n_ent_blocks; i_ent_block++) {
		line_ptr[0] = ptr_fgets_checked(line,sizeof(line),mesh_file);
		discard_line_values(line_ptr,3); // tagEntity, dimEntity, typeNode
		read_line_values_l(line_ptr,1,&n_nodes_entity,false);
		for (ptrdiff_t i_node_ent = 0; i_node_ent < n_nodes_entity; i_node_ent++) {
			line_ptr[0] = ptr_fgets_checked(line,sizeof(line),mesh_file);
			read_line_values_i(line_ptr,1,&node_id,true);
			bool own_node = find_val_Vector_i(const_node_nums_unique, node_id, true);
			if (own_node) {
				char* endptr = NULL;
				double* node_row = get_row_Matrix_d(i_node,nodes);
				i_node++;
				for (int i_dim = 0; i_dim < dim; i_dim++) {
					node_row[i_dim] = strtod(*line_ptr,&endptr);
					line_ptr[0] = endptr;
				}
			}
		}
	}
	for (int i = 0; i < node_nums_unique->ext_0; i++) {
		double* node_row = get_row_Matrix_d(i,nodes);
		for (int j = 0; j < dim; j++) {
			printf("%f ", node_row[j]);
		}
		printf("\n ");
	}
	return nodes;
}

static void mesh_reader_gmsh4 (const char*const mesh_name_full, const int d, struct Mesh_Data_l*const mesh_data_l)
{
	// Assume the mesh is given with the following mandatory tags in the following order
	// MeshFormat, Entities, (PartitionedEntities), Nodes, Elements, Periodic, (GhostElements)
	// Although each section may appear more than once, we assume it only appears once.

	int mpi_size; MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
	int mpi_rank; MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);

	mesh_data_l->nodes         = NULL;
	mesh_data_l->elem_data     = NULL;
	mesh_data_l->periodic_corr = NULL;

	char line[STRLEN_MAX];
	//char* line_ptr[1] = {line}; // Need to pass point to next functions
	//char* endptr = NULL;

	printf("Opening mesh named: %s\n", mesh_name_full);
	FILE* mesh_file = fopen_checked(mesh_name_full); // closed

	// Ensure that the mesh is in the MSH4 format and is not binary
	fgets_checked(line,sizeof(line),mesh_file);
	if (!strstr(line,"$MeshFormat")) EXIT_ERROR("Mandatory first file section $MeshFormat not present.");
	int mesh_format = get_mesh_format(mesh_file);
	if (mesh_format == 2) {
		fclose(mesh_file);
		mesh_reader_gmsh(mesh_name_full, d, mesh_data_l);
		return;
	}
	else if (mesh_format == 4) {
	}
	else {
		EXIT_ERROR("Invalid Gmsh file format version %d.\n", mesh_format);
	}

  	//ptrdiff_t n_elems_total = 0;
  	ptrdiff_t n_elems_proc  = 0;
  	ptrdiff_t n_elems_dom   = 0;
  	ptrdiff_t n_elems_ghost = 0;

  	// Pre-processing loop:
	//		Ensures that the number of partitions correspond to the number of processors
	// 		Read and construct the (Partitioned)Entities
	// 		Count number of elements that belongs to the processor should read
	//      Add number of ghost ghost elements to the count

	// Note about efficiency and memory limitations:
	// Chose to use push_back to add elements that the partition owns.
	// This takes O(n^2), which is ok since it only happens once.
	// Need to consider a more efficient push_back/resize function for larger grids (more than ~30k elements per partition)

	// Get entities from different sections depending on serial or parallel run.
	char entities_start_key[STRLEN_MAX];
	(mpi_size>1) ? sprintf(entities_start_key,"$PartitionedEntities") : sprintf(entities_start_key,"$Entities");

	bool hasGhosts = false;
	struct Elementary_Entity* entities = NULL;
	struct Vector_i* elements_id = NULL;
	struct Vector_i* ghost_elements_id = NULL;
	fpos_t entities_position, nodes_position, elements_position, ghostelements_position;
  	while (fgets(line,sizeof(line),mesh_file)) {
		if (strstr(line,entities_start_key)) {
			if(fgetpos(mesh_file,&entities_position) != 0) EXIT_UNSUPPORTED;
			entities = read_elementary_entities(mesh_file); // destructed
		}
  		else if (strstr(line,"$Nodes")) {
			if(fgetpos(mesh_file,&nodes_position) != 0) EXIT_UNSUPPORTED;
		}
  		else if (strstr(line,"$Elements")) {
			if(fgetpos(mesh_file,&elements_position) != 0) EXIT_UNSUPPORTED;
			if (entities == NULL) EXIT_ERROR("Found $Elements before %s. \n", entities_start_key);
			elements_id = get_elements_id(mesh_file,entities, mpi_rank); // NEED TO BE MOVED OR DESTROYED
		}
  		else if (strstr(line,"$GhostElements")) {
			if (fgetpos(mesh_file,&ghostelements_position) != 0) EXIT_UNSUPPORTED;
			if (mpi_size == 1) EXIT_ERROR("There should not be a ghost element section for serial run.\n");
			ghost_elements_id = get_ghost_elements_id(mesh_file, mpi_rank); // NEED TO BE MOVED OR DESTROYED
			hasGhosts = true;
		}
	}
	// Sort the vectors since we will be searching them when reading in the elements.
	if (elements_id == NULL) EXIT_ERROR("Elements section not found.\n");
	sort_Vector_i(elements_id);
	n_elems_dom = elements_id->ext_0;
	if (hasGhosts) {
		sort_Vector_i(ghost_elements_id);
		n_elems_ghost = ghost_elements_id->ext_0;
	}

	n_elems_proc = n_elems_dom+n_elems_ghost;
	struct Element_Data* elem_data = constructor_Element_Data(n_elems_proc,2); // keep

	ptrdiff_t index_first_ghost = elements_id->ext_0; UNUSED(index_first_ghost);
	//push_back_Vector_Vector_i(elements_id, ghost_elements_id);
	//copy_data_Vector_i_Vector_i(elements_id, elem_data->elem_global_id);

	// Go to elements section and read elements that belong to this partition
	if(fsetpos(mesh_file,&elements_position) != 0) EXIT_UNSUPPORTED;

	elem_data->ind_g = elements_id->ext_0;
	// Reads in elem_global_id, elem_owner_part, elem_types, node_nums
	read_elements_partition(mesh_file, elem_data, elements_id, ghost_elements_id, entities, mpi_rank);

	// Flatten node_nums, sort and remove repeated nodes
	struct Vector_i* node_nums_repeated = collapse_Multiarray_Vector_i(elem_data->node_nums); // destructed
	struct Vector_i* node_nums_unique = constructor_default_Vector_i();//(n_elems_proc); // destructed
	sort_Vector_i(node_nums_repeated);
	int last_node_id = -1;
	int current_node_id = node_nums_repeated->data[0];
	push_back_Vector_i(node_nums_unique, current_node_id, false, false);
	for (ptrdiff_t i_node = 0; i_node < node_nums_repeated->ext_0; i_node++) {
		last_node_id = current_node_id;
		current_node_id = node_nums_repeated->data[i_node];
		if (current_node_id == last_node_id) continue;
		push_back_Vector_i(node_nums_unique, current_node_id, false, false);
	}
	destructor_Vector_i(node_nums_repeated);

	// Go to nodes section and read nodes that appear in the elem_data->node_nums
	if(fsetpos(mesh_file,&nodes_position) != 0) EXIT_UNSUPPORTED;

	mesh_data_l->nodes = read_nodes_partition(mesh_file, node_nums_unique, d);
    mesh_data_l->elem_data = elem_data;
	

	// Read the elements information


	printf("DONE READING GMSH MSH4.\n");
	return;


	destructor_Vector_i(elements_id);
	destructor_Vector_i(ghost_elements_id);
	destructor_Elementary_Entity(entities);
    //UNUSED(n_elems);
    UNUSED(d);
	fclose(mesh_file);
	printf("DONE GMSH READ");
}
