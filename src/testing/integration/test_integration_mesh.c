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
/** \file
 */

#include <string.h>

#include "test_base.h"
#include "test_support_vector.h"
#include "test_support_matrix.h"
#include "test_support_multiarray.h"

#include "macros.h"
#include "definitions_mesh.h"
#include "definitions_alloc.h"

#include "multiarray.h"
#include "matrix.h"
#include "vector.h"

#include "intrusive.h"
#include "mesh.h"
#include "mesh_readers.h"
#include "mesh_connectivity.h"
#include "mesh_vertices.h"
#include "const_cast.h"

// Static function declarations ************************************************************************************* //

/** \brief Contructs a \ref Mesh_Input\*.
 *	\return Standard. */
static struct Mesh_Input* constructor_Mesh_Input
	(const char*const mesh_name ///< The test mesh name.
	);

/// \brief Destructs a \ref Mesh_Input\*.
static void destructor_Mesh_Input
	(struct Mesh_Input* mesh_input ///< Standard.
	);

/** \brief Compare members of the \ref Mesh with expected values.
 *  \return 1 if tests passed. */
static bool compare_members_Mesh
	(struct Test_Info*const test_info, ///< \ref Test_Info.
	 const char*const mesh_name_full,  ///< Defined in \ref Mesh_Input.
	 const struct Mesh*const mesh      ///< \ref Mesh.
	);

// Interface functions ********************************************************************************************** //

/** \test Performs integration testing for the mesh processing (\ref test_integration_mesh.c).
 *  \return 0 on success.
 *
 *  Compares members of the following containers with their expected values:
 *  - \ref Mesh_Data;
 *  - \ref Mesh_Connectivity;
 *  - \ref Mesh_Vertices.
 */
int main
	(int nargc,  ///< Standard.
	 char** argv ///< Standard.
	)
{
	assert_condition_message(nargc == 2,"Invalid number of input arguments");
	const char* mesh_name = argv[1];

	struct Test_Info test_info = { .n_warn = 0, };

	struct Mesh_Input* mesh_input = constructor_Mesh_Input(mesh_name); // destructed
	struct Mesh* mesh             = constructor_Mesh(mesh_input,NULL); // destructed

	const bool pass = compare_members_Mesh(&test_info,mesh_input->mesh_name_full,mesh);

	destructor_Mesh(mesh);
	destructor_Mesh_Input(mesh_input);

	assert_condition(pass);
	output_warning_count(&test_info);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Container for data of the \ref Mesh which is to be tested.
struct Mesh_Test_Data {
	struct Vector_i*            elem_per_dim;  ///< Defined in \ref Mesh_Data.
	struct Matrix_d*            nodes;         ///< Defined in \ref Mesh_Data.
	struct Vector_i*            elem_types;    ///< Defined in \ref Mesh_Data.
	struct Matrix_i*            elem_tags;     ///< Defined in \ref Mesh_Data.
	struct Multiarray_Vector_i* node_nums;     ///< Defined in \ref Mesh_Data.
	struct Matrix_i*            periodic_corr; ///< Defined in \ref Mesh_Data.

	struct Multiarray_Vector_i* v_to_v;  ///< Defined in \ref Mesh_Connectivity.
	struct Multiarray_Vector_i* v_to_lf; ///< Defined in \ref Mesh_Connectivity.

	struct Vector_i*            ve_curved;   ///< Defined in \ref Mesh_Vertices.
	struct Vector_i*            ve_boundary; ///< Defined in \ref Mesh_Vertices.
	struct Multiarray_Vector_i* ve_bc;       ///< Defined in \ref Mesh_Vertices.
};

/** \brief Set the members of the \ref Mesh_Input\*.
 *  \todo Check if `geom_spec` is being used and remove it if not. */
static void set_Mesh_Input
	(struct Mesh_Input*const mesh_input, ///< \ref Mesh_Input.
	 const int d,                        ///< Defined in \ref Mesh_Input.
	 const int domain_type,              ///< Defined in \ref Mesh_Input.
	 const bool mesh_unrealistic,        ///< Defined in \ref Mesh_Input.
	 const char* mesh_name,              ///< The mesh name without the relative path.
	 const char* geom_name_,             ///< Defined in \ref Mesh_Input.
	 const char* geom_spec_,             ///< Defined in \ref Mesh_Input.
	 const char* input_path_             ///< Defined in \ref Mesh_Input.
	);

/** \brief Constructs a \ref Mesh_Test_Data\*.
 *	\return Standard. */
static struct Mesh_Test_Data* constructor_Mesh_Test_Data
	(const char*const mesh_name_full ///< Defined in \ref Mesh_Input.
	);

/// \brief Destructor for \ref Mesh_Test_Data\*.
static void destructor_Mesh_Test_Data
	(struct Mesh_Test_Data* mesh_test_data ///< \ref Mesh_Test_Data.
	);

static struct Mesh_Input* constructor_Mesh_Input (const char*const mesh_name)
{
	struct Mesh_Input* mesh_input = calloc(1,sizeof *mesh_input); // free

	mesh_input->mesh_name_full = malloc(STRLEN_MAX * sizeof *mesh_input->mesh_name_full); // free
	mesh_input->geom_name      = malloc(STRLEN_MAX * sizeof *mesh_input->geom_name);      // free
	mesh_input->geom_spec      = malloc(STRLEN_MAX * sizeof *mesh_input->geom_spec);      // free
	mesh_input->input_path     = malloc(STRLEN_MAX * sizeof *mesh_input->input_path);     // free

	if (strstr(mesh_name,"curved_2d_mixed.msh")) {
		set_Mesh_Input(mesh_input,2,DOM_CURVED,true,mesh_name,"n-cylinder_hollow_section","",
		               "../input_files/euler/steady/supersonic_vortex/");
	} else if (strstr(mesh_name,"straight_2d_quad_periodic.msh")) {
		set_Mesh_Input(mesh_input,2,DOM_STRAIGHT,false,mesh_name,"","","");
	} else {
		EXIT_ERROR("Unsupported: %s\n",mesh_name);
	}

	return mesh_input;
}

static void destructor_Mesh_Input (struct Mesh_Input* mesh_input)
{
	free((void*)mesh_input->mesh_name_full);
	free((void*)mesh_input->geom_name);
	free((void*)mesh_input->geom_spec);
	free((void*)mesh_input->input_path);

	free(mesh_input);
}

static bool compare_members_Mesh
	(struct Test_Info*const test_info, const char*const mesh_name_full, const struct Mesh*const mesh)
{
	UNUSED(test_info);
	bool pass = true;

	struct Mesh_Test_Data* mesh_test_data = constructor_Mesh_Test_Data(mesh_name_full); // destructed

	// Mesh_Data
	struct Vector_i*            elem_per_dim  = (struct Vector_i*)            mesh->mesh_data->elem_per_dim;
	struct Matrix_d*            nodes         = (struct Matrix_d*)            mesh->mesh_data->nodes;
	struct Vector_i*            elem_types    = (struct Vector_i*)            mesh->mesh_data->elem_types;
	struct Matrix_i*            elem_tags     = (struct Matrix_i*)            mesh->mesh_data->elem_tags;
	struct Multiarray_Vector_i* node_nums     = (struct Multiarray_Vector_i*) mesh->mesh_data->node_nums;
	struct Matrix_i*            periodic_corr = (struct Matrix_i*)            mesh->mesh_data->periodic_corr;

	if ((diff_Vector_i(elem_per_dim,mesh_test_data->elem_per_dim) != 0)                    ||
	    (diff_Matrix_d(nodes,mesh_test_data->nodes,NODETOL_MESH) != 0)                     ||
	    (diff_Vector_i(elem_types,mesh_test_data->elem_types) != 0)                        ||
	    (diff_Matrix_i(elem_tags,mesh_test_data->elem_tags) != 0)                          ||
	    (diff_Multiarray_Vector_i(node_nums,mesh_test_data->node_nums) != 0)               ||
	    (periodic_corr && diff_Matrix_i(periodic_corr,mesh_test_data->periodic_corr) != 0))
	{
		pass = false;
		expect_condition(pass,"mesh data");

		print_diff_Vector_i(elem_per_dim,mesh_test_data->elem_per_dim);
		print_diff_Matrix_d(nodes,mesh_test_data->nodes,NODETOL_MESH);
		print_diff_Vector_i(elem_types,mesh_test_data->elem_types);
		print_diff_Matrix_i(elem_tags,mesh_test_data->elem_tags);
		print_diff_Multiarray_Vector_i(node_nums,mesh_test_data->node_nums);
		print_diff_Matrix_i(periodic_corr,mesh_test_data->periodic_corr);
	}

	// Mesh_Connectivity
	struct Multiarray_Vector_i* v_to_v  = (struct Multiarray_Vector_i*) mesh->mesh_conn->v_to_v;
	struct Multiarray_Vector_i* v_to_lf = (struct Multiarray_Vector_i*) mesh->mesh_conn->v_to_lf;

	if ((diff_Multiarray_Vector_i(v_to_v,mesh_test_data->v_to_v) != 0)   ||
	    (diff_Multiarray_Vector_i(v_to_lf,mesh_test_data->v_to_lf) != 0))
	{
		pass = false;
		expect_condition(pass,"mesh connectivity");

		print_diff_Multiarray_Vector_i(v_to_v,mesh_test_data->v_to_v);
		print_diff_Multiarray_Vector_i(v_to_lf,mesh_test_data->v_to_lf);
	}

	// Mesh_Vertices
	struct Vector_i* ve_curved        = (struct Vector_i*)            mesh->mesh_vert->ve_curved;
	struct Vector_i* ve_boundary      = (struct Vector_i*)            mesh->mesh_vert->ve_boundary;
	struct Multiarray_Vector_i* ve_bc = (struct Multiarray_Vector_i*) mesh->mesh_vert->ve_bc;

	if ((diff_Vector_i(ve_curved,mesh_test_data->ve_curved) != 0)     ||
	    (diff_Vector_i(ve_boundary,mesh_test_data->ve_boundary) != 0) ||
	    (diff_Multiarray_Vector_i(ve_bc,mesh_test_data->ve_bc) != 0))
	{
		pass = false;
		expect_condition(pass,"mesh vertices");

		print_diff_Vector_i(ve_curved,mesh_test_data->ve_curved);
		print_diff_Vector_i(ve_boundary,mesh_test_data->ve_boundary);
		print_diff_Multiarray_Vector_i(ve_bc,mesh_test_data->ve_bc);
	}

	destructor_Mesh_Test_Data(mesh_test_data);

	return pass;
}

// Level 1 ********************************************************************************************************** //

static void set_Mesh_Input
	(struct Mesh_Input*const mesh_input, const int d, const int domain_type, const bool mesh_unrealistic,
	 const char* mesh_name, const char* geom_name_, const char* geom_spec_, const char* input_path_)
{
	const_cast_i(&mesh_input->d,d);
	const_cast_i(&mesh_input->domain_type,domain_type);
	const_cast_b(&mesh_input->mesh_unrealistic,mesh_unrealistic);

	char* mesh_name_full = (char*) mesh_input->mesh_name_full;
	strcpy(mesh_name_full,"../testing/integration/mesh/");
	strcat(mesh_name_full,mesh_name);

	char* geom_name = (char*) mesh_input->geom_name;
	strcpy(geom_name,geom_name_);

	char* geom_spec = (char*) mesh_input->geom_spec;
	strcpy(geom_spec,geom_spec_);
	const_cast_c1(&mesh_input->geom_spec,geom_spec);

	char* input_path = (char*) mesh_input->input_path;
	strcpy(input_path,input_path_);
}

static struct Mesh_Test_Data* constructor_Mesh_Test_Data (const char*const mesh_name_full)
{
	struct Mesh_Test_Data* mesh_test_data = calloc(1,sizeof *mesh_test_data); // free

	char mesh_data_name_full[STRLEN_MAX];
	strcpy(mesh_data_name_full,mesh_name_full);
	strcat(mesh_data_name_full,".data");

	mesh_test_data->elem_per_dim = constructor_file_name_Vector_i("elem_per_dim",mesh_data_name_full); // destructed
	mesh_test_data->nodes        = constructor_file_name_Matrix_d("nodes",mesh_data_name_full);        // destructed
	mesh_test_data->elem_types   = constructor_file_name_Vector_i("elem_types",mesh_data_name_full);   // destructed
	mesh_test_data->elem_tags    = constructor_file_name_Matrix_i("elem_tags",mesh_data_name_full);    // destructed
	mesh_test_data->node_nums    =
		constructor_file_name_Multiarray_Vector_i("node_nums",mesh_data_name_full); // destructed
	if (strstr(mesh_name_full,"periodic"))
		mesh_test_data->periodic_corr =
			constructor_file_name_Matrix_i("periodic_corr",mesh_data_name_full); // destructed
	else
		mesh_test_data->periodic_corr = NULL;

	mesh_test_data->v_to_v  = constructor_file_name_Multiarray_Vector_i("v_to_v",mesh_data_name_full);  // destructed
	mesh_test_data->v_to_lf = constructor_file_name_Multiarray_Vector_i("v_to_lf",mesh_data_name_full); // destructed

	mesh_test_data->ve_curved   = constructor_file_name_Vector_i("ve_curved",mesh_data_name_full);        // destructed
	mesh_test_data->ve_boundary = constructor_file_name_Vector_i("ve_boundary",mesh_data_name_full);      // destructed
	mesh_test_data->ve_bc       = constructor_file_name_Multiarray_Vector_i("ve_bc",mesh_data_name_full); // destructed

	return mesh_test_data;
}

static void destructor_Mesh_Test_Data (struct Mesh_Test_Data* mesh_test_data)
{
	destructor_Vector_i(mesh_test_data->elem_per_dim);
	destructor_Matrix_d(mesh_test_data->nodes);
	destructor_Vector_i(mesh_test_data->elem_types);
	destructor_Matrix_i(mesh_test_data->elem_tags);
	destructor_Multiarray_Vector_i(mesh_test_data->node_nums);
	if (mesh_test_data->periodic_corr)
		destructor_Matrix_i(mesh_test_data->periodic_corr);

	destructor_Multiarray_Vector_i(mesh_test_data->v_to_v);
	destructor_Multiarray_Vector_i(mesh_test_data->v_to_lf);

	destructor_Vector_i(mesh_test_data->ve_curved);
	destructor_Vector_i(mesh_test_data->ve_boundary);
	destructor_Multiarray_Vector_i(mesh_test_data->ve_bc);

	free(mesh_test_data);
}
