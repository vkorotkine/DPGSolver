// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/**	\file
 */

#include "test_integration_mesh.h"

#include <string.h>

#include "Macros.h"
#include "intrusive.h"
#include "element.h"
#include "multiarray.h"
#include "mesh.h"
#include "const_cast.h"
#include "allocators.h"

#include "constants_mesh.h"
#include "constants_alloc.h"

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
 *	\return 1 if tests passed. */
static bool compare_members_Mesh
	(struct Test_Info*const test_info, ///< Defined in \ref test_integration_mesh.
	 const char*const mesh_name,       ///< Defined in \ref test_integration_mesh.
	 const struct Mesh*const mesh      ///< \ref Mesh.
	);

// Interface functions ********************************************************************************************** //

void test_integration_mesh (struct Test_Info*const test_info, const char*const mesh_name)
{
	struct Mesh_Input* mesh_input         = constructor_Mesh_Input(mesh_name);       // destructed
	struct const_Intrusive_List* elements = constructor_Element_List(mesh_input->d); // destructed
	struct Mesh* mesh                     = constructor_Mesh(mesh_input,elements);   // destructed

	const bool pass = compare_members_Mesh(test_info,mesh_name,mesh);

	char test_name[STRLEN_MAX];
	strcpy(test_name,"Mesh - ");
	strcat(test_name,mesh_name);

	test_increment_and_print(test_info, pass,test_name);

	destructor_Mesh(mesh);
	destructor_Elements((struct Intrusive_List*)elements);
	destructor_Mesh_Input(mesh_input);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Container for data of the \ref Mesh which is to be tested.
struct Mesh_Test_Data {
	struct Multiarray_Vector_i* v_to_v;  ///< Defined in \ref Mesh_Connectivity.
	struct Multiarray_Vector_i* v_to_lf; ///< Defined in \ref Mesh_Connectivity.
};

/** \brief Set the members of the \ref Mesh_Input\*.
 *	\todo Check if `geom_spec` is being used and remove it if not. */
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
	(const char*const mesh_name ///< The test mesh name.
	);

/// \brief Destructor for \ref Mesh_Test_Data\*.
static void destructor_Mesh_Test_Data
	(struct Mesh_Test_Data* mesh_test_data ///< \ref Mesh_Test_Data.
	);

static struct Mesh_Input* constructor_Mesh_Input (const char*const mesh_name)
{
	struct Mesh_Input* mesh_input = malloc(sizeof *mesh_input); // free

	mesh_input->mesh_name_full = mallocator(CHAR_T,1,STRLEN_MAX); // destructed
	mesh_input->geom_name      = mallocator(CHAR_T,1,STRLEN_MAX); // destructed
	mesh_input->geom_spec      = mallocator(CHAR_T,1,STRLEN_MAX); // destructed
	mesh_input->input_path     = mallocator(CHAR_T,1,STRLEN_MAX); // destructed

	if (strstr(mesh_name,"curved_2d_mixed.msh")) {
		set_Mesh_Input(mesh_input,2,DOM_CURVED,true,mesh_name,"n-cylinder_hollow_section","",
		               "../input_files/euler/internal/supersonic_vortex/");
	} else if (strstr(mesh_name,"straight_2d_quad_periodic.msh")) {
		set_Mesh_Input(mesh_input,2,DOM_STRAIGHT,false,mesh_name,"","","");
	} else {
		EXIT_UNSUPPORTED;
	}

	return mesh_input;
}

static void destructor_Mesh_Input (struct Mesh_Input* mesh_input)
{
	deallocator((void*)mesh_input->mesh_name_full,CHAR_T,1,STRLEN_MAX);
	deallocator((void*)mesh_input->geom_name,     CHAR_T,1,STRLEN_MAX);
	deallocator((void*)mesh_input->geom_spec,     CHAR_T,1,STRLEN_MAX);
	deallocator((void*)mesh_input->input_path,    CHAR_T,1,STRLEN_MAX);

	free(mesh_input);
}

static bool compare_members_Mesh
	(struct Test_Info*const test_info, const char*const mesh_name, const struct Mesh*const mesh)
{
	struct Mesh_Test_Data* mesh_test_data = constructor_Mesh_Test_Data(mesh_name); // destructed

	bool pass = 1;

	// Mesh_Data
	struct Vector_i* elem_per_dim = (struct Vector_i*) mesh->mesh_data->elem_per_dim;
	struct Matrix_d* nodes        = (struct Matrix_d*) mesh->mesh_data->nodes;

	test_print_warning(test_info,"Not yet testing Mesh_Data");

	// Mesh_Connectivity
	struct Multiarray_Vector_i* v_to_v  = (struct Multiarray_Vector_i*) mesh->mesh_conn->v_to_v;
	struct Multiarray_Vector_i* v_to_lf = (struct Multiarray_Vector_i*) mesh->mesh_conn->v_to_lf;

	if ((diff_Multiarray_Vector_i(v_to_v,mesh_test_data->v_to_v)   != 0) ||
	    (diff_Multiarray_Vector_i(v_to_lf,mesh_test_data->v_to_lf) != 0)) {
		test_print_failure(test_info,"Mesh Connectivity");
		pass = 0;
		print_diff_Multiarray_Vector_i(v_to_v,mesh_test_data->v_to_v);
		print_diff_Multiarray_Vector_i(v_to_lf,mesh_test_data->v_to_lf);
	}

	// Mesh_Vertices
	test_print_warning(test_info,"Not yet testing Mesh_Vertices");

	destructor_Mesh_Test_Data(mesh_test_data);
}

// Level 1 ********************************************************************************************************** //

/** \brief Construtor for \ref Mesh_Test_Data::v_to_v.
 *	\return Standard. */
struct Multiarray_Vector_i* constructor_v_to_v
	(const char*const mesh_name ///< The test mesh name.
	);

/** \brief Construtor for \ref Mesh_Test_Data::v_to_lf.
 *	\return Standard. */
struct Multiarray_Vector_i* constructor_v_to_lf
	(const char*const mesh_name ///< The test mesh name.
	);

static void set_Mesh_Input
	(struct Mesh_Input*const mesh_input, const int d, const int domain_type, const bool mesh_unrealistic,
	 const char* mesh_name, const char* geom_name_, const char* geom_spec_, const char* input_path_)
{
	const_cast_i(&mesh_input->d,d);
	const_cast_i(&mesh_input->domain_type,domain_type);
	const_cast_bool(&mesh_input->mesh_unrealistic,mesh_unrealistic);

	char* mesh_name_full = (char*) mesh_input->mesh_name_full;
	strcpy(mesh_name_full,"../testing/meshes/");
	strcat(mesh_name_full,mesh_name);

	char* geom_name = (char*) mesh_input->geom_name;
	strcpy(geom_name,geom_name_);

	char* geom_spec = (char*) mesh_input->geom_spec;
	strcpy(geom_spec,geom_spec_);
	const_cast_c1(&mesh_input->geom_spec,geom_spec);

	char* input_path = (char*) mesh_input->input_path;
	strcpy(input_path,input_path_);
}

static struct Mesh_Test_Data* constructor_Mesh_Test_Data (const char*const mesh_name)
{
	struct Mesh_Test_Data* mesh_test_data = malloc(sizeof *mesh_test_data); // free

	mesh_test_data->v_to_v  = constructor_v_to_v(mesh_name);  // destructed
	mesh_test_data->v_to_lf = constructor_v_to_lf(mesh_name); // destructed

	return mesh_test_data;
}

static void destructor_Mesh_Test_Data (struct Mesh_Test_Data* mesh_test_data)
{
	destructor_Multiarray_Vector_i(mesh_test_data->v_to_v);
	destructor_Multiarray_Vector_i(mesh_test_data->v_to_lf);

	free(mesh_test_data);
}

// Level 2 ********************************************************************************************************** //

struct Multiarray_Vector_i* constructor_v_to_v (const char*const mesh_name)
{
	if (strstr(mesh_name,"curved_2d_mixed.msh")) {
		const int n_el     = 12;
		const int ext_V[]  = {3,3,3,3,3,3,3,3,4,4,4,4};
		const int data_V[] =
			{  1, -1,  8,
			   4,  2,  0,
			   3, -1,  1,
			   6, -1,  2,
			   5,  1, 10,
			  -1,  6,  4,
			   7,  3,  5,
			  -1, -1,  6,
			   0,  9, -1, 10,
			   8, -1, -1, 11,
			   4, 11,  8, -1,
			  10, -1,  9, -1, };

		return constructor_copy_Multiarray_Vector_i_i(data_V,ext_V,1,n_el); // returned
	} else if (strstr(mesh_name,"straight_2d_quad_periodic.msh")) {
		const int n_el     = 4;
		const int ext_V[]  = {4,4,4,4};
		const int data_V[] =
			{  2,  2,  1,  1,
			   3,  3,  0,  0,
			   0,  0,  3,  3,
			   1,  1,  2,  2, };

		return constructor_copy_Multiarray_Vector_i_i(data_V,ext_V,1,n_el); // returned
	}

	printf("%s\n",mesh_name);
	EXIT_UNSUPPORTED;
}

struct Multiarray_Vector_i* constructor_v_to_lf (const char*const mesh_name)
{
	if (strstr(mesh_name,"curved_2d_mixed.msh")) {
		const int n_el     = 12;
		const int ext_V[]  = {3,3,3,3,3,3,3,3,4,4,4,4};
		const int data_V[] =
			{      2, 20002,     0,
			       1,     2,     0,
			       2, 20002,     1,
			       1, 10001,     0,
			       2,     0,     0,
			   30002,     2,     0,
			       2,     0,     1,
			   30002, 10001,     0,
			       2,     0, 20002,     2,
			       1, 10001, 20002,     2,
			       2,     0,     3, 30002,
			       1, 10001,     3, 30002, };

		return constructor_copy_Multiarray_Vector_i_i(data_V,ext_V,1,n_el); // returned
	} else if (strstr(mesh_name,"straight_2d_quad_periodic.msh")) {
		const int n_el     = 4;
		const int ext_V[]  = {4,4,4,4};
		const int data_V[] =
			{  1,  0,  3,  2,
			   1,  0,  3,  2,
			   1,  0,  3,  2,
			   1,  0,  3,  2, };

		return constructor_copy_Multiarray_Vector_i_i(data_V,ext_V,1,n_el); // returned
	}

	EXIT_UNSUPPORTED;
}
