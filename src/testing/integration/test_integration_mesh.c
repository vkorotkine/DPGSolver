// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/**	\file
 */

#include "test_integration_mesh.h"

#include <string.h>

#include "Macros.h"
#include "intrusive.h"
#include "element.h"
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

// Interface functions ********************************************************************************************** //

void test_integration_mesh (const char*const mesh_name)
{
	struct Mesh_Input* mesh_input         = constructor_Mesh_Input(mesh_name);       // destructed
	struct const_Intrusive_List* elements = constructor_Element_List(mesh_input->d); // destructed
	struct Mesh* mesh                     = constructor_Mesh(mesh_input,elements);   // destructed

EXIT_ERROR("Add comparison functions here");

	destructor_Mesh(mesh);
	destructor_Elements((struct Intrusive_List*)elements);
	destructor_Mesh_Input(mesh_input);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

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

// Level 1 ********************************************************************************************************** //

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
