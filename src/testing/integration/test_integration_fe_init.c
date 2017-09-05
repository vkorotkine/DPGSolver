// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/**	\file
 */

#include "test_integration_fe_init.h"

#include <string.h>

#include "test_support_intrusive.h"
#include "test_support_matrix.h"

#include "Macros.h"

#include "simulation.h"
#include "intrusive.h"
#include "element.h"
#include "mesh.h"
#include "volume.h"
#include "face.h"

#include "file_processing.h"
#include "constants_mesh.h"

// Static function declarations ************************************************************************************* //

/** \brief Compare members of the \ref Volume and \ref Face finite elements with expected values.
 *	\return 1 if tests passed. */
static bool compare_members_FE
	(struct Test_Info*const test_info, ///< Defined in \ref test_integration_mesh.
	 const struct Simulation*const sim ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

void test_integration_fe_init (struct Test_Info*const test_info, const char*const ctrl_name)
{
	struct Simulation*const sim = constructor_Simulation(ctrl_name); // destructed

	struct Mesh_Input mesh_input = set_Mesh_Input(sim);
	struct Mesh* mesh = constructor_Mesh(&mesh_input,sim->elements);

	sim->volumes = constructor_Volume_List(sim,mesh);
	sim->faces   = constructor_Face_List(sim,mesh);

	destructor_Mesh(mesh);

	const bool pass = compare_members_FE(test_info,sim);

	char* test_name_end = extract_name(ctrl_name,false); // free

	char test_name[STRLEN_MAX];
	strcpy(test_name,"FE initialization - ");
	strcat(test_name,test_name_end);

	free(test_name_end);

	test_increment_and_print(test_info,pass,test_name);

	destructor_Simulation(sim);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// Container for finite element test related data.
struct FE_Test_Data {
	struct Intrusive_List* volumes; ///< Defined in \ref Simulation.
	struct Intrusive_List* faces;   ///< Defined in \ref Simulation.
};

/** \brief Constructor for the name of the data file containing the finite element information.
 *	\return Standard. */
static char* constructor_data_name_fe
	(const char*const ctrl_name ///< The name of the control file.
	);

/**	\brief Constructor for the \ref FE_Test_Data.
 *	\return Standard. */
static struct FE_Test_Data* constructor_FE_Test_Data
	(const char*const data_name,                      ///< The name of the data file.
	 const struct const_Intrusive_List*const elements ///< \ref Simulation::elements.
	);

/// \brief Destructor for the \ref FE_Test_Data.
static void destructor_FE_Test_Data
	(struct FE_Test_Data* fe_test_data ///< \ref FE_Test_Data.
	);

static bool compare_members_FE
	(struct Test_Info*const test_info, const struct Simulation*const sim)
{
	bool pass = 1;

	char* data_name = constructor_data_name_fe(sim->ctrl_name_full); // destructed (free)

	struct FE_Test_Data* fe_test_data = constructor_FE_Test_Data(data_name,sim->elements); // destructed
	free(data_name);

	// Volumes
	for (const struct Intrusive_Link* curr = sim->volumes->first, *curr_test = fe_test_data->volumes->first;
		 curr || curr_test;
		 curr = curr->next, curr_test = curr_test->next)
	{
		if (!(curr && curr_test)) {
			pass = 0;
			test_print_failure(test_info,"Volumes (different number)");
			break;
		}

		struct Volume* volume      = (struct Volume*) curr,
		             * volume_test = (struct Volume*) curr_test;

		if ((volume->index != volume_test->index)                                       ||
		    (volume->boundary != volume_test->boundary)                                 ||
		    (volume->curved != volume_test->curved)                                     ||
		    (diff_const_Matrix_d(volume->xyz_ve,volume_test->xyz_ve,NODETOL_MESH) != 0) ||
		    (volume->element->type != volume_test->element->type))
		{
			test_print_failure(test_info,"Volume");
			pass = 0;
			printf("index: %d %d\n",volume->index,volume_test->index);
			printf("boundary: %d %d\n",volume->boundary,volume_test->boundary);
			printf("curved: %d %d\n",volume->curved,volume_test->curved);
			print_diff_const_Matrix_d(volume->xyz_ve,volume_test->xyz_ve,NODETOL_MESH);
			printf("elem_type: %d %d\n",volume->element->type,volume_test->element->type);
		}
	}

	destructor_FE_Test_Data(fe_test_data);

	return pass;
}

// Level 1 ********************************************************************************************************** //

static char* constructor_data_name_fe (const char*const ctrl_name)
{
	char* data_name_mid = extract_name(ctrl_name,true); // free

	char*const data_name = malloc(STRLEN_MAX * sizeof *data_name); // returned

	strcpy(data_name,"../testing/fe/");
	strcat(data_name,data_name_mid);
	strcat(data_name,".fe.data");

	free(data_name_mid);

	return data_name;
}

static struct FE_Test_Data* constructor_FE_Test_Data
	(const char*const data_name, const struct const_Intrusive_List*const elements)
{
	struct FE_Test_Data* fe_test_data = malloc(sizeof *fe_test_data); // returned

	fe_test_data->volumes = constructor_file_name_IL("Volume",data_name,elements); // destructed
//	fe_test_data->faces   = constructor_file_name_IL("Face",data_name,elements);   // destructed

	return fe_test_data;
}

static void destructor_FE_Test_Data (struct FE_Test_Data* fe_test_data)
{
	destructor_Volumes(fe_test_data->volumes);
//	destructor_Faces(fe_test_data->faces);

	free(fe_test_data);
}
