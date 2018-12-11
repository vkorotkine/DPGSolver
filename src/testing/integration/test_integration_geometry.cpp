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

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "petscsys.h"

#include "test_base.h"
#include "test_support_intrusive.h"
#include "test_support_multiarray.h"

#include "macros.h"
#include "definitions_geometry.h"
#include "definitions_intrusive.h"
#include "definitions_mesh.h"
#include "definitions_tol.h"
#include "definitions_visualization.h"

#include "computational_elements.h"
#include "face.h"
#include "volume.h"
#include "face_solver.h"

#include "file_processing.h"
#include "geometry.h"
#include "simulation.h"
#include "test_case.h"
#include "visualization.h"

// Static function declarations ************************************************************************************* //

/** \brief Compare \ref Solver_Volume_T::geom_coef and \ref Solver_Face_T::normals_fc finite members with their expected
 *         values.
 *  \return `true` if tests passed. */
static bool compare_members_geom
	(struct Test_Info*const test_info, ///< \ref Test_Info.
	 const struct Simulation*const sim ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

/** \test Performs integration testing for the geometry initialization (\ref test_integration_geometry.cpp).
 *  \return 0 on success.
 *
 *  Compares the members listed in \ref compare_members_geom with their expected values.
 */
int main
	(int argc,   ///< Standard.
	 char** argv ///< Standard.
	)
{
	PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);

	assert_condition_message(argc == 2,"Invalid number of input arguments");
	const char* ctrl_name = argv[1];

	struct Test_Info test_info = { .n_warn = 0, };

	struct Simulation*const sim = constructor_Simulation(ctrl_name); // destructed

	constructor_derived_computational_elements(sim,IL_SOLVER); // destructed
	constructor_derived_Elements(sim,IL_ELEMENT_SOLVER);       // destructed

	set_up_solver_geometry(sim);

	output_visualization(sim,VIS_GEOM_VOLUMES);
	output_visualization(sim,VIS_NORMALS);

	const bool pass = compare_members_geom(&test_info,sim);

	destructor_derived_computational_elements(sim,IL_BASE);
	destructor_derived_Elements(sim,IL_ELEMENT);

	test_print_warning(&test_info,"Not testing free-stream preservation.");

	destructor_Simulation(sim);

	assert_condition(pass);
	output_warning_count(&test_info);

	PetscFinalize();
	OUTPUT_SUCCESS;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// Container for finite element test related data.
struct Geom_Test_Data {
	struct Intrusive_List* faces; ///< Defined in \ref Simulation.
};

/** \brief Constructor for the \ref Geom_Test_Data.
 *  \return Standard. */
static struct Geom_Test_Data* constructor_Geom_Test_Data
	(const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Destructor for the \ref Geom_Test_Data.
static void destructor_Geom_Test_Data
	(struct Geom_Test_Data* geom_test_data ///< Standard.
	);

static bool compare_members_geom (struct Test_Info*const test_info, const struct Simulation*const sim)
{
	/// \todo Add comparison of Solver_Volume::geom_coef.
	UNUSED(test_info);
	bool pass = true;

	if (sim->domain_type == DOM_BLENDED) {
		struct Test_Case* test_case = (struct Test_Case*)sim->test_case_rc->tc;
		assert(test_case->geom_parametrization == GEOM_PRM_RADIAL_PROJ); // May be made flexible in future
	}

	struct Geom_Test_Data* geom_test_data = constructor_Geom_Test_Data(sim); // destructed

	assert(sim->faces->name == IL_FACE_SOLVER);

	for (const struct Intrusive_Link* curr = sim->faces->first, *curr_test = geom_test_data->faces->first;
	     curr || curr_test;
	     curr = curr->next, curr_test = curr_test->next)
	{
		if (!(curr && curr_test)) {
			pass = false;
			expect_condition(pass,"Faces (different number)");
			break;
		}

		struct Solver_Face* face      = (struct Solver_Face*) curr,
		                  * face_test = (struct Solver_Face*) curr_test;

		const double tol = 1e2*EPS;
		if (diff_const_Multiarray_d(face->normals_fc,face_test->normals_fc,tol) != 0) {
			pass = false;
			expect_condition(pass,"Face");

			print_diff_const_Multiarray_d(face->normals_fc,face_test->normals_fc,tol);
		}
	}

	destructor_Geom_Test_Data(geom_test_data);

	return pass;
}

// Level 1 ********************************************************************************************************** //

static struct Geom_Test_Data* constructor_Geom_Test_Data (const struct Simulation* sim)
{
	struct Geom_Test_Data* geom_test_data = calloc(1,sizeof *geom_test_data); // returned

	const char* data_name = set_data_file_name_integration(sim->ctrl_name_full,"fe");
	struct Intrusive_List* b_volumes = constructor_file_name_IL("Volume",data_name,sim->elements,NULL);    // destructed
	struct Intrusive_List* b_faces   = constructor_file_name_IL("Face",data_name,sim->elements,b_volumes); // destructed

	const size_t sizeof_base    = sizeof(struct Face),
	             sizeof_derived = sizeof(struct Solver_Face);

	struct Intrusive_List* faces = constructor_empty_IL(IL_FACE_SOLVER,b_faces); // moved
	for (struct Intrusive_Link* curr = b_faces->first; curr; curr = curr->next)
		push_back_IL(faces,constructor_derived_Intrusive_Link(curr,sizeof_base,sizeof_derived));

	for (struct Intrusive_Link* curr = faces->first; curr; curr = curr->next)
		constructor_derived_Solver_Face((struct Face*)curr,sim);

	data_name = set_data_file_name_integration(sim->ctrl_name_full,"geom");
	constructor_file_name_derived_Faces(faces,data_name);

	destructor_IL_base(faces);
	destructor_Volumes(b_volumes);

	geom_test_data->faces = faces;

	return geom_test_data;
}

static void destructor_Geom_Test_Data (struct Geom_Test_Data* geom_test_data)
{
	struct Intrusive_List* faces = geom_test_data->faces;
	for (struct Intrusive_Link* curr = faces->first; curr; curr = curr->next)
		destructor_derived_Solver_Face((struct Face*)curr);
	destructor_Faces(faces);

	free(geom_test_data);
}
