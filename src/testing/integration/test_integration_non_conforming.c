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
#include "petscsys.h"

#include "macros.h"
#include "definitions_adaptation.h"
#include "definitions_tol.h"

#include "test_base.h"
#include "test_integration.h"
#include "test_support_multiarray.h"

#include "element_geometry.h"
#include "element_solver.h"
#include "face_solver.h"
#include "volume_solver.h"

#include "multiarray.h"

#include "adaptation.h"
#include "intrusive.h"
#include "multiarray_operator.h"
#include "operator.h"
#include "simulation.h"

// Static function declarations ************************************************************************************* //

/// \brief Compare the values of the geometry nodes on all faces.
static void check_face_geometry
	(const struct Simulation*const sim ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

/** \test Performs integration testing for non-conforming meshes (both in h and in p)
 *        (\ref test_integration_non_conforming.c).
 *  \return 0 on success.
 *
 *  For the implementation as it is at the current time, that a mesh with internally curved non-conforming faces passes
 *  this test or not does not have any impact on the recovery of optimal convergence orders for isoparametric geometry
 *  representation (March 11, 2018). It may potentially be interesting to investigate if the mortar method of Kopriva
 *  \cite Kopriva1996 resolves this issue.
 *
 *  This checks:
 *  - That the values of the geometry nodes on non-conforming faces match to machine precision.
 */
int main
	(int argc,   ///< Standard.
	 char** argv ///< Standard.
	)
{
	PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);

	assert_condition_message(argc == 2,"Invalid number of input arguments");

	const char* ctrl_name = argv[1];

	struct Integration_Test_Info* int_test_info = constructor_Integration_Test_Info(ctrl_name);

	const int p  = int_test_info->p_ref[0],
	          ml = int_test_info->ml[0],
	          p_prev  = p-1,
	          ml_prev = ml-1;

	const int adapt_type = int_test_info->adapt_type;
	const char*const ctrl_name_curr = set_file_name_curr(adapt_type,p,ml,false,ctrl_name);

	struct Simulation* sim = NULL;
	structor_simulation(&sim,'c',adapt_type,p,ml,p_prev,ml_prev,ctrl_name_curr,'r'); // destructed

	adapt_initial_mesh_if_required(sim);
	check_face_geometry(sim);

	structor_simulation(&sim,'d',ADAPT_0,p,ml,p_prev,ml_prev,NULL,'r');

	destructor_Integration_Test_Info(int_test_info);

	PetscFinalize();
	OUTPUT_SUCCESS;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Print information to terminal for non-matching geom values from either side of the face.
static void print_diff_geom
	(const struct Face*const face,                     ///< \ref Face.
	 const struct const_Multiarray_d*const geom_fg[2], ///< Face Geometry values interpolated from each volume.
	 const double tol                                  ///< The tolerance for the difference check.
	);

static void check_face_geometry (const struct Simulation*const sim)
{
	bool pass = true;
	const double tol = 4e3*EPS;
	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
		const struct Face*const face = (struct Face*) curr;
		if (face->boundary || !face->curved)
			continue;

		const struct Solver_Face*const s_face = (struct Solver_Face*) curr;

		const struct const_Multiarray_d* geom_fg[2] = {NULL,};
		for (int i = 0; i < 2; ++i)
			geom_fg[i] = constructor_geom_fg(i,0,s_face,false,false); // destructed

		if (diff_const_Multiarray_d(geom_fg[0],geom_fg[1],tol)) {
			pass = false;
			print_diff_geom(face,geom_fg,tol);
		}

		for (int i = 0; i < 2; ++i)
			destructor_const_Multiarray_d(geom_fg[i]);
	}
	assert_condition(pass);
}

// Level 1 ********************************************************************************************************** //

static void print_diff_geom
	(const struct Face*const face, const struct const_Multiarray_d*const geom_fg[2], const double tol)
{
	for (int i = 0; i < 2; ++i) {
		const struct Volume*const vol          = face->neigh_info[i].volume;
		const struct Solver_Volume*const s_vol = (struct Solver_Volume*) vol;
		printf("volume %d: (ind: %d, p: %d, ml: %d).\n",i,vol->index,s_vol->p_ref,s_vol->ml);
	}
	print_const_Multiarray_d(geom_fg[0]);
	print_const_Multiarray_d(geom_fg[1]);
	print_diff_const_Multiarray_d(geom_fg[0],geom_fg[1],tol);
}
