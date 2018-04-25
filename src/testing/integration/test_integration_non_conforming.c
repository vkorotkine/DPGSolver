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
#include "geometry.h"
#include "intrusive.h"
#include "multiarray_operator.h"
#include "operator.h"
#include "simulation.h"

// Static function declarations ************************************************************************************* //

/** \brief Compare the values of the geometry related parameters on all faces.
 *
 *  The non-conformity correction is tested based on matching geometry nodes on all faces when interpolated from the
 *  volume on either side of the face. The currently stored normal vector is compared with the recomputed normal vector
 *  to ensure that conforming faces were updated if required.
 */
static void check_face_geometry
	(const struct Simulation*const sim ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

/** \test Performs integration testing for non-conforming meshes (both in h and in p)
 *        (\ref test_integration_non_conforming.c).
 *  \return 0 on success.
 *
 *  Note that meshes which do not pass this test resulted in the loss of optimal convergence for both the linear
 *  advection and Euler equations for isoparametric geometry representation (1/2 to 1 order). This condition thus
 *  appears to be necessary for optimal convergence in these cases.
 *
 *  Also note that the non-conforming treatment implemented here is different from that of Kopriva \cite Kopriva1996
 *  which uses mortar elements between volumes. As the current implementation does not require a basis on the face to
 *  represent the normal flux, it seems likely that there is a lower chance of aliasing. Given that all three of
 *  Kopriva's conditions are satisfied:
 *  1. Conservative treatment of numerical fluxes;
 *  2. Satisfaction of the outflow condition (Upwind values used to compute flux and unaffected by downwind values for
 *     cases where characteristics come only from the upwind direction);
 *  3. Same geometry representation used for both non-conforming elements on the face.
 *
 *  there does not seem to be any motivation to adopt the mortar method.
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
	structor_simulation(&sim,'c',adapt_type,p,ml,p_prev,ml_prev,ctrl_name_curr,'r',false); // destructed

	adapt_initial_mesh_if_required(sim);
	check_face_geometry(sim);

	structor_simulation(&sim,'d',ADAPT_0,p,ml,p_prev,ml_prev,NULL,'r',false);

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
	 const double tol,                                 ///< The tolerance for the difference check.
	 const char*const diff_name                        ///< Name to be displayed while printing differences.
	);

static void check_face_geometry (const struct Simulation*const sim)
{
	bool pass = true;
	const double tol = 5e3*EPS;
	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
		const struct Face*const face    = (struct Face*) curr;
		struct Solver_Face*const s_face = (struct Solver_Face*) curr;

		// Check normals
		const struct const_Multiarray_d* normals_fc[2] =
			{ constructor_copy_const_Multiarray_d(s_face->normals_fc), // destructed
			  NULL, };
		compute_geometry_face(s_face,sim);
		normals_fc[1] = s_face->normals_fc;

		if (diff_const_Multiarray_d(normals_fc[0],normals_fc[1],tol)) {
			pass = false;
			print_diff_geom(face,normals_fc,tol,"normals");
		}
		destructor_const_Multiarray_d(normals_fc[0]);

		// Check matching non-conforming geometry
		if (face->boundary || !face->curved)
			continue;

		const struct const_Multiarray_d* geom_fg[2] = {NULL,};
		for (int i = 0; i < 2; ++i)
			geom_fg[i] = constructor_geom_fg(i,0,s_face,false,false); // destructed

		if (diff_const_Multiarray_d(geom_fg[0],geom_fg[1],tol)) {
			pass = false;
			print_diff_geom(face,geom_fg,tol,"geom");
		}

		for (int i = 0; i < 2; ++i)
			destructor_const_Multiarray_d(geom_fg[i]);
	}
	assert_condition(pass);
}

// Level 1 ********************************************************************************************************** //

static void print_diff_geom
	(const struct Face*const face, const struct const_Multiarray_d*const geom_fg[2], const double tol,
	 const char*const diff_name)
{
	printf("\n\n\ndiff: %s\n\n",diff_name);
	for (int i = 0; i < 2; ++i) {
		const struct Volume*const vol          = face->neigh_info[i].volume;
		const struct Solver_Volume*const s_vol = (struct Solver_Volume*) vol;
		printf("volume %d: (ind: %d, p: %d, ml: %d).\n",i,vol->index,s_vol->p_ref,s_vol->ml);
	}
	print_const_Multiarray_d(geom_fg[0]);
	print_const_Multiarray_d(geom_fg[1]);
	print_diff_const_Multiarray_d(geom_fg[0],geom_fg[1],tol);
}
