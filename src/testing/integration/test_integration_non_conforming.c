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

#include "test_base.h"
#include "test_integration.h"

#include "element_geometry.h"
#include "element_solver.h"
#include "face_solver.h"
#include "volume_solver.h"

#include "multiarray.h"

#include "intrusive.h"
#include "multiarray_operator.h"
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
 *  This checks:
 *  - That the values of the geometry nodes on non-conforming faces match exactly.
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
EXIT_ADD_SUPPORT; // Ensure that this is working with Bezier geometry as well.

	structor_simulation(&sim,'d',adapt_type,p,ml,p_prev,ml_prev,NULL,'r');

	destructor_Integration_Test_Info(int_test_info);

	PetscFinalize();
	OUTPUT_SUCCESS;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Get the pointer to the appropriate \ref Geometry_Element::cv0_vgc_fgc operator.
 *  \return See brief. */
const struct Operator* get_operator__cv0_vgc_fgc
	(const int side_index,                 ///< The index of the side of the face under consideration.
	 const struct Solver_Face*const s_face ///< The current \ref Solver_Face.
	);

static void check_face_geometry (const struct Simulation*const sim)
{
	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
		const struct Face*const face = (struct Face*) curr;
		if (face->boundary)
			continue;

		const struct Solver_Volume*const s_vol[] = { (struct Solver_Volume*) face->neigh_info[0].volume,
		                                             (struct Solver_Volume*) face->neigh_info[1].volume, };

print_const_Multiarray_d(s_vol[0]->geom_coef);
print_const_Multiarray_d(s_vol[1]->geom_coef);
		const struct Solver_Face*const s_face = (struct Solver_Face*) curr;
		const struct Operator* cv0_vgc_fgc_op = get_operator__cv0_vgc_fgc(0,s_face);
UNUSED(cv0_vgc_fgc_op);

	}
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		const struct Solver_Volume*const s_vol = (struct Solver_Volume*) curr;
struct Volume* vol = (struct Volume*) curr;
printf("%d %d %d\n",vol->index,s_vol->ml,s_vol->p_ref);
//print_const_Multiarray_d(vol->xyz_ve);
	}
EXIT_ADD_SUPPORT;
}

// Level 1 ********************************************************************************************************** //

const struct Operator* get_operator__cv0_vgc_fgc (const int side_index, const struct Solver_Face*const s_face)
{
	const struct Face*const face            = (struct Face*) s_face;
	const struct Volume*const vol           = face->neigh_info[side_index].volume;
	const struct Solver_Volume*const s_vol  = (struct Solver_Volume*) vol;
	const struct Geometry_Element*const g_e = &((struct Solver_Element*)vol->element)->g_e;

	const int ind_lf   = face->neigh_info[side_index].ind_lf,
	          ind_href = face->neigh_info[side_index].ind_href;
	const int p_v = s_vol->p_ref,
	          p_f = s_face->p_ref;

	const int curved = ( (s_face->cub_type == 's') ? 0 : 1 );
	assert(curved);

	return get_Multiarray_Operator(g_e->cv0_vgc_fgc,(ptrdiff_t[]){ind_lf,ind_href,0,p_f,p_v});
}
