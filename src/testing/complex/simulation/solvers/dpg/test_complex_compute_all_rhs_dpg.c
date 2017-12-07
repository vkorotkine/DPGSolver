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

#include "test_complex_compute_all_rhs_dpg.h"

#include "complex_matrix.h"
#include "complex_multiarray.h"
#include "complex_vector.h"

#include "test_complex_flux.h"
#include "test_complex_compute_face_rhs.h"
#include "test_complex_compute_volume_rhs.h"
#include "test_complex_numerical_flux.h"
#include "test_complex_operators.h"
#include "test_complex_test_case.h"

#include "test_complex_face_solver_dpg.h"
#include "test_complex_volume_solver_dpg.h"


#include <assert.h>
#include "petscmat.h"

#include "macros.h"
#include "definitions_test_integration.h"

#include "element_solver_dpg.h"

#include "matrix.h"
#include "multiarray.h"
#include "vector.h"

#include "compute_all_rlhs_dpg.h"
#include "compute_face_rlhs.h"
#include "compute_volume_rlhs.h"
#include "intrusive.h"
#include "multiarray_operator.h"
#include "operator.h"
#include "simulation.h"
#include "solve.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/** \brief Version of \ref add_to_petsc_Mat_Vec_dpg setting a single column of Solver_Storage_Implicit::A to the values
 *         computed using the complex step rhs. */
static void add_to_petsc_Mat_dpg_c
	(const struct Solver_Volume_c* s_vol,  ///< See brief.
	 const struct const_Vector_c* rhs_neg, ///< See brief.
	 struct Solver_Storage_Implicit* ssi   ///< See brief.
	);

// Interface functions ********************************************************************************************** //

#include "def_templates_type_dc.h"
#include "compute_all_rlhs_dpg_T.c"

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void add_to_petsc_Mat_dpg_c
	(const struct Solver_Volume_c* s_vol, const struct const_Vector_c* rhs_neg, struct Solver_Storage_Implicit* ssi)
{
	const ptrdiff_t ext_0 = rhs_neg->ext_0;

	const struct const_Vector_i* idxm = constructor_petsc_idxm_dpg_c(ext_0,s_vol); // destructed.

	PetscScalar rhs_c_data[ext_0];
	for (int i = 0; i < ext_0; ++i)
		rhs_c_data[i] = cimag((-rhs_neg->data[i])/CX_STEP);

	MatSetValues(ssi->A,(PetscInt)ext_0,idxm->data,1,&ssi->col,rhs_c_data,ADD_VALUES);

	destructor_const_Vector_i(idxm);
}
