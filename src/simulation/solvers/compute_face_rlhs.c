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

#include "compute_face_rlhs.h"

#include "face_solver.h"
#include "volume_solver.h"
#include "element_solver.h"

#include "matrix.h"
#include "multiarray.h"
#include "vector.h"

#include "boundary.h"
#include "flux.h"
#include "multiarray_operator.h"
#include "numerical_flux.h"
#include "operator.h"
#include "simulation.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/// \brief Version of \ref scale_by_Jacobian_fptr_T scaling \ref Numerical_Flux_T::nnf.
static void scale_by_Jacobian_nnf
	(struct Numerical_Flux*const num_flux, ///< See brief.
	 const struct Solver_Face*const s_face ///< See brief.
		);

/// \brief Version of \ref scale_by_Jacobian_fptr_T scaling \ref Numerical_Flux_T::Neigh_Info_NF_T::dnnf_ds.
static void scale_by_Jacobian_dnnf_ds
	(struct Numerical_Flux*const num_flux, ///< See brief.
	 const struct Solver_Face*const s_face ///< See brief.
		);

/// \brief Version of \ref scale_by_Jacobian_fptr_T scaling \ref Numerical_Flux_T::Neigh_Info_NF_T::dnnf_dg.
static void scale_by_Jacobian_dnnf_dg
	(struct Numerical_Flux*const num_flux, ///< See brief.
	 const struct Solver_Face*const s_face ///< See brief.
		);

// Interface functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "compute_face_rlhs_T.c"
#include "undef_templates_type.h"

#include "def_templates_type_dc.h"
#include "compute_face_rlhs_T.c"
#include "undef_templates_type.h"

void scale_by_Jacobian_i1 (struct Numerical_Flux*const num_flux, const struct Solver_Face*const s_face)
{
	scale_by_Jacobian_nnf(num_flux,s_face);
	scale_by_Jacobian_dnnf_ds(num_flux,s_face);
}

void scale_by_Jacobian_i2 (struct Numerical_Flux*const num_flux, const struct Solver_Face*const s_face)
{
	scale_by_Jacobian_nnf(num_flux,s_face);
	scale_by_Jacobian_dnnf_dg(num_flux,s_face);
}

void scale_by_Jacobian_i12 (struct Numerical_Flux*const num_flux, const struct Solver_Face*const s_face)
{
	scale_by_Jacobian_nnf(num_flux,s_face);
	scale_by_Jacobian_dnnf_ds(num_flux,s_face);
	scale_by_Jacobian_dnnf_dg(num_flux,s_face);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void scale_by_Jacobian_nnf (struct Numerical_Flux*const num_flux, const struct Solver_Face*const s_face)
{
	const struct const_Vector_d jacobian_det_fc = interpret_const_Multiarray_as_Vector_d(s_face->jacobian_det_fc);
	scale_Multiarray_by_Vector_d('L',1.0,(struct Multiarray_d*)num_flux->nnf,&jacobian_det_fc,false);
}

static void scale_by_Jacobian_dnnf_ds (struct Numerical_Flux*const num_flux, const struct Solver_Face*const s_face)
{
	const struct Face*const face = (struct Face*) s_face;
	assert((!face->boundary && num_flux->neigh_info[1].dnnf_ds != NULL) ||
	       ( face->boundary && num_flux->neigh_info[1].dnnf_ds == NULL));

	const struct const_Vector_d jacobian_det_fc = interpret_const_Multiarray_as_Vector_d(s_face->jacobian_det_fc);
	scale_Multiarray_by_Vector_d('L',1.0,(struct Multiarray_d*)num_flux->neigh_info[0].dnnf_ds,&jacobian_det_fc,false);
	if (!face->boundary)
		scale_Multiarray_by_Vector_d('L',1.0,(struct Multiarray_d*)num_flux->neigh_info[1].dnnf_ds,&jacobian_det_fc,false);
}

static void scale_by_Jacobian_dnnf_dg (struct Numerical_Flux*const num_flux, const struct Solver_Face*const s_face)
{
	const struct Face*const face = (struct Face*) s_face;
	assert((!face->boundary && num_flux->neigh_info[1].dnnf_dg != NULL) ||
	       ( face->boundary && num_flux->neigh_info[1].dnnf_dg == NULL));

	const struct const_Vector_d jacobian_det_fc = interpret_const_Multiarray_as_Vector_d(s_face->jacobian_det_fc);
	scale_Multiarray_by_Vector_d('L',1.0,(struct Multiarray_d*)num_flux->neigh_info[0].dnnf_dg,&jacobian_det_fc,false);
	if (!face->boundary)
		scale_Multiarray_by_Vector_d('L',1.0,(struct Multiarray_d*)num_flux->neigh_info[1].dnnf_dg,&jacobian_det_fc,false);
}
