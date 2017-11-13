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

#include <assert.h>

#include "face_solver.h"
#include "element_solver.h"
#include "volume_solver.h"

#include "matrix.h"
#include "multiarray.h"
#include "vector.h"

#include "multiarray_operator.h"
#include "operator.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

const struct Operator* get_operator__tw0_vt_fc (const int side_index, const struct Solver_Face* s_face)
{
	const struct Face* face           = (struct Face*) s_face;
	const struct Volume* vol          = face->neigh_info[side_index].volume;
	const struct Solver_Volume* s_vol = (struct Solver_Volume*) vol;
	const struct Solver_Element* e    = (struct Solver_Element*) vol->element;

	const int ind_lf   = face->neigh_info[side_index].ind_lf,
	          ind_href = face->neigh_info[side_index].ind_href;
	const int p_v = s_vol->p_ref,
	          p_f = s_face->p_ref;

	const int curved = ( (s_face->cub_type == 's') ? 0 : 1 );
	return get_Multiarray_Operator(e->tw0_vt_fc[curved],(ptrdiff_t[]){ind_lf,ind_href,0,p_f,p_v});
}

void permute_Matrix_d_fc
	(struct Matrix_d* data, const char perm_layout, const int side_index_dest, const struct Solver_Face* s_face)
{
	assert(perm_layout == 'R');
	assert(data->layout == 'R');

	const struct const_Vector_i* nc_fc = get_operator__nc_fc(side_index_dest,s_face);
	permute_Matrix_d_V(data,nc_fc);
}

const struct const_Vector_i* get_operator__nc_fc (const int side_index_dest, const struct Solver_Face* s_face)
{
	const struct Neigh_Info* neigh_info = &((struct Face*)s_face)->neigh_info[side_index_dest];

	struct Volume* vol = neigh_info->volume;
	const struct Solver_Element* e = (struct Solver_Element*) vol->element;

	const int ind_ord = neigh_info->ind_ord,
	          ind_e   = get_face_element_index((struct Face*)s_face),
	          p_f     = s_face->p_ref;
	const int curved = ( (s_face->cub_type == 's') ? 0 : 1 );

	return get_const_Multiarray_Vector_i(e->nc_fc[curved],(ptrdiff_t[]){ind_ord,ind_e,ind_e,0,0,p_f,p_f});
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
