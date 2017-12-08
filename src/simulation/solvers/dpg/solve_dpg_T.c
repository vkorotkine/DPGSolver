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

#include "solve_dpg_T.h"

#include <assert.h>

#include "macros.h"


#include "def_templates_solve_dpg.h"

#include "def_templates_face_solver.h"
#include "def_templates_volume_solver.h"

#include "def_templates_multiarray.h"

#include "def_templates_test_case.h"

// Static function declarations ************************************************************************************* //

/// \brief Increment the appropriate rows of nnz based on the off-diagonal terms.
static void increment_nnz_off_diag
	(struct Vector_i* nnz,              ///< Vector of 'n'umber of 'n'on-'z'ero entries per row.
	 const struct Solver_Face_T* s_face ///< The current face.
	);

// Interface functions ********************************************************************************************** //

void update_ind_dof_dpg_T (const struct Simulation* sim)
{
	ptrdiff_t dof = 0;
	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
		struct Solver_Face_T* s_face = (struct Solver_Face_T*) curr;

		struct Multiarray_T* nf_coef = s_face->nf_coef;
		const ptrdiff_t size = compute_size(nf_coef->order,nf_coef->extents);
		if (size == 0) {
			const_cast_ptrdiff(&s_face->ind_dof,-1);
			continue;
		}

		const_cast_ptrdiff(&s_face->ind_dof,dof);
		dof += size;
	}

	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Solver_Volume_T* s_vol = (struct Solver_Volume_T*) curr;

		const_cast_ptrdiff(&s_vol->ind_dof,dof);

		struct Multiarray_T* sol_coef = s_vol->sol_coef;
		dof += compute_size(sol_coef->order,sol_coef->extents);
	}
	assert(dof == compute_dof(sim));
}

struct Vector_i* constructor_nnz_dpg_T (const struct Simulation* sim)
{
	struct Test_Case_T* test_case = (struct Test_Case_T*) sim->test_case_rc->tc;
	assert(test_case->has_2nd_order == false); // Add support.

	const ptrdiff_t dof = compute_dof(sim);
	struct Vector_i* nnz = constructor_zero_Vector_i(dof); // returned

	// Volume contribution (Diagonal)
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Solver_Volume_T* s_vol = (struct Solver_Volume_T*) curr;

		struct Multiarray_T* sol_coef = s_vol->sol_coef;
		const ptrdiff_t size = compute_size(sol_coef->order,sol_coef->extents);
		increment_nnz(nnz,s_vol->ind_dof,size,size);
	}

	// Face contributions (Diagonal and Off-diagonal)
	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
		struct Face* face = (struct Face*) curr;
		if (face->boundary)
			continue;

		struct Solver_Face_T* s_face = (struct Solver_Face_T*) curr;

		// Diagonal
		struct Multiarray_T* nf_coef = s_face->nf_coef;
		const ptrdiff_t size_nf = compute_size(nf_coef->order,nf_coef->extents);
		increment_nnz(nnz,s_face->ind_dof,size_nf,size_nf);

		// Off-diagonal
		increment_nnz_off_diag(nnz,s_face);
	}
	return nnz;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Increment the appropriate rows of nnz based on the volume off-diagonal terms.
static void increment_nnz_off_diag_v
	(struct Vector_i* nnz,              ///< Defined for \ref increment_nnz_off_diag.
	 const int ind_neigh,               ///< The index of the neighbouring volume.
	 const ptrdiff_t n_rows,            ///< The number of rows to increment.
	 const struct Solver_Face_T* s_face ///< The current face.
	);

/// \brief Increment the appropriate rows of nnz based on the face off-diagonal terms.
static void increment_nnz_off_diag_f
	(struct Vector_i* nnz,              ///< Defined for \ref increment_nnz_off_diag.
	 const int ind_neigh,               ///< The index of the neighbouring volume.
	 const ptrdiff_t n_rows,            ///< The number of rows to increment.
	 const struct Solver_Face_T* s_face ///< The current face.
	);

static void increment_nnz_off_diag (struct Vector_i* nnz, const struct Solver_Face_T* s_face)
{
	const struct Face* face = (struct Face*) s_face;
	struct Multiarray_T* nf_coef = s_face->nf_coef;
	const ptrdiff_t size_nf = compute_size(nf_coef->order,nf_coef->extents);

	for (int n = 0; n < 2; ++n) {
		if (n == 1 && face->boundary)
			continue;

		increment_nnz_off_diag_v(nnz,n,size_nf,s_face);
		increment_nnz_off_diag_f(nnz,n,size_nf,s_face);
	}
}

// Level 1 ********************************************************************************************************** //

static void increment_nnz_off_diag_v
	(struct Vector_i* nnz, const int ind_neigh, const ptrdiff_t n_rows, const struct Solver_Face_T* s_face)
{
	const struct Face* face             = (struct Face*) s_face;
	const struct Solver_Volume_T* s_vol = (struct Solver_Volume_T*) face->neigh_info[ind_neigh].volume;
	struct Multiarray_T* sol_coef = s_vol->sol_coef;
	ptrdiff_t size_sol = compute_size(sol_coef->order,sol_coef->extents);

	increment_nnz(nnz,s_face->ind_dof,n_rows,size_sol);
	increment_nnz(nnz,s_vol->ind_dof,size_sol,n_rows);
}

static void increment_nnz_off_diag_f
	(struct Vector_i* nnz, const int ind_neigh, const ptrdiff_t n_rows, const struct Solver_Face_T* s_face)
{
	const struct Face* face  = (struct Face*) s_face;
	const struct Volume* vol = (struct Volume*) face->neigh_info[ind_neigh].volume;

	for (int i = 0; i < NFMAX;    ++i) {
	for (int j = 0; j < NSUBFMAX; ++j) {
		const struct Face* face_n = vol->faces[i][j];
		if (!face_n || face_n->boundary)
			continue;

		if (face == face_n)
			continue;

		const struct Solver_Face_T* s_face_n = (struct Solver_Face_T*) face_n;

		struct Multiarray_T* nf_coef = s_face_n->nf_coef;
		const ptrdiff_t size_nf_n = compute_size(nf_coef->order,nf_coef->extents);
		increment_nnz(nnz,s_face->ind_dof,n_rows,size_nf_n);
	}}
}
