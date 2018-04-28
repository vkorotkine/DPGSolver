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
#include <stdio.h>
#include <string.h>

#include "macros.h"
#include "definitions_intrusive.h"

#include "def_templates_restart.h"

#include "def_templates_matrix.h"
#include "def_templates_multiarray.h"

#include "def_templates_volume_solver.h"

#include "def_templates_solution.h"

// Static function declarations ************************************************************************************* //

/** \brief Return a \ref Multiarray_T\* container holding the solution values at the input coordinates.
 *  \return See brief. */
static struct Multiarray_T* constructor_sol_restart
	(const struct const_Multiarray_R* xyz, ///< xyz coordinates at which to evaluate the solution.
	 const struct Simulation* sim          ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

void set_sol_restart_T (const struct Simulation*const sim, struct Solution_Container_T sol_cont)
{
	sol_cont.using_restart = true;
	const struct const_Multiarray_R* xyz = constructor_xyz_sol_T(sim,&sol_cont); // destructed
	struct Multiarray_T* sol = constructor_sol_restart(xyz,sim); // destructed
	destructor_const_Multiarray_R(xyz);
	sol_cont.using_restart = false;

	update_Solution_Container_sol_T(&sol_cont,sol,sim);
	destructor_Multiarray_T(sol);
}

const struct const_Multiarray_T* constructor_const_sol_restart_T
	(const struct const_Multiarray_R*const xyz, const struct Simulation*const sim)
{
	EXIT_ADD_SUPPORT; // Ensure that the straight high-order coordinates are being passed here.
	struct Multiarray_T* sol = constructor_sol_restart(xyz,sim); // returned
	return (const struct const_Multiarray_T*) sol;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Container for information relating to the solution from which to restart.
 *
 *  The background portion of \ref Restart_Info::sss stores the centroids of the volumes from the restart mesh.
 */
struct Restart_Info {
	struct Simulation* sim; ///< A minimalist \ref Simulation set up from the restart file.

	const ptrdiff_t n_v;                    ///< The 'n'umber of 'v'olumes.
	struct Solver_Volume_T** solver_volume; ///< The array of solver volumes.

	struct SSS_ANN* sss; ///< \ref SSS_ANN.
};

/** \brief Get the pointer to a \ref Restart_Info with members set from the restart file.
 *  \return See brief.
 *
 *  \note The container is constructed the first time that this function is entered and then simply returned for every
 *        subsequent call to this function. This results in an expected, fixed size memory leak.
 */
static struct Restart_Info get_Restart_Info
	(const struct Simulation*const sim ///< Standard.
	);

static struct Multiarray_T* constructor_sol_restart
	(const struct const_Multiarray_R* xyz, const struct Simulation* sim)
{
	struct Restart_Info ri = get_Restart_Info(sim);

/// \todo Make external
	const struct const_Matrix_R xyz_M = interpret_const_Multiarray_as_Matrix_R(xyz);

	const bool requires_transpose = ( xyz_M.layout == 'R' ? false : true );
	if (requires_transpose)
		transpose_Matrix_R((struct Matrix_R*)&xyz_M,true);

	// sorted input nodes.
print_const_Matrix_R(&xyz_M);
	const struct Nodes_Sorted_ANN* ns = constructor_Nodes_Sorted_ANN(&xyz_M);

	// find nearest neighbors in background mesh
	struct Input_ANN ann_i = { .nodes_s = ns->nodes, };
	constructor_SSS_ANN_s(&ann_i,ri.sss); // destructed

	const struct const_Vector_i*const ind_ann = constructor_ann_indices_from_sss(ri.sss); // destructed
print_const_Vector_i(ind_ann);
	destructor_const_Vector_i(ind_ann);

	destructor_SSS_ANN_s(ri.sss); // destructed

	if (requires_transpose)
		transpose_Matrix_R((struct Matrix_R*)&xyz_M,true);

print_const_Matrix_R(ns->nodes);
print_const_Vector_i(ns->indices);

/* 1. Make array of background volume pointers to avoid linked list iteration.
 * 2.
 *
 */

EXIT_UNSUPPORTED;
	destructor_SSS_ANN_s(ri.sss);



}

// Level 1 ********************************************************************************************************** //

/// \brief Set the list of volumes in \ref Restart_Info.
static void constructor_volume_list
	(struct Restart_Info*const ri ///< Standard.
	);

/// \brief Set the required volume members from the restart file.
static void initialize_volumes_restart
	(const struct Simulation*const sim ///< Standard.
	);

/// \brief Initialize the background node portion of \ref Restart_Info::sss.
static void initialize_ann_background
	(struct Restart_Info*const ri ///< Standard.
	);

static struct Restart_Info get_Restart_Info (const struct Simulation*const sim)
{
	static bool needs_computation = true;
	static struct Restart_Info ri;
	if (needs_computation) {
		needs_computation = false;
		ri.sim = constructor_Simulation_restart(sim); // leaked (static)
		constructor_derived_computational_elements(ri.sim,IL_SOLVER); // leaked (static)

		constructor_volume_list(&ri);
		initialize_volumes_restart(ri.sim);
		initialize_ann_background(&ri);

		assert(ri.sim->elements->name == IL_ELEMENT);
		destructor_const_Elements(ri.sim->elements);
	}
	return ri;
}

// Level 2 ********************************************************************************************************** //

/// \brief Set the values of \ref Solver_Volume_T::sol_coef for all volumes from the input restart file.
static void initialize_volumes_sol_coef
	(FILE* file,                       ///< Pointer to the restart file.
	 const struct Simulation*const sim ///< Standard.
	);

static void constructor_volume_list (struct Restart_Info*const ri)
{
	const ptrdiff_t n_v = compute_n_volumes(ri->sim);

	struct Solver_Volume_T**const s_vols = malloc((size_t) n_v * sizeof *s_vols); // keep

	int n = 0;
	for (struct Intrusive_Link* curr = ri->sim->volumes->first; curr; curr = curr->next)
		s_vols[n++] = (struct Solver_Volume_T*) curr;

	ri->solver_volume = s_vols;
}

static void initialize_volumes_restart (const struct Simulation*const sim)
{
	FILE* file = fopen_checked(get_restart_name()); // closed

	char line[STRLEN_MAX];
	while (fgets(line,sizeof(line),file)) {
		if (strstr(line,"$SolutionCoefficients"))
			initialize_volumes_sol_coef(file,sim);
	}
	fclose(file);
}

static void initialize_ann_background (struct Restart_Info*const ri)
{
	const struct const_Matrix_d*const xyz_min_max = constructor_volume_xyz_min_max(ri->sim); // destructed
	const struct const_Matrix_d*const centroids   = constructor_volume_centroids(ri->sim);   // destructed
print_const_Matrix_d(centroids);

	ri->sss = calloc(1,sizeof *ri->sss); // keep

	struct Input_ANN ann_i = { .nodes_s = NULL, };

	ann_i.nodes_b = xyz_min_max;
	constructor_SSS_ANN_xyz(&ann_i,ri->sss); // keep
	destructor_const_Matrix_d(xyz_min_max);

	ann_i.nodes_b = centroids;
	constructor_SSS_ANN_b(&ann_i,ri->sss); // keep
	destructor_const_Matrix_d(centroids);
}

// Level 3 ********************************************************************************************************** //

#define LINELEN_MAX (10*10*10*5*25) ///< Line length which has enough space for p9, 5 variable, HEX solution data.

/// \brief Read the \ref Solver_Volume_T::sol_coef from the current line (assumed to be stored in the bezier basis).
static void read_sol_coef_bezier
	(struct Solver_Volume_T*const s_vol, ///< The current volume.
	 char* line                          ///< The current line of the file.
	);

static void initialize_volumes_sol_coef (FILE* file, const struct Simulation*const sim)
{
	char line[LINELEN_MAX];
	char* endptr = NULL;

	fgets_checked(line,sizeof(line),file);
	ptrdiff_t n_v = strtol(line,&endptr,10);

	struct Intrusive_Link* curr = sim->volumes->first;

	ptrdiff_t row = 0;
	while (fgets(line,sizeof(line),file) != NULL) {
		if (strstr(line,"$EndSolutionCoefficients"))
			break;
		if (row++ == n_v)
			EXIT_ERROR("Reading too many rows.\n");
		if (curr == NULL)
			EXIT_ERROR("Reading too many volumes.\n");

		read_sol_coef_bezier((struct Solver_Volume_T*)curr,line);
		curr = curr->next;
	}
}

// Level 4 ********************************************************************************************************** //

static void read_sol_coef_bezier (struct Solver_Volume_T*const s_vol, char* line)
{
	enum { n_i = 3, order = 2, };

	discard_line_values(&line,1);

	int data_i[n_i] = { 0, };
	read_line_values_i(&line,n_i,data_i,false);

	const_cast_i(&s_vol->p_ref,data_i[0]);

	const ptrdiff_t exts[] = { data_i[1], data_i[2], };
	struct Multiarray_T*const s_coef = s_vol->sol_coef;
	resize_Multiarray_T(s_coef,order,exts);

#if TYPE_RC == TYPE_REAL
	read_line_values_d(&line,compute_size(order,exts),s_coef->data);
#else
	EXIT_ADD_SUPPORT;
#endif
}
