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
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "macros.h"
#include "definitions_intrusive.h"
#include "definitions_h_ref.h"
#include "definitions_tol.h"

#include "def_templates_restart.h"

#include "def_templates_matrix.h"
#include "def_templates_multiarray.h"

#include "def_templates_face_solver.h"
#include "def_templates_volume_solver.h"

#include "def_templates_geometry.h"
#include "def_templates_math_functions.h"
#include "def_templates_solution.h"
#include "def_templates_test_case.h"

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

/** \brief Constructor for a \ref Vector_T of indices of the background volumes which contain the input nodes.
 *  \return See brief. */
static const struct const_Vector_i* constructor_indices_vol_background
	(const struct Restart_Info*const ri,       ///< Standard.
	 const struct Nodes_Sorted_ANN*const ns,   ///< Standard.
	 const struct const_Vector_i*const ind_ann ///< Background volume Indices of the approximate nearest neighbors.
	);

/** \brief Constructor for the \ref Multiarray_T\* of solution data from the stored solution of the background mesh.
 *  \return See brief. */
static struct Multiarray_T* constructor_sol_restart_from_background
	(const struct Restart_Info*const ri,          ///< Standard.
	 const struct Nodes_Sorted_ANN*const ns,      ///< Standard.
	 const struct const_Vector_i*const ind_vol_b, ///< Indices of the background volumes containing the nodes.
	 const struct Simulation*const sim            ///< Standard.
	);

static struct Multiarray_T* constructor_sol_restart
	(const struct const_Multiarray_R* xyz, const struct Simulation* sim)
{
	struct Restart_Info ri = get_Restart_Info(sim);

	const struct const_Matrix_R xyz_M = interpret_const_Multiarray_as_Matrix_R(xyz);
	const struct Nodes_Sorted_ANN*const ns = constructor_Nodes_Sorted_ANN_with_trans((struct Matrix_R*)&xyz_M); // d

	struct Input_ANN ann_i = { .nodes_s = ns->nodes, };
	constructor_SSS_ANN_s(&ann_i,ri.sss); // destructed

	const struct const_Vector_i*const ind_ann   = constructor_ann_indices_from_sss(ri.sss);           // destructed
	const struct const_Vector_i*const ind_vol_b = constructor_indices_vol_background(&ri,ns,ind_ann); // destructed
	destructor_const_Vector_i(ind_ann);

	struct Multiarray_T*const sol = constructor_sol_restart_from_background(&ri,ns,ind_vol_b,sim); // returned
	destructor_Nodes_Sorted_ANN(ns);
	destructor_SSS_ANN_s(ri.sss);
	destructor_const_Vector_i(ind_vol_b);

	return sol;
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

/** \brief If required, update the index of the background volume to that of the neighbouring volume which is closest
 *         to the search node. */
static void update_ind_vol_b
	(int*const ind_vol_b,                ///< Pointer to memory location holding the index of the background volume.
	 const struct Restart_Info*const ri, ///< Standard.
	 const Real*const data_node          ///< Pointer to the xyz coordinate data of the search node.
	);

static struct Restart_Info get_Restart_Info (const struct Simulation*const sim)
{
	static bool needs_computation = true;
	static struct Restart_Info ri;
	if (needs_computation) {
		printf("\tSetting up restart.\n");

		needs_computation = false;
		ri.sim = constructor_Simulation_restart(sim); // leaked (static)
		constructor_derived_Elements(ri.sim,IL_ELEMENT_SOLVER); // destructed
		constructor_derived_computational_elements(ri.sim,IL_SOLVER); // leaked (static)

		set_up_solver_geometry_T(ri.sim);
		set_up_solver_geometry_p1_T(ri.sim);
		initialize_volumes_restart(ri.sim);
		initialize_ann_background(&ri);
		constructor_volume_list(&ri);

		destructor_derived_Elements(ri.sim,IL_ELEMENT);
//		destructor_const_Elements(ri.sim->elements); // leaked (static)
		printf("\tFinished Setting up restart.\n\n");
	}
	return ri;
}

static const struct const_Vector_i* constructor_indices_vol_background
	(const struct Restart_Info*const ri, const struct Nodes_Sorted_ANN*const ns,
	 const struct const_Vector_i*const ind_ann)
{
	const struct const_Matrix_R*const nodes = ns->nodes;
	assert(nodes->layout == 'R');

	const ptrdiff_t n_n = ns->nodes->ext_0;

	struct Vector_i*const ind_vol_b = (struct Vector_i*) constructor_copy_const_Vector_i(ind_ann); // returned
	for (int n = 0; n < n_n; ++n) {
		const Real*const data_node = get_row_const_Matrix_R(n,nodes);
		update_ind_vol_b(&(ind_vol_b->data[n]),ri,data_node);
	}
	return (struct const_Vector_i*) ind_vol_b;
}

static struct Multiarray_T* constructor_sol_restart_from_background
	(const struct Restart_Info*const ri, const struct Nodes_Sorted_ANN*const ns,
	 const struct const_Vector_i*const ind_vol_b, const struct Simulation*const sim)
{
	const struct const_Matrix_R*const nodes = ns->nodes;
	assert(nodes->layout == 'R');

	const ptrdiff_t n_n = nodes->ext_0;
	assert(DIM == nodes->ext_1);

	struct Test_Case_T*const test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
	const int n_var = test_case->n_var;

	struct Multiarray_T* sol = constructor_empty_Multiarray_T('R',2,(ptrdiff_t[]){n_n,n_var}); // returned
	struct Matrix_T sol_row = { .layout = 'R', .ext_0 = 1, .ext_1 = n_var, .data = NULL, };

	// Note: Changes are required below if ext_0 != 1 as it is assumed that the matrix is interchangeably
	//       interpretted as being either row or column major oriented.
	struct Matrix_R xyz_m = { .layout = 'R', .ext_0 = 1, .ext_1 = DIM, .owns_data = false, .data = NULL, };
	const struct const_Matrix_R*const xyz = (struct const_Matrix_R*) &xyz_m;

	for (int n = 0; n < n_n; ++n) {
		const struct Solver_Volume_T* s_vol_b = ri->solver_volume[ind_vol_b->data[n]];
		const struct Volume*const vol_b = (struct Volume*) s_vol_b;

		const int e_type = vol_b->element->type;
		const struct const_Matrix_R xyz_ve_M = interpret_const_Multiarray_as_Matrix_R(vol_b->xyz_ve);
		xyz_m.data = (Real*) get_row_const_Matrix_R(n,nodes);

		struct Matrix_R*const rst = constructor_inverse_mapping_mutable(e_type,&xyz_ve_M,xyz); // destructed
		rst->layout = 'C';

		const int s_type = vol_b->element->s_type;
		constructor_basis_fptr constructor_basis = get_constructor_basis_bezier_by_super_type(s_type);

		const int p_b = s_vol_b->p_ref;
		const struct const_Matrix_d*const cv_rst = constructor_basis(p_b,(struct const_Matrix_R*)rst); // dest.
		destructor_Matrix_R(rst);

		sol_row.data = get_row_Multiarray_T(ns->indices->data[n],sol);
		struct Matrix_T sol_coef_b = interpret_Multiarray_as_Matrix_T(s_vol_b->sol_coef);

		mm_RTT('N','N',1.0,0.0,cv_rst,(struct const_Matrix_T*)&sol_coef_b,&sol_row);
		destructor_const_Matrix_R(cv_rst);
	}

	transpose_Multiarray_T(sol,true);
	assert(sol->layout == 'C');
	return sol;
}

// Level 2 ********************************************************************************************************** //

/// \brief Set the values of \ref Solver_Volume_T::sol_coef for all volumes from the input restart file.
static void initialize_volumes_sol_coef
	(FILE* file,                       ///< Pointer to the restart file.
	 const struct Simulation*const sim ///< Standard.
	);

/** \brief Return a pointer to a statically allocated array holding the vector of the difference between the input node
 *         coordinate data and the xyz vertex data for a vertex lying on the face under consideration.
 *  \return See brief. */
static const Real* get_xyz_node_m_ve
	(const int e_type,                            ///< \ref Element::type.
	 const int ind_lf,                            ///< Index of the face under consideration.
	 const Real*const xyz_node,                   ///< xyz coordinate data for the search node.
	 const struct const_Multiarray_R*const xyz_ve ///< \ref Volume::xyz_ve.
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

	ri->sss = calloc(1,sizeof *ri->sss); // keep

	struct Input_ANN ann_i = { .nodes_s = NULL, };

	ann_i.nodes_b = xyz_min_max;
	constructor_SSS_ANN_xyz(&ann_i,ri->sss); // keep
	destructor_const_Matrix_d(xyz_min_max);

	ann_i.nodes_b = centroids;
	constructor_SSS_ANN_b(&ann_i,ri->sss); // keep
	destructor_const_Matrix_d(centroids);
}

static void update_ind_vol_b
	(int*const ind_vol_b, const struct Restart_Info*const ri, const Real*const data_node)
{
	bool found = false;
	while (!found) {
		const struct Solver_Volume_T* s_vol_b = ri->solver_volume[*ind_vol_b];
		const struct Volume*const vol_b = (struct Volume*) s_vol_b;

		const int e_type = vol_b->element->type;
		const struct const_Multiarray_R*const xyz_ve = vol_b->xyz_ve;

		found = true;
		Real max_n_dot_diff = 0.0;
		for (int i = 0; i < NFMAX;    ++i) {
		for (int j = 0; j < NSUBFMAX; ++j) {
			const struct Face*const face_b = vol_b->faces[i][j];
			if (face_b == NULL || face_b->boundary)
				continue;
			const struct Solver_Face_T*const s_face_b = (struct Solver_Face_T*) face_b;

			const int si = compute_side_index_face(face_b,vol_b);
			const double scale_n = ( si == 0 ? 1.0 : -1.0 );

			const Real*const data_n_p1      = s_face_b->normals_p1->data;
			const Real*const data_node_m_ve = get_xyz_node_m_ve(e_type,i,data_node,xyz_ve);
			const Real n_dot_diff = scale_n*dot_R(DIM,data_n_p1,data_node_m_ve);

			if (n_dot_diff > max_n_dot_diff) {
				*ind_vol_b = get_volume_neighbour(vol_b,face_b)->index;
				max_n_dot_diff = n_dot_diff;
				if (abs_R(max_n_dot_diff) > EPS)
					found = false;
			}
		}}
	}
}


// Level 3 ********************************************************************************************************** //

#define LINELEN_MAX (10*10*10*5*25) ///< Line length which has enough space for p9, 5 variable, HEX solution data.

/// \brief Read the \ref Solver_Volume_T::sol_coef from the current line (assumed to be stored in the bezier basis).
static void read_sol_coef_bezier
	(struct Solver_Volume_T*const s_vol, ///< The current volume.
	 char* line                          ///< The current line of the file.
	);

/** \brief Return the index of a vertex on the element face under consideration.
 *  \return See brief. */
static int get_ind_ve
	(const int e_type, ///< \ref Element::type.
	 const int ind_lf  ///< The index of the face.
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

static const Real* get_xyz_node_m_ve
	(const int e_type, const int ind_lf, const Real*const xyz_node, const struct const_Multiarray_R*const xyz_ve)
{
	static Real node_m_ve[DIM] = { 0.0, };

	const int ind_ve = get_ind_ve(e_type,ind_lf);

	const Real*const ve_data = get_row_const_Multiarray_R(ind_ve,xyz_ve);
	for (int d = 0; d < DIM; ++d)
		node_m_ve[d] = xyz_node[d]-ve_data[d];
	return node_m_ve;
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
	struct Multiarray_R*const s_coef_R = constructor_empty_Multiarray_R(s_coef->layout,order,exts); // destructed
	read_line_values_d(&line,compute_size(order,exts),s_coef_R->data);
	copy_into_Multiarray_T_from_R(s_coef,(struct const_Multiarray_R*)s_coef_R);
	destructor_Multiarray_R(s_coef_R);
#endif
}

static int get_ind_ve (const int e_type, const int ind_lf)
{
	switch (e_type) {
	case LINE:
		switch (ind_lf) {
			case H_LINE1_F0: return 0; break;
			case H_LINE1_F1: return 1; break;
			default: EXIT_ERROR("Unsupported: %d\n",ind_lf); break;
		}
		break;
	case TRI:
		switch (ind_lf) {
			case H_TRI1_F0: return 1; break;
			case H_TRI1_F1: return 0; break;
			case H_TRI1_F2: return 0; break;
			default: EXIT_ERROR("Unsupported: %d\n",ind_lf); break;
		}
		break;
	case QUAD:
		switch (ind_lf) {
			case H_QUAD1_F0: return 0; break;
			case H_QUAD1_F1: return 3; break;
			case H_QUAD1_F2: return 0; break;
			case H_QUAD1_F3: return 3; break;
			default: EXIT_ERROR("Unsupported: %d\n",ind_lf); break;
		}
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",e_type);
		break;
	}
}

#include "undef_templates_restart.h"

#include "undef_templates_matrix.h"
#include "undef_templates_multiarray.h"

#include "undef_templates_face_solver.h"
#include "undef_templates_volume_solver.h"

#include "undef_templates_geometry.h"
#include "undef_templates_math_functions.h"
#include "undef_templates_solution.h"
#include "undef_templates_test_case.h"
