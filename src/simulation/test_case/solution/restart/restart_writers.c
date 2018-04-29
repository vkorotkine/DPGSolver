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
/// \file

#include "restart_writers.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include "macros.h"
#include "definitions_elements.h"
#include "definitions_mesh.h"
#include "definitions_tol.h"

#include "matrix.h"
#include "multiarray.h"
#include "vector.h"

#include "face.h"
#include "element.h"
#include "volume_solver.h"

#include "file_processing.h"
#include "intrusive.h"
#include "mesh_readers.h"
#include "simulation.h"
#include "solution.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

///\{ \name Supported restart file types.
#define RESTART_GMSH 1
///\}

/// \brief Write a restart file with the mesh component in gmsh format.
static void output_restart_gmsh
	(const struct Simulation*const sim ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

void output_restart (const struct Simulation*const sim)
{
	if (!outputting_restart()) {
		printf("Restart outputting is currently disabled in the test_case data file.\n");
		printf("Returning without outputting.\n");
		return;
	}

	switch (sim->domain_type) {
	case DOM_STRAIGHT:
	case DOM_PARAMETRIC:
		break; // Do nothing
	default:
		// Only affine restart files are currently supported.
/// \todo Think about whether blended domains could fit in this framework as only the boundary volumes are curved.
		EXIT_ERROR("Unsupported: %d\n",sim->domain_type);
		break;
	}

	if (strstr(sim->mesh_name_full,".msh"))
		output_restart_gmsh(sim);
	else
		EXIT_ERROR("Unsupported: %s\n",sim->mesh_name_full);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Container for information related to the restart mesh nodes.
struct Nodes_Info {
	const struct const_Matrix_d* nodes;       ///< The xyz node coordinates.
	const struct const_Vector_i* inds_sorted; ///< The indices of the sorted coordinates.
	const struct const_Vector_i* inds_unique; ///< The indices of the unique coordinates.
};

/// \brief Container for information related to the restart mesh elements.
struct Elements_Info {
	const struct const_Vector_i* v_types; ///< The volume element type indices.

	const struct const_Vector_i* f_types;              ///< The face element type indices.
	const struct const_Vector_i* f_bcs;                ///< The face boundary condition type indices.
	const struct const_Multiarray_Vector_i* ve_inds_f; ///< Vertex indices of the faces.
};

/** \brief Constructor for the \ref Nodes_Info container.
 *  \return See brief. */
static const struct Nodes_Info* constructor_Nodes_Info
	(const struct Simulation*const sim ///< \ref Simulation.
	);

/// \brief Destructor for the \ref Nodes_Info container.
static void destructor_Nodes_Info
	(const struct Nodes_Info*const nodes_info ///< Standard.
	);

/** \brief Constructor for the \ref Elements_Info container.
 *  \return See brief. */
static const struct Elements_Info* constructor_Elements_Info
	(const struct Simulation*const sim ///< \ref Simulation.
	);

/// \brief Destructor for the \ref Elements_Info container.
static void destructor_Elements_Info
	(const struct Elements_Info*const elements_info ///< Standard.
	);

/** \brief Return the name of the restart file.
 *  \return See brief. */
static const char* compute_restart_name
	(const int restart_type,           ///< The type of restart file. See options above.
	 const struct Simulation*const sim ///< \ref Simulation.
	);

/// \brief Print the format format to the gmsh file.
static void fprintf_format_gmsh
	(FILE* file ///< The file.
	);

/// \brief Print the nodes to the gmsh file.
static void fprintf_nodes_gmsh
	(FILE* file,                      ///< The file.
	 const struct Nodes_Info*const ni ///< Standard.
	);

/// \brief Print the elements to the gmsh file.
static void fprintf_elements_gmsh
	(FILE* file,                         ///< The file.
	 const struct Nodes_Info*const ni,   ///< Standard.
	 const struct Elements_Info*const ei ///< Standard.
	);

/** \brief Print the solution coefficients in the Bezier basis to the gmsh file.
 *
 *  The solution is outputted as Bezier basis coefficients (as opposed to nodal Lagrange basis values, for example) so
 *  that the operator required to obtain the values at the nodes of a non-nested mesh can be computed directly from the
 *  basis functions (without requiring the standard construction in \ref constructor_operators).
 */
static void fprintf_solution_bezier_gmsh
	(FILE* file,                       ///< The file.
	 const struct Simulation*const sim ///< Standard.
	);

static void output_restart_gmsh (const struct Simulation*const sim)
{
	const struct Nodes_Info*const ni    = constructor_Nodes_Info(sim);    // destructed
	const struct Elements_Info*const ei = constructor_Elements_Info(sim); // destructed

	const char*const f_name = compute_restart_name(RESTART_GMSH,sim);

	FILE* file = fopen_create_dir(f_name); // closed

	fprintf_format_gmsh(file);
	fprintf_nodes_gmsh(file,ni);
	fprintf_elements_gmsh(file,ni,ei);
	fprintf_solution_bezier_gmsh(file,sim);

	fclose(file);

	destructor_Nodes_Info(ni);
	destructor_Elements_Info(ei);
}

// Level 1 ********************************************************************************************************** //

static const struct Nodes_Info* constructor_Nodes_Info (const struct Simulation*const sim)
{
	const double eps_node = 1e-10;
	struct Nodes_Info* nodes_info = calloc(1,sizeof(*nodes_info)); // free

	struct Multiarray_d*const xyz_ve = constructor_empty_Multiarray_d('R',2,(ptrdiff_t[]){0,DIM}); // destructed
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		const struct Volume*const vol = (struct Volume*) curr;
		push_back_Multiarray_d(xyz_ve,vol->xyz_ve);
	}
	struct Matrix_d xyz_ve_tmp = interpret_Multiarray_as_Matrix_d(xyz_ve);
	nodes_info->inds_sorted = row_sort_DIM_Matrix_d(&xyz_ve_tmp,true,eps_node); // destructed
	nodes_info->inds_unique = make_unique_row_Multiarray_d(xyz_ve,eps_node,true); // destructed

	const ptrdiff_t* e = xyz_ve->extents;
	struct Matrix_d*const xyz_ve_M = constructor_move_Matrix_d_d('R',e[0],e[1],true,xyz_ve->data); // destructed
	xyz_ve->owns_data = false;
	nodes_info->nodes  = (struct const_Matrix_d*) xyz_ve_M;

	destructor_Multiarray_d(xyz_ve);
	return nodes_info;
}

static void destructor_Nodes_Info (const struct Nodes_Info*const nodes_info)
{
	destructor_const_Matrix_d(nodes_info->nodes);
	destructor_const_Vector_i(nodes_info->inds_sorted);
	destructor_const_Vector_i(nodes_info->inds_unique);
	free((void*)nodes_info);
}

static const struct Elements_Info* constructor_Elements_Info (const struct Simulation*const sim)
{
	struct Elements_Info* elements_info = calloc(1,sizeof(*elements_info)); // free

	struct Vector_i*const v_types = constructor_empty_Vector_i(0); // destructed
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		const struct Volume*const vol = (struct Volume*) curr;
		push_back_Vector_i(v_types,vol->element->type,false,false);
	}
	elements_info->v_types = (struct const_Vector_i*) v_types;

	struct Multiarray_Vector_i* ve_inds_v = constructor_empty_Multiarray_Vector_i(false,1,&v_types->ext_0); // dest.

	int ind_ve_v = 0;
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		const struct Volume*const vol = (struct Volume*) curr;

		const ptrdiff_t ext_0 = vol->xyz_ve->extents[0];
		struct Vector_i* ve_inds_lv = constructor_empty_Vector_i(ext_0); // moved
		for (int i = 0; i < ext_0; ++i)
			ve_inds_lv->data[i] = ind_ve_v++;
		ve_inds_v->data[vol->index] = ve_inds_lv;
	}

	struct Vector_i*const f_types = constructor_empty_Vector_i(0), // destructed
	               *const f_bcs   = constructor_empty_Vector_i(0); // destructed
	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
		const struct Face*const face = (struct Face*) curr;
		if (!face->boundary)
			continue;
		push_back_Vector_i(f_types,face->element->type,false,false);
		push_back_Vector_i(f_bcs,face->bc,false,false);
	}
	elements_info->f_types = (struct const_Vector_i*) f_types;
	elements_info->f_bcs   = (struct const_Vector_i*) f_bcs;

	struct Multiarray_Vector_i* ve_inds_f = constructor_empty_Multiarray_Vector_i(false,1,&f_types->ext_0); // dest.
	int ind_f = 0;
	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
		const struct Face*const face = (struct Face*) curr;
		if (!face->boundary)
			continue;
		const struct Volume*const vol = face->neigh_info[0].volume;
		const int ind_lf = face->neigh_info[0].ind_lf;

		const struct const_Vector_i*const f_ve_lf = vol->element->f_ve->data[ind_lf];
		const ptrdiff_t ext_0 = f_ve_lf->ext_0;
		struct Vector_i* ve_inds_lf = constructor_empty_Vector_i(ext_0); // moved
		for (int i = 0; i < ext_0; ++i)
			ve_inds_lf->data[i] = ve_inds_v->data[vol->index]->data[f_ve_lf->data[i]];
		ve_inds_f->data[ind_f] = ve_inds_lf;
		++ind_f;
	}
	destructor_Multiarray_Vector_i(ve_inds_v);
	elements_info->ve_inds_f = (const struct const_Multiarray_Vector_i*) ve_inds_f;

	return elements_info;
}

static void destructor_Elements_Info (const struct Elements_Info*const elements_info)
{
	destructor_const_Vector_i(elements_info->v_types);

	destructor_const_Vector_i(elements_info->f_types);
	destructor_const_Vector_i(elements_info->f_bcs);
	destructor_const_Multiarray_Vector_i(elements_info->ve_inds_f);
	free((void*)elements_info);
}

static const char* compute_restart_name
	(const int restart_type, const struct Simulation*const sim)
{
	static const char* name_path = "../restart/";

	char ml_p_spec[STRLEN_MIN];
	sprintf(ml_p_spec,"%s%d%s%d","__ml",sim->ml_p_curr[0],"__p",sim->ml_p_curr[1]);

	char name_ext[STRLEN_MIN];
	switch (restart_type) {
		case RESTART_GMSH: strcpy(name_ext,".msh");                      break;
		default:           EXIT_ERROR("Unsupported: %d\n",restart_type); break;
	}

	static char output_name[5*STRLEN_MAX];
	sprintf(output_name,"%s%s%c%s%c%s%s%s",
	        name_path,sim->pde_name,'/',sim->pde_spec,'/',"restart",ml_p_spec,name_ext);
	return output_name;
}

static void fprintf_format_gmsh (FILE* file)
{
	fprintf(file,"$MeshFormat\n");
	fprintf(file,"2.2 0 8\n");
	fprintf(file,"$EndMeshFormat\n");
}

static void fprintf_nodes_gmsh (FILE* file, const struct Nodes_Info*const ni)
{
	const ptrdiff_t n_n = ni->nodes->ext_0;
	fprintf(file,"$Nodes\n");
	fprintf(file,"%td\n",n_n);
	for (int n = 0; n < n_n; ++n) {
		const double*const data_n = get_row_const_Matrix_d(n,ni->nodes);
		fprintf(file,"%d",n+1);
		for (int d = 0; d < DMAX; ++d)
			fprintf(file," %g",( d < DIM ? data_n[d] : 0.0 ));
		fprintf(file,"\n");
	}
	fprintf(file,"$EndNodes\n");
}

static void fprintf_elements_gmsh (FILE* file, const struct Nodes_Info*const ni, const struct Elements_Info*const ei)
{
	const struct const_Vector_i*const inds_s = ni->inds_sorted;
	const struct const_Vector_i*const inds_u = ni->inds_unique;

	int inds_data[NVEMAX] = { 0, };
	struct Vector_i inds_ve = { .ext_0 = 0, .data = inds_data, };

	int ind_e = 0;
	// faces
	const ptrdiff_t n_ef = ei->f_types->ext_0,
	                n_ev = ei->v_types->ext_0;
	fprintf(file,"$Elements\n");
	fprintf(file,"%td\n",n_ef+n_ev);

	for (int e = 0; e < n_ef; ++e) {
		const int e_type = ei->f_types->data[e];
		fprintf(file,"%d %d",++ind_e,e_type);

		fprintf(file," %d %d",GMSH_N_TAGS,ei->f_bcs->data[e]);
		for (int i = 1; i < GMSH_N_TAGS; ++i)
			fprintf(file," %d",0);

		const struct const_Vector_i*const ve_inds_f = ei->ve_inds_f->data[e];
		const ptrdiff_t n_n = ve_inds_f->ext_0;
		for (int n = 0; n < n_n; ++n)
			inds_ve.data[n] = inds_u->data[inds_s->data[ve_inds_f->data[n]]]+1;
		inds_ve.ext_0 = n_n;
		reorder_nodes_gmsh(e_type,&inds_ve);
		for (int n = 0; n < n_n; ++n)
			fprintf(file," %d",inds_ve.data[n]);
		fprintf(file,"\n");
	}

	// volumes
	for (int e = 0, ind_n = 0; e < n_ev; ++e) {
		const int e_type = ei->v_types->data[e];
		fprintf(file,"%d %d",++ind_e,e_type);

		fprintf(file," %d",GMSH_N_TAGS);
		for (int i = 0; i < GMSH_N_TAGS; ++i)
			fprintf(file," %d",0);

		const int n_n = get_n_nodes(e_type);
		for (int n = 0; n < n_n; ++n) {
			inds_ve.data[n] = inds_u->data[inds_s->data[ind_n]]+1;
			++ind_n;
		}
		inds_ve.ext_0 = n_n;
		reorder_nodes_gmsh(e_type,&inds_ve);
		for (int n = 0; n < n_n; ++n)
			fprintf(file," %d",inds_ve.data[n]);
		fprintf(file,"\n");
	}
	fprintf(file,"$EndElements\n");
}

static void fprintf_solution_bezier_gmsh (FILE* file, const struct Simulation*const sim)
{
	enum { n_dec = 15 };
	const ptrdiff_t n_v = compute_n_volumes(sim);

	fprintf(file,"$SolutionCoefficients\n");
	fprintf(file,"%td\n",n_v);
	int v = 0;
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		const struct Solver_Volume*const s_vol = (struct Solver_Volume*) curr;
		struct Multiarray_d*const s_coef = constructor_s_coef_bezier(s_vol,sim);
		assert(s_coef->layout == 'C');

		const int p_ref = s_vol->p_ref;
		const ptrdiff_t n_dof = s_coef->extents[0],
		                n_var = s_coef->extents[1];

		fprintf(file,"%d %d %td %td",v+1,p_ref,n_dof,n_var);
		for (int i = 0; i < n_dof*n_var; ++i)
			fprintf(file," % .*e",n_dec,s_coef->data[i]);
		fprintf(file,"\n");

		++v;
	}
	fprintf(file,"$EndSolutionCoefficients\n");
}
