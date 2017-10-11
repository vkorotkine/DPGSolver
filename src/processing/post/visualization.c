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

#include "visualization.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "macros.h"
#include "definitions_alloc.h"
#include "definitions_intrusive.h"
#include "definitions_visualization.h"

#include "volume.h"
#include "solver_volume.h"

#include "multiarray.h"
#include "vector.h"

#include "computational_elements.h"
#include "element_plotting.h"
#include "file_processing.h"
#include "multiarray_operator.h"
#include "nodes_plotting.h"
#include "operator.h"
#include "simulation.h"

// Static function declarations ************************************************************************************* //

///\{ \name Available visualization softwares.
#define VIS_SOFTWARE_PARAVIEW 11 ///< Paraview.
///\}

/// \brief Output the visualization of the specified output in a format suitable for Paraview.
static void output_visualization_paraview
	(const struct Simulation* sim, ///< \ref Simulation.
	 const int vis_type            ///< The type of visualization. Options: see \ref definitions_visualization.h.
	);

// Interface functions ********************************************************************************************** //

void output_visualization (struct Simulation* sim, const int vis_type)
{
	assert(sim->volumes->name == IL_SOLVER_VOLUME);
	assert(sim->faces->name   == IL_SOLVER_FACE);

	constructor_derived_Elements(sim,IL_PLOTTING_ELEMENT);

	output_visualization_paraview(sim,vis_type);

/// \todo change input here to desired base list.
	destructor_derived_Elements(sim,IL_PLOTTING_ELEMENT);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Output the visualization of the specified computational element geometry in vtk xml format.
static void output_visualization_vtk_geom
	(const char geom_type,        ///< The type of geometry feature to output. Options: 'v'olumes, 'e'dges.
	 const struct Simulation* sim ///< \ref Simulation.
	);

static void output_visualization_paraview (const struct Simulation* sim, const int vis_type)
{
	switch (vis_type) {
	case VIS_GEOM_VOLUMES: // fallthrough
	case VIS_GEOM_EDGES:
		output_visualization_vtk_geom('v',sim);
		output_visualization_vtk_geom('e',sim);
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",vis_type);
		break;
	}
}

// Level 1 ********************************************************************************************************** //

/** \brief Set the output file name specific to the visualization under consideration.
 *  \return The output name (no free needed). */
const char* set_output_name
	(const int vis_software,    /**< The software to be used for the visualization.
	                             *   Options: see definitions_visualization.h. */
	 const char*const name_spec ///< The specific name.
	);

/// \brief Print the input string to the file with the specified number of tabs.
static void fprintf_tn
	(FILE* file,               ///< The file.
	 const int n_tabs,         ///< The number of tabs
	 const char*const string_i ///< The input string.
	);

/// \brief Print the vtk header for the output file.
static void fprint_vtk_header_footer
	(FILE* file,             ///< The file.
	 const bool is_parallel, ///< Flag for whether the output file is the parallel or normal file.
	 const char hf_type      ///< Type indicator for whether the 'h'eader or 'f'ooter should be printed.
	);

/// \brief Print the 's'tart/'e'nd of a piece to the file.
static void fprint_vtk_piece
	(FILE* file,                                ///< The file.
	 const char sp_type,                        ///< Type indicator for 's'erial or 'p'arallel.
	 const char se_type,                        ///< Type indicator for 's'tart or 'e'nd.
	 const char geom_type,                      ///< Defined for \ref output_visualization_vtk_geom.
	 const struct const_Multiarray_d* xyz,      ///< Coordinates of the points (required for 's'erial, 's'tart).
	 const struct const_Plotting_Nodes* p_nodes ///< \ref Plotting_Nodes (required for 's'erial, 's'tart).
	);

static void output_visualization_vtk_geom (const char geom_type, const struct Simulation* sim)
{
	assert(sim->elements->name == IL_PLOTTING_ELEMENT);

	static char output_part[STRLEN_MIN] = { 0, };
	sprintf(output_part,"%s%c","geom_",geom_type);

	const char*const output_name = set_output_name(VIS_SOFTWARE_PARAVIEW,output_part);

	if (sim->mpi_rank == 0) {
		static char parallel_name[STRLEN_MAX] = { 0, };
		sprintf(parallel_name,"%s%s",output_name,".pvtu");

		FILE* p_file = fopen_create_dir(parallel_name);

		fprint_vtk_header_footer(p_file,true,'h');

		fprint_vtk_piece(p_file,'p','s',geom_type,NULL,NULL);

		for (int i = 0; i < sim->mpi_size; ++i)
			fprintf(p_file,"<Piece Source=\"%s_%d.vtu\"/>\n",output_part,i);
		fprintf(p_file,"\n");

		fprint_vtk_header_footer(p_file,true,'f');

		fclose(p_file);
	}

	static char serial_name[STRLEN_MAX] = { 0, };
	sprintf(serial_name,"%s%s%d%s",output_name,"_",sim->mpi_rank,".vtu");

	FILE* s_file = fopen_create_dir(serial_name);

	fprint_vtk_header_footer(s_file,false,'h');

	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Volume* base_volume   = (struct Volume*)curr;
		struct Solver_Volume* volume = (struct Solver_Volume*)curr;

		const struct const_Plotting_Element* element = (const struct const_Plotting_Element*)base_volume->element;

		const int p = volume->p_ref;
		const struct Operator* cv0_vg_vp =
			(!base_volume->curved ? get_Multiarray_Operator(element->cv0_vgs_vps,(ptrdiff_t[]){0,0,p,1})
			                      : get_Multiarray_Operator(element->cv0_vgc_vpc,(ptrdiff_t[]){0,0,p,p}) );

		print_const_Multiarray_d(volume->geom_coef);

		const struct const_Multiarray_d* g_coef = volume->geom_coef;
		const struct const_Multiarray_d* xyz_p =
			constructor_mm_NN1_Operator_const_Multiarray_d(cv0_vg_vp,g_coef,'R','d',g_coef->order,NULL); // destructed

		fprint_vtk_piece(s_file,'s','s',geom_type,xyz_p,element->p_nodes[p]);
		fprint_vtk_piece(s_file,'s','e',geom_type,NULL,NULL);

		destructor_const_Multiarray_d(xyz_p);
	}

	fprint_vtk_header_footer(s_file,false,'f');

	fclose(s_file);
}

// Level 2 ********************************************************************************************************** //

/** \brief Print a \ref const_Multiarray_d to a file with the input number of tabs before each row, and padding rows
 *         with zeroes until they have 3 entries. */
void fprint_const_Multiarray_d_vtk_point
	(FILE* file,                        ///< The file.
	 const int n_tab,                   ///< The number of tabs.
	 const struct const_Multiarray_d* a ///< Standard.
	);

/** \brief Print the array of cummulative sum of `ext_0` of the \ref Vector_i\*s of the input \ref
 * const_Multiarray_Vector_i to a file with the input number of tabs before each row. */
void fprint_const_Multiarray_Vector_i_offsets
	(FILE* file,                               ///< The file.
	 const int n_tab,                          ///< The number of tabs.
	 const struct const_Multiarray_Vector_i* a ///< Standard.
	);

const char* set_output_name (const int vis_software, const char*const name_spec)
{
	static char output_name[STRLEN_MAX] = { 0, };

	strcpy(output_name,"../output/");
	switch (vis_software) {
	case VIS_SOFTWARE_PARAVIEW:
		strcat(output_name,"paraview/");
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",vis_software);
		break;
	}

	if (strcmp(name_spec,"geom_v") == 0 ||
	    strcmp(name_spec,"geom_e") == 0) {
		strcat(output_name,name_spec);
	} else {
		EXIT_ERROR("Unsupported: %s\n",name_spec);
	}

	return output_name;
}

static void fprintf_tn (FILE* file, const int n_tabs, const char*const string_i)
{
	for (int i = 0; i < n_tabs; i++)
		fprintf(file,"\t");
	fprintf(file,"%s\n",string_i);
}

static void fprint_vtk_header_footer (FILE* file, const bool is_parallel, const char hf_type)
{
	static char vtk_type[STRLEN_MIN] = { 0, },
	            string_i[STRLEN_MAX] = { 0, };

	if (is_parallel)
		strcpy(vtk_type,"PUnstructuredGrid");
	else
		strcpy(vtk_type,"UnstructuredGrid");


	if (hf_type == 'h') {
		fprintf_tn(file,0,"<?xml version=\"1.0\"?>");

		sprintf(string_i,"%s%s%s","<VTKFile type=\"",vtk_type,"\" version=\"0.1\" byte_order=\"LittleEndian\">");
		fprintf_tn(file,0,string_i);

		if (is_parallel)
			sprintf(string_i,"%s%s%s","<",vtk_type," GhostLevel=\"0\">");
		else
			sprintf(string_i,"%s%s%s","<",vtk_type,">");
		fprintf_tn(file,0,string_i);
	} else if (hf_type == 'f') {
		sprintf(string_i,"%s%s%s","</",vtk_type,">");
		fprintf_tn(file,0,string_i);
		fprintf_tn(file,0,"</VTKFile>");
	} else {
		EXIT_ERROR("Unsupported: %c\n",hf_type);
	}
}

static void fprint_vtk_piece
	(FILE* file, const char sp_type, const char se_type, const char geom_type, const struct const_Multiarray_d* xyz,
	const struct const_Plotting_Nodes* p_nodes)
{
	// Note: Points **must** have 3 values.

	assert(sp_type == 's' || sp_type == 'p');
	assert(se_type == 's' || se_type == 'e');
	assert(geom_type == 'v' || geom_type == 'e');

	if (sp_type == 'p') {
		if (se_type == 's') {
			fprintf_tn(file,1,"<PPoints>");
				fprintf(file,"\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\"/>\n");
			fprintf_tn(file,1,"</PPoints>");

			fprintf_tn(file,1,"<PCells>");
				fprintf_tn(file,2,"<PDataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\"/>");
				fprintf_tn(file,2,"<PDataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\"/>");
				fprintf_tn(file,2,"<PDataArray type=\"UInt8\" Name=\"types\" format=\"ascii\"/>");
			fprintf_tn(file,1,"</PCells>");
			fprintf(file,"\n");
		} else if (se_type == 'e') {
			EXIT_UNSUPPORTED;
		}
	} else if (sp_type == 's') {
		if (se_type == 's') {
			const struct const_Multiarray_Vector_i* connect = NULL;
			const struct const_Vector_i* vtk_types          = NULL;
			if (geom_type == 'v') {
				connect = p_nodes->connect;
				vtk_types = p_nodes->vtk_types;
			} else if (geom_type == 'e') {
				connect = p_nodes->connect_e;
				vtk_types = p_nodes->vtk_types_e;
			}

			fprintf(file,"\n<Piece NumberOfPoints=\"%td\" NumberOfCells=\"%td\">\n",
			        xyz->extents[0],connect->extents[0]);
			fprintf_tn(file,1,"<Points>");
				fprintf(file,"\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n");
				fprint_const_Multiarray_d_vtk_point(file,2,xyz);
				fprintf_tn(file,2,"</DataArray>");
			fprintf_tn(file,1,"</Points>");

			fprintf_tn(file,1,"<Cells>");
				fprintf_tn(file,2,"<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">");
				fprint_const_Multiarray_Vector_i(file,2,connect);
				fprintf_tn(file,2,"</DataArray>");

				fprintf_tn(file,2,"<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">");
				fprint_const_Multiarray_Vector_i_offsets(file,2,connect);
				fprintf_tn(file,2,"</DataArray>");

				fprintf_tn(file,2,"<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">");
				fprint_const_Vector_i(file,2,vtk_types);
				fprintf_tn(file,2,"</DataArray>");
			fprintf_tn(file,1,"</Cells>");
		} else if (se_type == 'e') {
			fprintf_tn(file,0,"</Piece>\n");
		}
	}
}

// Level 3 ********************************************************************************************************** //

void fprint_const_Multiarray_d_vtk_point (FILE* file, const int n_tab, const struct const_Multiarray_d* a)
{
	const int order               = a->order;
	const ptrdiff_t*const extents = a->extents;

	assert(order == 2);

	const bool transpose_Ma = ( a->layout == 'R' ? false : true );
	if (transpose_Ma)
		transpose_Multiarray_d((struct Multiarray_d*)a,true);

	const ptrdiff_t ext_0 = extents[0],
	                ext_1 = extents[1];

	for (ptrdiff_t i = 0; i < ext_0; ++i) {
		const double* data = get_row_const_Multiarray_d(i,a);

		for (int j = 0; j < n_tab; ++j)
			fprintf(file,"\t");
		for (ptrdiff_t j = 0; j < 3; ++j) {
			if (j < ext_1)
				fprintf(file," % .8e",data[j]);
			else
				fprintf(file," %d",0);
		}
		fprintf(file,"\n");
	}

	if (transpose_Ma)
		transpose_Multiarray_d((struct Multiarray_d*)a,true);
}

void fprint_const_Multiarray_Vector_i_offsets (FILE* file, const int n_tab, const struct const_Multiarray_Vector_i* a)
{
	const ptrdiff_t size = compute_size(a->order,a->extents);

	struct Vector_i* a_V = constructor_empty_Vector_i(size); // destructed
	ptrdiff_t sum = 0;
	for (ptrdiff_t i = 0; i < size; ++i) {
		a_V->data[i] = sum + a->data[i]->ext_0;
		sum = a_V->data[i];
	}
	fprint_Vector_i(file,n_tab,a_V);
	destructor_Vector_i(a_V);
}
