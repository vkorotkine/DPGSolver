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
#include "definitions_test_case.h"
#include "definitions_visualization.h"

#include "volume.h"
#include "volume_solver.h"
#include "face_solver.h"

#include "multiarray.h"
#include "vector.h"

#include "computational_elements.h"
#include "const_cast.h"
#include "element_plotting.h"
#include "file_processing.h"
#include "multiarray_operator.h"
#include "nodes_plotting.h"
#include "operator.h"
#include "simulation.h"
#include "solution_euler.h"
#include "test_case.h"

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
	assert(sim->volumes->name == IL_VOLUME_SOLVER);
	assert(sim->faces->name   == IL_FACE_SOLVER);

	constructor_derived_Elements(sim,IL_ELEMENT_PLOTTING);

	output_visualization_paraview(sim,vis_type);

	destructor_derived_Elements(sim,IL_ELEMENT);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Output the visualization of the specified computational element geometry in vtk xml format.
static void output_visualization_vtk_geom
	(const char geom_type,         ///< The type of geometry feature to output. Options: 'v'olumes, 'e'dges.
	 const struct Simulation* sim, ///< \ref Simulation.
	 const bool use_test_case_path ///< Flag for whether the output should be placed in the test case directory.
	);

/// \brief Output the visualization of the face normals in vtk xml format.
static void output_visualization_vtk_normals
	(const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Output the visualization of the computational element solution in vtk xml format.
static void output_visualization_vtk_sol
	(const struct Simulation* sim ///< \ref Simulation.
	);

static void output_visualization_paraview (const struct Simulation* sim, const int vis_type)
{
	switch (vis_type) {
	case VIS_GEOM_VOLUMES:
		output_visualization_vtk_geom('v',sim,false);
		output_visualization_vtk_geom('e',sim,false);
		break;
	case VIS_GEOM_EDGES:
		output_visualization_vtk_geom('e',sim,true);
		break;
	case VIS_NORMALS:
		output_visualization_vtk_normals(sim);
		break;
	case VIS_SOLUTION:
		output_visualization_vtk_sol(sim);
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
	(FILE* file,                    ///< The file.
	 const bool is_parallel,        ///< Flag for whether the output file is the parallel or normal file.
	 const char hf_type,            ///< Type indicator for whether the 'h'eader or 'f'ooter should be printed.
	 const char*const vtk_type_part ///< The partial vtk file type.
	);

/// \brief Print the 's'tart/'e'nd of a geometry piece to the file.
static void fprint_vtk_piece_geom
	(FILE* file,                                ///< The file.
	 const char sp_type,                        ///< Type indicator for 's'erial or 'p'arallel.
	 const char se_type,                        ///< Type indicator for 's'tart or 'e'nd.
	 const char geom_type,                      ///< Defined for \ref output_visualization_vtk_geom.
	 const struct const_Multiarray_d* xyz,      ///< Coordinates of the points (required for 's'erial, 's'tart).
	 const struct const_Plotting_Nodes* p_nodes ///< \ref Plotting_Nodes (required for 's'erial, 's'tart).
	);

/// \brief Print the 's'tart/'e'nd of a normals piece to the file.
static void fprint_vtk_piece_normals
	(FILE* file,                              ///< The file.
	 const char sp_type,                      ///< Type indicator for 's'erial or 'p'arallel.
	 const char se_type,                      ///< Type indicator for 's'tart or 'e'nd.
	 const struct const_Multiarray_d* xyz,    ///< Coordinates of the points (required for 's'erial, 's'tart).
	 const struct const_Multiarray_d* normals ///< Normal vector coordinates at each of the coordinates.
	);

/// \brief Print the 's'tart/'e'nd of a solution piece to the file.
static void fprint_vtk_piece_sol
	(FILE* file,                                 ///< The file.
	 const char sp_type,                         ///< Type indicator for 's'erial or 'p'arallel.
	 const char se_type,                         ///< Type indicator for 's'tart or 'e'nd.
	 const struct const_Multiarray_d* xyz,       ///< Coordinates of the points (required for 's'erial, 's'tart).
	 const struct const_Multiarray_d* sol,       ///< The solution variables    (required for 's'erial, 's'tart).
	 const struct const_Plotting_Nodes* p_nodes, ///< \ref Plotting_Nodes (required for 's'erial, 's'tart).
	 const int pde_index                         ///< \ref Test_Case::pde_index.
	);

static void output_visualization_vtk_geom
	(const char geom_type, const struct Simulation* sim, const bool use_test_case_path)
{
	assert(sim->elements->name == IL_ELEMENT_PLOTTING);

	static char output_part[STRLEN_MAX] = { 0, };
	if (!use_test_case_path) {
		sprintf(output_part,"%s%c","geom_",geom_type);
	} else {
		sprintf(output_part,"%s%c%s%c%s%s",
		        sim->pde_name,'/',sim->pde_spec,'/',"geom_e__",extract_name(sim->ctrl_name_full,true));
	}

	const char*const output_name = set_output_name(VIS_SOFTWARE_PARAVIEW,output_part);

	static char* extension_part = "vtu";
	if (sim->mpi_rank == 0) {
		FILE* p_file = fopen_sp_output_file('p',output_name,extension_part,sim->mpi_rank); // closed

		fprint_vtk_header_footer(p_file,true,'h',"UnstructuredGrid");

		fprint_vtk_piece_geom(p_file,'p','s',geom_type,NULL,NULL);

		for (int i = 0; i < sim->mpi_size; ++i)
			fprintf(p_file,"<Piece Source=\"%s_%d.vtu\"/>\n",extract_name(output_part,false),i);
		fprintf(p_file,"\n");

		fprint_vtk_header_footer(p_file,true,'f',"UnstructuredGrid");

		fclose(p_file);
	}

	FILE* s_file = fopen_sp_output_file('s',output_name,extension_part,sim->mpi_rank); // closed

	fprint_vtk_header_footer(s_file,false,'h',"UnstructuredGrid");
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Volume* base_volume   = (struct Volume*)curr;
		struct Solver_Volume* volume = (struct Solver_Volume*)curr;

		const struct const_Plotting_Element* element = (struct const_Plotting_Element*)base_volume->element;

		const int p = volume->p_ref;
		const struct Operator* cv0_vg_vp =
			(!base_volume->curved ? get_Multiarray_Operator(element->cv0_vgs_vp,(ptrdiff_t[]){0,0,p,1})
			                      : get_Multiarray_Operator(element->cv0_vgc_vp,(ptrdiff_t[]){0,0,p,p}) );

		const struct const_Multiarray_d* g_coef = volume->geom_coef;
		const struct const_Multiarray_d* xyz_p =
			constructor_mm_NN1_Operator_const_Multiarray_d(cv0_vg_vp,g_coef,'R','d',g_coef->order,NULL); // destructed

		fprint_vtk_piece_geom(s_file,'s','s',geom_type,xyz_p,element->p_nodes[p]);
		fprint_vtk_piece_geom(s_file,'s','e',geom_type,NULL,NULL);

		destructor_const_Multiarray_d(xyz_p);
	}
	fprint_vtk_header_footer(s_file,false,'f',"UnstructuredGrid");

	fclose(s_file);
}

static void output_visualization_vtk_normals (const struct Simulation* sim)
{
	static const char* output_part = "normals";
	const char*const output_name = set_output_name(VIS_SOFTWARE_PARAVIEW,output_part);

	static char* extension_part = "vtp";
	if (sim->mpi_rank == 0) {
		FILE* p_file = fopen_sp_output_file('p',output_name,extension_part,sim->mpi_rank); // closed

		fprint_vtk_header_footer(p_file,true,'h',"PolyData");
		fprint_vtk_piece_normals(p_file,'p','s',NULL,NULL);

		for (int i = 0; i < sim->mpi_size; ++i)
			fprintf(p_file,"<Piece Source=\"%s_%d.vtp\"/>\n",output_part,i);
		fprintf(p_file,"\n");

		fprint_vtk_header_footer(p_file,true,'f',"PolyData");

		fclose(p_file);
	}

	FILE* s_file = fopen_sp_output_file('s',output_name,extension_part,sim->mpi_rank); // closed

	fprint_vtk_header_footer(s_file,false,'h',"PolyData");
	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
		struct Solver_Face* face = (struct Solver_Face*)curr;

		fprint_vtk_piece_normals(s_file,'s','s',face->xyz_fc,face->normals_fc);
		fprint_vtk_piece_normals(s_file,'s','e',NULL,NULL);
	}
	fprint_vtk_header_footer(s_file,false,'f',"PolyData");

	fclose(s_file);
}

static void output_visualization_vtk_sol (const struct Simulation* sim)
{
/// \todo Add output for trace unknowns also if present.
	assert(sim->elements->name == IL_ELEMENT_PLOTTING);

	static char output_part[STRLEN_MAX] = { 0, };
	sprintf(output_part,"%s%c%s%c%s%s",
	        sim->pde_name,'/',sim->pde_spec,'/',"sol_v__",extract_name(sim->ctrl_name_full,true));

	const char*const output_name = set_output_name(VIS_SOFTWARE_PARAVIEW,output_part);

	struct Test_Case* test_case = (struct Test_Case*)sim->test_case_rc->tc;

	static char* extension_part = "vtu";
	if (sim->mpi_rank == 0) {
		FILE* p_file = fopen_sp_output_file('p',output_name,extension_part,sim->mpi_rank); // closed

		fprint_vtk_header_footer(p_file,true,'h',"UnstructuredGrid");

		fprint_vtk_piece_sol(p_file,'p','s',NULL,NULL,NULL,test_case->pde_index);

		for (int i = 0; i < sim->mpi_size; ++i)
			fprintf(p_file,"<Piece Source=\"%s_%d.vtu\"/>\n",extract_name(output_part,false),i);
		fprintf(p_file,"\n");

		fprint_vtk_header_footer(p_file,true,'f',"UnstructuredGrid");

		fclose(p_file);
	}

	FILE* s_file = fopen_sp_output_file('s',output_name,extension_part,sim->mpi_rank); // closed

	fprint_vtk_header_footer(s_file,false,'h',"UnstructuredGrid");
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Volume* base_volume   = (struct Volume*)curr;
		struct Solver_Volume* volume = (struct Solver_Volume*)curr;

		const struct const_Plotting_Element* element = (struct const_Plotting_Element*)base_volume->element;

		const int p = volume->p_ref;
		const struct Operator* cv0_vg_vp =
			(!base_volume->curved ? get_Multiarray_Operator(element->cv0_vgs_vp,(ptrdiff_t[]){0,0,p,1})
			                      : get_Multiarray_Operator(element->cv0_vgc_vp,(ptrdiff_t[]){0,0,p,p}) );

		const struct const_Multiarray_d* g_coef = volume->geom_coef;
		const struct const_Multiarray_d* xyz_p =
			constructor_mm_NN1_Operator_const_Multiarray_d(cv0_vg_vp,g_coef,'R','d',g_coef->order,NULL); // destructed

		const struct Operator* cv0_vs_vp = get_Multiarray_Operator(element->cv0_vs_vp,(ptrdiff_t[]){0,0,p,p});

		const struct const_Multiarray_d* s_coef = (const struct const_Multiarray_d*)volume->sol_coef;
		const struct const_Multiarray_d* sol_p =
			constructor_mm_NN1_Operator_const_Multiarray_d(cv0_vs_vp,s_coef,'C','d',s_coef->order,NULL); // destructed

		fprint_vtk_piece_sol(s_file,'s','s',xyz_p,sol_p,element->p_nodes[p],test_case->pde_index);
		fprint_vtk_piece_sol(s_file,'s','e',NULL,NULL,NULL,test_case->pde_index);

		destructor_const_Multiarray_d(xyz_p);
		destructor_const_Multiarray_d(sol_p);
	}
	fprint_vtk_header_footer(s_file,false,'f',"UnstructuredGrid");

	fclose(s_file);
}

// Level 2 ********************************************************************************************************** //

/** \brief Print the array of cummulative sum of `ext_0` of the \ref Vector_T\*s of the input
 *         \ref const_Multiarray_Vector_T to a file with the input number of tabs before each row. */
void fprint_const_Multiarray_Vector_i_offsets
	(FILE* file,                               ///< The file.
	 const int n_tab,                          ///< The number of tabs.
	 const struct const_Multiarray_Vector_i* a ///< Standard.
	);

/// \brief Print a vtk Points entry to the file.
void fprint_vtk_Points
	(FILE* file,                          ///< The file.
	 const char sp_type,                  ///< Type indicator for 's'erial or 'p'arallel.
	 const struct const_Multiarray_d* xyz ///< Coordinates of the points (required for 's'erial).
	);

/// \brief Print a vtk DataArray entry for a `double` Vector/Scalar to the file.
void fprint_vtk_DataArray_d
	(FILE* file,                            ///< The file.
	 const char sp_type,                    ///< Type indicator for 's'erial or 'p'arallel.
	 const char*const data_name,            ///< The name of the data.
	 const struct const_Multiarray_d* data, ///< Coordinates of the points (required for 's'erial).
	 const bool include_hf,                 ///< Flag for whether the hearder/footer should be included.
	 const char data_type                   ///< The type of data. Options: 'v'ector, 's'calar.
	);

/// \brief Print the solution data of an advection piece to the file.
static void fprint_vtk_piece_sol_advection
	(FILE* file,                          ///< Defined for \ref fprint_vtk_piece_sol.
	 const char sp_type,                  ///< Defined for \ref fprint_vtk_piece_sol.
	 const struct const_Multiarray_d* sol ///< Defined for \ref fprint_vtk_piece_sol.
	);

/// \brief Print the solution data of an euler piece to the file.
static void fprint_vtk_piece_sol_euler
	(FILE* file,                          ///< Defined for \ref fprint_vtk_piece_sol.
	 const char sp_type,                  ///< Defined for \ref fprint_vtk_piece_sol.
	 const struct const_Multiarray_d* sol ///< Defined for \ref fprint_vtk_piece_sol.
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

	strcat(output_name,name_spec);

	return output_name;
}

static void fprintf_tn (FILE* file, const int n_tabs, const char*const string_i)
{
	for (int i = 0; i < n_tabs; i++)
		fprintf(file,"\t");
	fprintf(file,"%s\n",string_i);
}

static void fprint_vtk_header_footer
	(FILE* file, const bool is_parallel, const char hf_type, const char*const vtk_type_part)
{
	static char vtk_type[STRLEN_MIN] = { 0, },
	            string_i[STRLEN_MAX] = { 0, };

	if (is_parallel)
		sprintf(vtk_type,"%c%s",'P',vtk_type_part);
	else
		strcpy(vtk_type,vtk_type_part);

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

static void fprint_vtk_piece_geom
	(FILE* file, const char sp_type, const char se_type, const char geom_type, const struct const_Multiarray_d* xyz,
	const struct const_Plotting_Nodes* p_nodes)
{
	assert(sp_type == 's' || sp_type == 'p');
	assert(se_type == 's' || se_type == 'e');
	assert(geom_type == 'v' || geom_type == 'e');

	if (sp_type == 'p') {
		if (se_type == 's') {
			fprint_vtk_Points(file,sp_type,NULL);
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
			fprint_vtk_Points(file,sp_type,xyz);

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

static void fprint_vtk_piece_normals
	(FILE* file, const char sp_type, const char se_type, const struct const_Multiarray_d* xyz,
	 const struct const_Multiarray_d* normals)
{
	assert(sp_type == 's' || sp_type == 'p');
	assert(se_type == 's' || se_type == 'e');

	if (sp_type == 'p') {
		if (se_type == 's') {
			fprint_vtk_Points(file,sp_type,NULL);
			fprint_vtk_DataArray_d(file,sp_type,"Normals",normals,true,'v');
			fprintf(file,"\n");
		} else if (se_type == 'e') {
			EXIT_UNSUPPORTED;
		}
	} else if (sp_type == 's') {
		if (se_type == 's') {
			fprintf(file,"\n<Piece NumberOfPoints=\"%td\" NumberOfVerts=\"0\" NumberOfLines=\"0\" "
			             "NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n",xyz->extents[0]);
			fprint_vtk_Points(file,sp_type,xyz);
			fprint_vtk_DataArray_d(file,sp_type,"Normals",normals,true,'v');
		} else if (se_type == 'e') {
			fprintf_tn(file,0,"</Piece>\n");
		}
	}
}

static void fprint_vtk_piece_sol
	(FILE* file, const char sp_type, const char se_type, const struct const_Multiarray_d* xyz,
	 const struct const_Multiarray_d* sol, const struct const_Plotting_Nodes* p_nodes, const int pde_index)
{
	assert(sp_type == 's' || sp_type == 'p');
	assert(se_type == 's' || se_type == 'e');

	if (sp_type == 'p') {
		if (se_type == 's') {
			fprint_vtk_Points(file,sp_type,NULL);
			fprintf_tn(file,1,"<PCells>");
				fprintf_tn(file,2,"<PDataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\"/>");
				fprintf_tn(file,2,"<PDataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\"/>");
				fprintf_tn(file,2,"<PDataArray type=\"UInt8\" Name=\"types\" format=\"ascii\"/>");
			fprintf_tn(file,1,"</PCells>");

			fprintf_tn(file,1,"<PPointData Scalars=\"Scalars\" Vectors=\"Vectors\">");
			switch (pde_index) {
				case PDE_ADVECTION: fprint_vtk_piece_sol_advection(file,sp_type,sol); break;
				case PDE_EULER:     fprint_vtk_piece_sol_euler(file,sp_type,sol);     break;
				default: EXIT_ERROR("Unsupported: %d\n",pde_index); break;
			}
			fprintf_tn(file,1,"</PPointData>");

			fprintf(file,"\n");
		} else if (se_type == 'e') {
			EXIT_UNSUPPORTED;
		}
	} else if (sp_type == 's') {
		if (se_type == 's') {
			const struct const_Multiarray_Vector_i* connect = p_nodes->connect;
			const struct const_Vector_i* vtk_types          = p_nodes->vtk_types;

			fprintf(file,"\n<Piece NumberOfPoints=\"%td\" NumberOfCells=\"%td\">\n",
			             xyz->extents[0],connect->extents[0]);
			fprint_vtk_Points(file,sp_type,xyz);

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

			fprintf_tn(file,1,"<PointData Scalars=\"Scalars\" Vectors=\"Vectors\">");
			switch (pde_index) {
				case PDE_ADVECTION: fprint_vtk_piece_sol_advection(file,sp_type,sol); break;
				case PDE_EULER:      fprint_vtk_piece_sol_euler(file,sp_type,sol);    break;
				default: EXIT_ERROR("Unsupported: %d\n",pde_index); break;
			}
			fprintf_tn(file,1,"</PointData>");
		} else if (se_type == 'e') {
			fprintf_tn(file,0,"</Piece>\n");
		}
	}
}

// Level 3 ********************************************************************************************************** //

/** \brief Print a \ref const_Multiarray_T to a file with the input number of tabs before each row, and padding rows
 *         with zeroes until they have 3 entries. */
void fprint_const_Multiarray_d_vtk_point
	(FILE* file,                         ///< The file.
	 const int n_tab,                    ///< The number of tabs.
	 const struct const_Multiarray_d* a, ///< Standard.
	 const char data_type                ///< Defined for \ref fprint_vtk_DataArray_d.
	);

void fprint_vtk_Points (FILE* file, const char sp_type, const struct const_Multiarray_d* xyz)
{
	// Note: Points **must** have 3 values.
	assert(sp_type == 's' || sp_type == 'p');

	static char points_name[STRLEN_MIN] = { 0, };
	if (sp_type == 's')
		strcpy(points_name,"Points");
	else
		sprintf(points_name,"%c%s",'P',"Points");

	fprintf(file,"\t<%s>\n",points_name);
	fprintf_tn(file,2,"<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">");
	if (sp_type == 's')
		fprint_const_Multiarray_d_vtk_point(file,2,xyz,'v');
	fprintf_tn(file,2,"</DataArray>");
	fprintf(file,"\t</%s>\n",points_name);
}

void fprint_vtk_DataArray_d
	(FILE* file, const char sp_type, const char*const data_name, const struct const_Multiarray_d* data,
	 const bool include_hf, const char data_type)
{
	assert(sp_type == 's' || sp_type == 'p');
	assert(data_type == 'v' || data_type == 's');

	static char pointdata_name[STRLEN_MIN] = { 0, };
	if (include_hf) {
		if (sp_type == 's')
			strcpy(pointdata_name,"PointData");
		else
			sprintf(pointdata_name,"%c%s",'P',"PointData");
		fprintf(file,"\t<%s Vectors=\"%s\">\n",pointdata_name,data_name);
	}

	if (data_type == 'v') {
		// Note: Vector Points **must** have 3 values.
		fprintf(file,"\t\t<DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"3\" format=\"ascii\">\n",
		        data_name);
	} else if (data_type == 's') {
		fprintf(file,"\t\t<DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\">\n",data_name);
	}
	if (sp_type == 's')
		fprint_const_Multiarray_d_vtk_point(file,2,data,data_type);
	fprintf_tn(file,2,"</DataArray>");

	if (include_hf)
		fprintf(file,"\t</%s>\n",pointdata_name);
}

void fprint_const_Multiarray_Vector_i_offsets (FILE* file, const int n_tab, const struct const_Multiarray_Vector_i* a)
{
	const ptrdiff_t size = compute_size(a->order,a->extents);

	struct Vector_i* a_V = constructor_empty_Vector_i(size); // destructed
	ptrdiff_t sum = 0;
	for (ptrdiff_t i = 0; i < size; ++i) {
		a_V->data[i] = (int)(sum + a->data[i]->ext_0);
		sum = a_V->data[i];
	}
	fprint_Vector_i(file,n_tab,a_V);
	destructor_Vector_i(a_V);
}

static void fprint_vtk_piece_sol_euler (FILE* file, const char sp_type, const struct const_Multiarray_d* sol)
{
	if (sp_type == 'p') {
		fprint_vtk_DataArray_d(file,sp_type,"rho",NULL,false,'s');
		fprint_vtk_DataArray_d(file,sp_type,"velocity",NULL,false,'v');
		fprint_vtk_DataArray_d(file,sp_type,"p",NULL,false,'s');
		fprint_vtk_DataArray_d(file,sp_type,"s",NULL,false,'s');
		fprint_vtk_DataArray_d(file,sp_type,"mach",NULL,false,'s');
	} else if (sp_type == 's') {
		convert_variables((struct Multiarray_d*)sol,'c','p');

		const ptrdiff_t ext_0 = sol->extents[0],
		                n_var = sol->extents[1];

		struct Multiarray_d* var = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){ext_0,1}); // destructed.
		double* data = var->data;

		// Primitive variables
		var->extents[1] = 1;
		var->data = (double*) get_col_const_Multiarray_d(0,sol);
		fprint_vtk_DataArray_d(file,sp_type,"rho",(struct const_Multiarray_d*)var,false,'s');

		var->extents[1] = n_var-2;
		var->data = (double*) get_col_const_Multiarray_d(1,sol);
		fprint_vtk_DataArray_d(file,sp_type,"velocity",(struct const_Multiarray_d*)var,false,'v');

		var->extents[1] = 1;
		var->data = (double*) get_col_const_Multiarray_d(n_var-1,sol);
		fprint_vtk_DataArray_d(file,sp_type,"p",(struct const_Multiarray_d*)var,false,'s');

		// Additional variables
		var->data = data;

		var->extents[1] = 1;
		compute_entropy(var,sol,'p');
		fprint_vtk_DataArray_d(file,sp_type,"s",(struct const_Multiarray_d*)var,false,'s');

		compute_mach(var,sol,'p');
		fprint_vtk_DataArray_d(file,sp_type,"mach",(struct const_Multiarray_d*)var,false,'s');

		destructor_Multiarray_d(var);
	} else {
		EXIT_ERROR("Unsupported: %c\n",sp_type);
	}
}

static void fprint_vtk_piece_sol_advection (FILE* file, const char sp_type, const struct const_Multiarray_d* sol)
{
	if (sp_type == 'p') {
		fprint_vtk_DataArray_d(file,sp_type,"u",NULL,false,'s');
	} else if (sp_type == 's') {
		assert(sol->extents[1] == 1);

		fprint_vtk_DataArray_d(file,sp_type,"u",sol,false,'s');
	} else {
		EXIT_ERROR("Unsupported: %c\n",sp_type);
	}
}

// Level 4 ********************************************************************************************************** //

void fprint_const_Multiarray_d_vtk_point
	(FILE* file, const int n_tab, const struct const_Multiarray_d* a, const char data_type)
{
	assert(data_type == 'v' || data_type == 's');

	const int order               = a->order;
	const ptrdiff_t*const extents = a->extents;

	assert(order == 2);

	const ptrdiff_t ext_0 = extents[0],
	                ext_1 = extents[1];

	if (data_type == 'v') {
		const bool transpose_Ma = ( a->layout == 'R' ? false : true );
		if (transpose_Ma)
			transpose_Multiarray_d((struct Multiarray_d*)a,true);

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
	} else if (data_type == 's') {
		assert(ext_1 == 1);

		const double* data = a->data;

		bool new_line = true;
		for (ptrdiff_t i = 0; i < ext_0; ++i) {
			if (new_line) {
				for (int j = 0; j < n_tab; ++j)
					fprintf(file,"\t");
				new_line = false;
			}

			fprintf(file," % .8e",data[i]);

			if ((i+1)%8 == 0 || (i == ext_0-1)) {
				fprintf(file,"\n");
				new_line = true;
			}
		}
	} else {
		EXIT_ERROR("Unsupported: %c.\n",data_type);
	}
}
