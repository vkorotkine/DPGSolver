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

#include "test_integration.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "macros.h"
#include "definitions_adaptation.h"
#include "definitions_alloc.h"
#include "definitions_core.h"
#include "definitions_intrusive.h"
#include "definitions_simulation.h"

#include "test_base.h"
#include "test_complex_computational_elements.h"
#include "test_complex_geometry.h"
#include "test_complex_solution.h"
#include "test_complex_test_case.h"
#include "test_support_multiarray.h"
#include "test_support_vector.h"

#include "multiarray.h"
#include "vector.h"

#include "adaptation.h"
#include "computational_elements.h"
#include "compute_error.h"
#include "const_cast.h"
#include "file_processing.h"
#include "geometry.h"
#include "simulation.h"
#include "solution.h"

// Static function declarations ************************************************************************************* //

/** Flag used to specify that the simulation should be set up as if the adaptation type was \ref ADAPT_0 but with
 *  potential h-refinement for the initial mesh. */
#define ADAPT_0_FOR_H 10

// Interface functions ********************************************************************************************** //

struct Integration_Test_Info* constructor_Integration_Test_Info (const char*const ctrl_name)
{
	struct Integration_Test_Info* int_test_info = calloc(1,sizeof *int_test_info); // returned
	int_test_info->ctrl_name = ctrl_name;

	// Read information
	char line[STRLEN_MAX];

	const char* ctrl_name_full = set_ctrl_name_full(ctrl_name);
	FILE* ctrl_file = fopen_checked(ctrl_name_full);

	int ml_ctrl[2];
	int p_ref_ctrl[2];

	while (fgets(line,sizeof(line),ctrl_file)) {
		if (strstr(line,"mesh_level")) read_skip_const_i_1(line,1,ml_ctrl,2);
		if (strstr(line,"p_ref"))      read_skip_const_i_1(line,1,p_ref_ctrl,2);

		if (strstr(line,"ml_range_test")) read_skip_const_i_1(line,1,int_test_info->ml,2);
		if (strstr(line,"p_range_test"))  read_skip_const_i_1(line,1,int_test_info->p_ref,2);
	}
	fclose(ctrl_file);

	if (p_ref_ctrl[0] != int_test_info->p_ref[0])
		EXIT_ERROR("Mismatch in p: %d != %d.\n",p_ref_ctrl[0],int_test_info->p_ref[0]);
	if (ml_ctrl[0] != int_test_info->ml[0])
		EXIT_ERROR("Mismatch in ml: %d != %d.\n",ml_ctrl[0],int_test_info->ml[0]);

	const_cast_i(&int_test_info->adapt_type,compute_adapt_type(p_ref_ctrl,ml_ctrl));

	return int_test_info;
}

void destructor_Integration_Test_Info (struct Integration_Test_Info* int_test_info)
{
	free(int_test_info);
}

void structor_simulation
	(struct Simulation** sim, const char mode, const int adapt_type, const int p, const int ml, const int p_prev,
	 const int ml_prev, const char*const ctrl_name, const char type_rc)
{
	assert(mode == 'c' || mode == 'd');
	assert(type_rc == 'r' || type_rc == 'c');

	switch (adapt_type) {
	case ADAPT_0: { // fallthrough
	} case ADAPT_0_FOR_H: {
		if (mode == 'c') {
			if (*sim != NULL)
				structor_simulation(sim,'d',ADAPT_0,p,ml,p_prev,ml_prev,ctrl_name,type_rc);
			*sim = constructor_Simulation(ctrl_name); // destructed
			constructor_derived_Elements(*sim,IL_ELEMENT_SOLVER); // destructed
			switch (type_rc) {
			case 'r':
				constructor_derived_computational_elements(*sim,IL_SOLVER); // destructed
				set_up_solver_geometry(*sim);
				set_initial_solution(*sim);
				break; // dest.
			case 'c':
				constructor_derived_computational_elements_c(*sim,IL_SOLVER); // destructed
				convert_to_Test_Case_rc(*sim,'c'); // converted back
				set_up_solver_geometry_c(*sim);
				set_initial_solution_c(*sim);
				break;
			default:
				EXIT_ERROR("Unsupported: %c\n",type_rc);
				break;
			}
			if (adapt_type == ADAPT_0_FOR_H)
				adapt_initial_mesh_if_required(*sim);
		} else if (mode == 'd') {
			switch (type_rc) {
			case 'r':
				destructor_derived_computational_elements(*sim,IL_BASE); // destructed
				break;
			case 'c':
				convert_to_Test_Case_rc(*sim,'r');
				destructor_derived_computational_elements_c(*sim,IL_BASE); // destructed
				break;
			default:
				EXIT_ERROR("Unsupported: %c\n",type_rc);
				break;
			}
			destructor_derived_Elements(*sim,IL_ELEMENT);
			destructor_Simulation(*sim);
		}
		break;
	} case ADAPT_P:
		assert(mode == 'c');
		if (ml != ml_prev)
			structor_simulation(sim,mode,ADAPT_0,p,ml,p_prev,ml_prev,ctrl_name,type_rc);
		else
			adapt_hp(*sim,ADAPT_S_P_REFINE,NULL);
		break;
	case ADAPT_H:
		assert(mode == 'c');
		static bool entered = false;
		if (!entered) {
			structor_simulation(sim,mode,ADAPT_0_FOR_H,p,ml,p_prev,ml_prev,ctrl_name,type_rc);
			entered = true;
			return;
		}

		if (ml > ml_prev) {
			assert(ml-ml_prev == 1);
			adapt_hp(*sim,ADAPT_S_H_REFINE,NULL);
		} else {
			EXIT_ERROR("Should only be performing h-refinement.\n");
		}
		break;
	case ADAPT_HP: {
		assert(mode == 'c');
		static bool entered = false;
		if (!entered) {
			structor_simulation(sim,mode,ADAPT_0_FOR_H,p,ml,p_prev,ml_prev,ctrl_name,type_rc);
			entered = true;
			return;
		}

		if (ml > ml_prev) {
			assert(ml-ml_prev == 1);
			adapt_hp(*sim,ADAPT_S_H_REFINE,NULL);
		} else if (ml < ml_prev) {
			for (int ml_curr = ml; ml_curr != ml; --ml_curr)
				adapt_hp(*sim,ADAPT_S_H_COARSE,NULL);
		}
		if (p > p_prev) {
			assert(p-p_prev == 1);
			adapt_hp(*sim,ADAPT_S_P_REFINE,NULL);
		} else if (p < p_prev) {
			for (int p_curr = p_prev; p_curr != p; --p_curr)
				adapt_hp(*sim,ADAPT_S_P_COARSE,NULL);
		}
		break;
	} default:
		EXIT_ERROR("Unsupported: %d\n",adapt_type);
		break;
	}

	if (mode == 'c')
		set_ml_p_curr(ml,p,*sim);
}

const char* set_file_name_curr
	(const int adapt_type, const int p, const int ml, const bool add_missing, const char*const file_name)
{
	static char file_name_curr[STRLEN_MAX] = { 0, };
	strcpy(file_name_curr,file_name);

	assert(p  >= 0 && p  <= 9); // Constraints are only imposed because of the currently limited string processing.
	assert(ml >= 0 && ml <= 9);

	char* index = NULL;
	switch (adapt_type) {
	case ADAPT_0:
		if (add_missing)
			correct_file_name_ml_p(ml,p,file_name_curr);
		index = strstr(file_name_curr,"__ml");
		assert(index);
		index[4] = (char)('0'+ml);
		assert(!isdigit(index[5]));
		index = strstr(file_name_curr,"__p");
		assert(index);
		index[3] = (char)('0'+p);
		assert(!isdigit(index[4]));
		break;
	case ADAPT_P:
		assert(strstr(file_name_curr,"__p") == NULL);
		index = strstr(file_name_curr,"__ml");
		if (index) {
			index[4] = (char)('0'+ml);
			assert(!isdigit(index[5]));
		}
		break;
	case ADAPT_H:
		assert(strstr(file_name_curr,"__ml") == NULL);
		index = strstr(file_name_curr,"__p");
		if (index) {
			index[3] = (char)('0'+p);
			assert(!isdigit(index[4]));
		}
		break;
	case ADAPT_HP:
		assert(strstr(file_name_curr,"__ml") == NULL);
		assert(strstr(file_name_curr,"__p") == NULL);
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",adapt_type);
		break;
	}

	return file_name_curr;
}

void adapt_initial_mesh_if_required (struct Simulation*const sim)
{
	const int* count_to_find = (int[]) {2,3};
	int count_found = 0;

	struct Adaptation_Data adapt_data =
		{ .adapt_h = { 0, 0, },
		  .xyz_ve_refine = NULL,
		  .xyz_ve_ml = NULL,
		  .xyz_ve_p  = NULL,
		};

	char line[STRLEN_MAX];
	FILE* input_file = fopen_input('t',NULL,NULL); // closed
	while (fgets(line,sizeof(line),input_file)) {
		read_skip_name_i("n_adapt_h_refine",line,&adapt_data.adapt_h[0]);
		if (strstr(line,"xyz_ve_refine")) {
			++count_found;
			adapt_data.xyz_ve_refine = constructor_file_const_Multiarray_d(input_file,true); // destructed
		}
		if (strstr(line,"xyz_ve_mesh_level")) {
			++count_found;
			adapt_data.xyz_ve_ml = constructor_file_const_Vector_i(input_file,true); // destructed
		}
		if (strstr(line,"xyz_ve_polynomial_order")) {
			++count_found;
			adapt_data.xyz_ve_p = constructor_file_const_Vector_i(input_file,true); // destructed
		}
	}
	fclose(input_file);

	if (adapt_data.xyz_ve_refine == NULL) {
		assert(adapt_data.xyz_ve_ml == NULL || adapt_data.xyz_ve_p == NULL);
		return;
	}

	if (count_found < count_to_find[0] || count_found > count_to_find[1])
		EXIT_ERROR("Did not find the required number of variables (Found: %d/[%d,%d]).\n",
		           count_found,count_to_find[0],count_to_find[1]);

	adapt_hp(sim,ADAPT_S_XYZ_VE,&adapt_data);
	for (int i = 0; i < adapt_data.adapt_h[0]; ++i)
		adapt_hp(sim,ADAPT_S_H_REFINE,NULL);
	for (int i = 0; i < adapt_data.adapt_h[1]; ++i)
		adapt_hp(sim,ADAPT_S_H_COARSE,NULL);

	destructor_const_Multiarray_d(adapt_data.xyz_ve_refine);
	destructor_conditional_const_Vector_i(adapt_data.xyz_ve_ml);
	destructor_conditional_const_Vector_i(adapt_data.xyz_ve_p);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
