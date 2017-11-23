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
#include "definitions_alloc.h"
#include "definitions_core.h"
#include "definitions_intrusive.h"
#include "simulation/solvers/adaptation/definitions_adaptation.h"

#include "test_base.h"

#include "adaptation.h"
#include "computational_elements.h"
#include "compute_error.h"
#include "const_cast.h"
#include "file_processing.h"
#include "simulation.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

struct Integration_Test_Info* constructor_Integration_Test_Info (const char*const ctrl_name)
{
	struct Integration_Test_Info* int_test_info = malloc(sizeof *int_test_info); // returned
	int_test_info->ctrl_name = ctrl_name;

	// Read information
	const char* ctrl_name_full = set_ctrl_name_full(ctrl_name);
	FILE *ctrl_file = fopen_checked(ctrl_name_full);

	int ml_ctrl[2];
	int p_ref_ctrl[2];

	char line[STRLEN_MAX];
	while (fgets(line,sizeof(line),ctrl_file)) {
		if (strstr(line,"mesh_level")) read_skip_const_i_1(line,1,ml_ctrl,2);
		if (strstr(line,"p_ref"))      read_skip_const_i_1(line,1,p_ref_ctrl,2);

		if (strstr(line,"ml_range_test")) read_skip_const_i_1(line,1,int_test_info->ml,2);
		if (strstr(line,"p_range_test"))  read_skip_const_i_1(line,1,int_test_info->p_ref,2);
	}
	fclose(ctrl_file);

	const_cast_i(&int_test_info->adapt_type,compute_adapt_type(p_ref_ctrl,ml_ctrl));

	return int_test_info;
}

void destructor_Integration_Test_Info (struct Integration_Test_Info* int_test_info)
{
	free(int_test_info);
}

void structor_simulation
	(struct Simulation** sim, const char mode, const int adapt_type, const int p, const int ml, const int p_prev,
	 const int ml_prev, const char*const ctrl_name)
{
	assert(mode == 'c' || mode == 'd');

	switch (adapt_type) {
	case ADAPT_0:
		if (mode == 'c') {
			*sim = constructor_Simulation(ctrl_name); // destructed
			constructor_derived_computational_elements(*sim,IL_SOLVER); // destructed
		} else if (mode == 'd') {
			destructor_derived_computational_elements(*sim,IL_BASE);
			destructor_Simulation(*sim);
		}
		break;
	case ADAPT_P:
		if (ml != ml_prev) {
			structor_simulation(sim,mode,ADAPT_0,p,ml,p_prev,ml_prev,ctrl_name);
		} else {
			if (mode == 'c')
				adapt_hp(*sim,ADAPT_S_P_REFINE);
		}
		break;
	case ADAPT_H:
		EXIT_ADD_SUPPORT;
		if (p != p_prev) {
			structor_simulation(sim,mode,ADAPT_0,p,ml,p_prev,ml_prev,ctrl_name);
		} else {
			; // h-adapt
		}
		break;
	case ADAPT_HP:
		EXIT_ADD_SUPPORT;
		if (ml != ml_prev)
			{ ; } // h-adapt
		if (p != p_prev)
			{ ; } // p-adapt
		break;
	default:
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
		index[4] = '0'+ml;
		assert(!isdigit(index[5]));
		index = strstr(file_name_curr,"__p");
		index[3] = '0'+p;
		assert(!isdigit(index[4]));
		break;
	case ADAPT_P:
		assert(strstr(file_name_curr,"__p") == NULL);
		index = strstr(file_name_curr,"__ml");
		index[4] = '0'+ml;
		assert(!isdigit(index[5]));
		break;
	case ADAPT_H:
		assert(strstr(file_name_curr,"__ml") == NULL);
		index = strstr(file_name_curr,"__p");
		index[3] = '0'+p;
		assert(!isdigit(index[4]));
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

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
