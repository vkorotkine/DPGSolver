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

#include "test_base.h"
#include "test_support.h"
#include "test_support_matrix.h"

#include "macros.h"
#include "definitions_elements.h"
#include "definitions_tol.h"

#include "matrix.h"

#include "inverse_mapping.h"

// Static function declarations ************************************************************************************* //

/// \brief Container for inverse mapping data.
struct Inv_Map_Data {
	// Read from file.
	const struct const_Matrix_d* ve_rand; ///< Random vertices.
	const struct const_Matrix_d* ve_std;  ///< Vertices of the standard reference element.
	const struct const_Matrix_d* ve_symm; ///< Vertices of the symmetric reference element.
	const struct const_Matrix_d* rst_std; ///< rst reference coordinates for the standard reference element.

	// Computed
	int e_type; ///< The element type.
	const struct const_Matrix_d* xyz_rand; ///< xyz coordinates corresponding to \ref Inv_Map_Data::ve_rand.
	const struct const_Matrix_d* xyz_std;  ///< xyz coordinates corresponding to \ref Inv_Map_Data::ve_std.
	const struct const_Matrix_d* xyz_symm; ///< xyz coordinates corresponding to \ref Inv_Map_Data::ve_symm.
};

/** \brief Constructor for a \ref Inv_Map_Data container.
 *  \return See brief. */
static const struct Inv_Map_Data* constructor_file_Inv_Map_Data
	(const char*const data_name ///< String containing the name of the data file (without path).
	);

/// \brief Destructor for a \ref Inv_Map_Data container.
static void destructor_Inv_Map_Data
	(const struct Inv_Map_Data*const imd ///< Standard.
	);

// Interface functions ********************************************************************************************** //

/** \test Performs unit testing for the inverse p1 mapping (from physical to reference coordinate) functionality
 *        (\ref test_unit_inverse_mapping.c).
 *  \return 0 on success.
 *
 *  The test checks that the physical coordinates obtained from the p1 mapping of input reference coordinates and
 *  physical vertices are recovered exactly by the inverse mapping algorithm.
 */
int main
	(int argc,   ///< Standard.
	 char** argv ///< Standard.
	)
{
	assert_condition_message(argc == 2,"Invalid number of input arguments");

	const struct Inv_Map_Data*const imd = constructor_file_Inv_Map_Data(argv[1]); // destructed

	const int e = imd->e_type;
	const struct const_Matrix_d*const rst_rand = constructor_inverse_mapping(e,imd->ve_rand,imd->xyz_rand); // dest.
	const struct const_Matrix_d*const rst_std  = constructor_inverse_mapping(e,imd->ve_std ,imd->xyz_std ); // dest.
	const struct const_Matrix_d*const rst_symm = constructor_inverse_mapping(e,imd->ve_symm,imd->xyz_symm); // dest.

	bool pass = true;
	const double tol = 1e3*EPS;
	const bool* differences = (bool[])
		{ diff_const_Matrix_d(imd->xyz_symm,rst_rand,tol),
		  diff_const_Matrix_d(imd->xyz_symm,rst_std, tol),
		  diff_const_Matrix_d(imd->xyz_symm,rst_symm,tol),
		};
	if (check_diff(3,differences,&pass)) {
		if (differences[0])
			print_diff_const_Matrix_d(imd->xyz_symm,rst_rand,tol);
		if (differences[1])
			print_diff_const_Matrix_d(imd->xyz_symm,rst_std, tol);
		if (differences[2])
			print_diff_const_Matrix_d(imd->xyz_symm,rst_symm,tol);
	}
	assert_condition(pass);

	destructor_Inv_Map_Data(imd);
	destructor_const_Matrix_d(rst_rand);
	destructor_const_Matrix_d(rst_std);
	destructor_const_Matrix_d(rst_symm);

	OUTPUT_SUCCESS;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Compute and return the element index based on the input name.
 *  \return See brief. */
static int compute_e_type_from_name
	(const char*const name ///< The input string.
	);

/** \brief Constructor for physical xyz coordinates from standard reference element reference coordinates.
 *  \return See brief. */
static const struct const_Matrix_d* constructor_xyz_from_rst_std
	(const int e_type,                         ///< The element type.
	 const struct const_Matrix_d*const ve_xyz, ///< The xyz coordinates of the vertices.
	 const struct const_Matrix_d*const rst_std ///< The standard reference element coordinates.
	);

static const struct Inv_Map_Data* constructor_file_Inv_Map_Data (const char*const data_name)
{
	char d_name[STRLEN_MAX];
	sprintf(d_name,"%s%s","inverse_mapping/",data_name);
	const char*const f_name = set_data_file_name_unit(d_name);

	struct Inv_Map_Data*const imd = calloc(1,sizeof(*imd)); // free

	imd->ve_rand = constructor_file_name_const_Matrix_d("ve_random",f_name); // destructed
	imd->ve_std  = constructor_file_name_const_Matrix_d("ve_std",   f_name); // destructed
	imd->ve_symm = constructor_file_name_const_Matrix_d("ve_symm",  f_name); // destructed
	imd->rst_std = constructor_file_name_const_Matrix_d("rst_std",  f_name); // destructed

	const int e_type = compute_e_type_from_name(data_name);

	imd->e_type   = e_type;
	imd->xyz_rand = constructor_xyz_from_rst_std(e_type,imd->ve_rand,imd->rst_std); // destructed
	imd->xyz_std  = constructor_xyz_from_rst_std(e_type,imd->ve_std ,imd->rst_std); // destructed
	imd->xyz_symm = constructor_xyz_from_rst_std(e_type,imd->ve_symm,imd->rst_std); // destructed

	return imd;
}

static void destructor_Inv_Map_Data (const struct Inv_Map_Data*const imd)
{
	destructor_const_Matrix_d(imd->ve_rand);
	destructor_const_Matrix_d(imd->ve_std);
	destructor_const_Matrix_d(imd->ve_symm);
	destructor_const_Matrix_d(imd->rst_std);

	destructor_const_Matrix_d(imd->xyz_rand);
	destructor_const_Matrix_d(imd->xyz_std);
	destructor_const_Matrix_d(imd->xyz_symm);
	free((void*)imd);
}

// Level 1 ********************************************************************************************************** //

static int compute_e_type_from_name (const char*const name)
{
	if (strstr(name,"line") || strstr(name,"LINE"))
		return LINE;
	else if (strstr(name,"tri") || strstr(name,"TRI"))
		return TRI;
	else if (strstr(name,"quad") || strstr(name,"QUAD"))
		return QUAD;
	else if (strstr(name,"tet") || strstr(name,"TET"))
		return TET;
	else if (strstr(name,"hex") || strstr(name,"HEX"))
		return HEX;
	else if (strstr(name,"wedge") || strstr(name,"WEDGE"))
		return WEDGE;
	else if (strstr(name,"pyr") || strstr(name,"PYR"))
		return PYR;
	EXIT_ERROR("Did not find the name: %s\n",name);
}

static const struct const_Matrix_d* constructor_xyz_from_rst_std
	(const int e_type, const struct const_Matrix_d*const ve_xyz, const struct const_Matrix_d*const rst_std)
{
	const struct const_Matrix_d*const phi_rst = constructor_basis_std_p1(e_type,rst_std); // destructed

	const struct const_Matrix_d*const xyz_std =
		constructor_mm_const_Matrix_d('N','N',1.0,phi_rst,ve_xyz,'R'); // returned
	destructor_const_Matrix_d(phi_rst);

	return xyz_std;
}
