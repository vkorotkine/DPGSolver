// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/** \file
 */

#include "test_unit_bases.h"

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

#include "test_base.h"

// Static function declarations ************************************************************************************* //

/// \brief Provides unit tests for the tensor-product basis functions.
static void test_unit_basis_tensor_product
	(struct Test_Info*const test_info ///< \ref Test_Info.
	);

// Interface functions ********************************************************************************************** //

void test_unit_bases (struct Test_Info*const test_info)
{
	test_unit_basis_tensor_product(test_info);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void test_unit_basis_tensor_product (struct Test_Info*const test_info)
{
	char* test_name = "Bases - TP";
	bool pass = true;

	pass = false;

	test_increment_and_print(test_info,pass,test_name);
}
