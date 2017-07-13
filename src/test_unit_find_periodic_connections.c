// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "test_unit_find_periodic_connections.h"

#include <stdlib.h>
#include <stdio.h>

#include "Parameters.h"
#include "Test.h"

#include "test_support.h"
#include "array_norm.h"
#include "gmsh_reader.h"

/*
 *	Purpose:
 *		Test correctness of implementation of find_periodic_connections.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

void test_unit_find_periodic_connections(void)
{
	unsigned int pass;

	/*
	 *	find_periodic_connections:
	 *
	 *		Input:
	 *
	 *			PVe = [1 0
	 *			       2 8
	 *			       1 3
	 *			       8 3
	 *			       3 2
	 *			       8 1
	 *			       5 8]
	 *
	 *			pve = 7
	 *
	 *		Expected progression/output:
	 *
	 *			Pass 1:
	 *
	 *				PVeMatches:
	 *				1 0 0 0 0 0 0 0
	 *				0 3 8 0 0 0 0 0
	 *				3 8 0 0 0 0 0 0
	 *				1 2 8 0 0 0 0 0
	 *				8 0 0 0 0 0 0 0
	 *				1 2 3 5 0 0 0 0
	 *
	 *				Updated PVe (Sorted):
	 *				0 1
	 *				0 3
	 *				0 8
	 *				1 2
	 *				1 3
	 *				1 5
	 *				1 8
	 *				2 3
	 *				2 5
	 *				2 8
	 *				3 5
	 *				3 8
	 *				5 8
	 *
	 *			Pass 2:
	 *
	 *				PVeMatches:
	 *				1 3 8 0 0 0 0 0
	 *				0 2 3 5 8 0 0 0
	 *				1 3 5 8 0 0 0 0
	 *				0 1 2 5 8 0 0 0
	 *				1 2 3 8 0 0 0 0
	 *				0 1 2 3 5 0 0 0
	 *
	 *				Updated PVe (Sorted):
	 *				0 1
	 *				0 2
	 *				0 3
	 *				0 5
	 *				0 8
	 *				1 2
	 *				1 3
	 *				1 5
	 *				1 8
	 *				2 3
	 *				2 5
	 *				2 8
	 *				3 5
	 *				3 8
	 *				5 8
	 *
	 *			Pass 3:
	 *
	 *				PVeMatches:
	 *				1 2 3 5 8 0 0 0
	 *				0 2 3 5 8 0 0 0
	 *				0 1 3 5 8 0 0 0
	 *				0 1 2 5 8 0 0 0
	 *				0 1 2 3 8 0 0 0
	 *				0 1 2 3 5 0 0 0
	 *
	 *				Updated PVe (Sorted):
	 *				No change => Done!
	 */

	unsigned int *PVe, pve,
	              PVeFinal[30] = { 0, 1, 0, 2, 0, 3, 0, 5, 0, 8, 1, 2, 1, 3, 1,
	                               5, 1, 8, 2, 3, 2, 5, 2, 8, 3, 5, 3, 8, 5, 8 };

	PVe = malloc(6*6*2 * sizeof *PVe); // free

	PVe[0] = 1;  PVe[1] = 0;
	PVe[2] = 2;  PVe[3] = 8;
	PVe[4] = 1;  PVe[5] = 3;
	PVe[6] = 8;  PVe[7] = 3;
	PVe[8] = 3;  PVe[9] = 2;
	PVe[10] = 8; PVe[11] = 1;
	PVe[12] = 5; PVe[13] = 8;

	pve = 7;

	find_periodic_connections(PVe,&pve,9);

	pass = 0;
	if (array_norm_diff_ui(30,PVeFinal,PVe,"Inf") < EPS)
		pass = 1;

	test_print2(pass,"find_periodic_connections:");

	free(PVe);
}
