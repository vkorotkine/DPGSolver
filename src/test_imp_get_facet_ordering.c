// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>

#include "test.h"
#include "parameters.h"
#include "functions.h"

/*
 *	Purpose:
 *		Test correctness of implementation of get_facet_ordering.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

void test_imp_get_facet_ordering(void)
{
	unsigned int pass;

	/*
	 *	Input:
	 *
	 *		Nn, d, IndOrd, FType
	 *
	 *	Expected Output (d = 1, Trivial):
	 *
	 *		nOrd = [ 0 ]
	 *
	 */

	unsigned int i, Nn, d, FType, IndOrd,
	             *nOrd;

	// d == 1 (FType = POINT, All P)
	d     = 1;
	FType = POINT;

	// All P
	Nn    = 1;
	unsigned int nOrd10[1] = { 0 };

	nOrd = malloc(Nn * sizeof *nOrd); // free

	// case 0
	IndOrd = 0;

	get_facet_ordering(d,IndOrd,FType,Nn,0,NULL,nOrd);

	pass = 0;
	if (array_norm_diff_ui(Nn,nOrd,nOrd10,"Inf") < EPS)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("get_facet_ordering (d = 1, case 0):              ");
	test_print(pass);

	free(nOrd);

	/*
	 *	Expected Output (d = 2):
	 *
	 *		P2:
	 *			nOrd = [ see below ]
	 *
	 *			Possible node positions (with indices as in case 0):
	 *
	 *			0 - 2 - 1   1 - 2 - 0
	 *
	 *			case 0      case 1
	 *
	 *		P3:
	 *			nOrd = [ see below ]
	 *
	 *			Possible node positions (with indices as in case 0):
	 *
	 *			0 - 2 - 3 - 1   1 - 3 - 2 - 0
	 *
	 *			case 0          case 1
	 */

	// d == 2 (FType = LINE)
	d     = 2;
	FType = LINE;

	// P == 2
	Nn    = 3;
	unsigned int nOrd2P2[2][3] = { { 0, 1, 2}, { 1, 0, 2}};

	nOrd = malloc(Nn * sizeof *nOrd); // free

	// case 0-1
	for (IndOrd = 0; IndOrd <= 1; IndOrd++) {
		for (i = 0; i < Nn; i++)
			nOrd[i] = 0.0;

		get_facet_ordering(d,IndOrd,FType,Nn,0,NULL,nOrd);

		pass = 0;
		if (array_norm_diff_ui(Nn,nOrd,nOrd2P2[IndOrd],"Inf") < EPS)
			pass = 1, TestDB.Npass++;

		//     0         10        20        30        40        50
		printf("                   (d = 2, case %d (P2)):         ",IndOrd);
		test_print(pass);
	}
	free(nOrd);

	// P == 3
	Nn    = 4;
	unsigned int nOrd2P3[2][4] = { { 0, 1, 2, 3}, { 1, 0, 3, 2}};

	nOrd = malloc(Nn * sizeof *nOrd); // free

	// case 0-1
	for (IndOrd = 0; IndOrd <= 1; IndOrd++) {
		for (i = 0; i < Nn; i++)
			nOrd[i] = 0.0;

		get_facet_ordering(d,IndOrd,FType,Nn,0,NULL,nOrd);

		pass = 0;
		if (array_norm_diff_ui(Nn,nOrd,nOrd2P3[IndOrd],"Inf") < EPS)
			pass = 1, TestDB.Npass++;

		//     0         10        20        30        40        50
		printf("                   (d = 2, case %d (P3)):         ",IndOrd);
		test_print(pass);
	}
	free(nOrd);


	/*
	 *	Expected Output (d = 3, QUAD):
	 *
	 *		P2:
	 *			nOrd = [ see below ]
	 *
	 *			Possible node positions (with indices as in case 0):
	 *
	 *			3 - 5 - 4   4 - 5 - 3   0 - 2 - 1   1 - 2 - 0   1 - 7 - 4   4 - 7 - 1   0 - 6 - 3   3 - 6 - 0
	 *			|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |
	 *			6 - 8 - 7   7 - 8 - 6   6 - 8 - 7   7 - 8 - 6   2 - 8 - 5   5 - 8 - 2   2 - 8 - 5   5 - 8 - 2
	 *			|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |
	 *			0 - 2 - 1   1 - 2 - 0   3 - 5 - 4   4 - 5 - 3   0 - 6 - 3   3 - 6 - 0   1 - 7 - 4   4 - 7 - 1
	 *
	 *			case 0      case 1      case 2      case 3      case 4      case 5      case 6      case 7
	 *
	 *		P3:
	 *			nOrd = [ see below ]
	 *
	 *			Possible node positions (with indices as in case 0):
	 *
	 *			004 - 006 - 007 - 005   005 - 007 - 006 - 004   000 - 002 - 003 - 001   001 - 003 - 002 - 000
	 *			 |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |
	 *			012 - 014 - 015 - 013   013 - 015 - 014 - 012   008 - 010 - 011 - 009   009 - 011 - 010 - 008
	 *			 |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |
	 *			008 - 010 - 011 - 009   009 - 011 - 010 - 008   012 - 014 - 015 - 013   013 - 015 - 014 - 012
	 *			 |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |
	 *			000 - 002 - 003 - 001   001 - 003 - 002 - 000   004 - 006 - 007 - 005   005 - 007 - 006 - 004
	 *
	 *			case 0                  case 1                  case 2                  case 3
	 *
	 *
	 *			001 - 009 - 013 - 005   005 - 013 - 009 - 001   000 - 008 - 012 - 004   004 - 012 - 008 - 000
	 *			 |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |
	 *			003 - 011 - 015 - 007   007 - 015 - 011 - 003   002 - 010 - 014 - 006   006 - 014 - 010 - 002
	 *			 |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |
	 *			002 - 010 - 014 - 006   006 - 014 - 010 - 002   003 - 011 - 015 - 007   007 - 015 - 011 - 003
	 *			 |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |
	 *			000 - 008 - 012 - 004   004 - 012 - 008 - 000   001 - 009 - 013 - 005   005 - 013 - 009 - 001
	 *
	 *			case 4                  case 5                  case 6                  case 7
	 */

	// d == 3 (FType = QUAD)
	d     = 3;
	FType = QUAD;

	// P == 2
	Nn    = 9;
	unsigned int nOrd3P2Q[8][9] = {{ 0, 1, 2, 3, 4, 5, 6, 7, 8},
	                               { 1, 0, 2, 4, 3, 5, 7, 6, 8},
	                               { 3, 4, 5, 0, 1, 2, 6, 7, 8},
	                               { 4, 3, 5, 1, 0, 2, 7, 6, 8},
	                               { 0, 3, 6, 1, 4, 7, 2, 5, 8},
	                               { 3, 0, 6, 4, 1, 7, 5, 2, 8},
	                               { 1, 4, 7, 0, 3, 6, 2, 5, 8},
	                               { 4, 1, 7, 3, 0, 6, 5, 2, 8}};

	nOrd = malloc(Nn * sizeof *nOrd); // free

	// case 0-7
	for (IndOrd = 0; IndOrd <= 7; IndOrd++) {
		for (i = 0; i < Nn; i++)
			nOrd[i] = 0.0;

		get_facet_ordering(d,IndOrd,FType,Nn,0,NULL,nOrd);

		pass = 0;
		if (array_norm_diff_ui(Nn,nOrd,nOrd3P2Q[IndOrd],"Inf") < EPS)
			pass = 1, TestDB.Npass++;

		//     0         10        20        30        40        50
		printf("                   (d = 3 (QUAD), case %d (P2)):  ",IndOrd);
		test_print(pass);
	}
	free(nOrd);

	// P == 3
	Nn    = 16;
	unsigned int nOrd3P3Q[8][16] = {{  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15},
	                                {  1,  0,  3,  2,  5,  4,  7,  6,  9,  8, 11, 10, 13, 12, 15, 14},
	                                {  4,  5,  6,  7,  0,  1,  2,  3, 12, 13, 14, 15,  8,  9, 10, 11},
	                                {  5,  4,  7,  6,  1,  0,  3,  2, 13, 12, 15, 14,  9,  8, 11, 10},
	                                {  0,  4,  8, 12,  1,  5,  9, 13,  2,  6, 10, 14,  3,  7, 11, 15},
	                                {  4,  0, 12,  8,  5,  1, 13,  9,  6,  2, 14, 10,  7,  3, 15, 11},
	                                {  1,  5,  9, 13,  0,  4,  8, 12,  3,  7, 11, 15,  2,  6, 10, 14},
	                                {  5,  1, 13,  9,  4,  0, 12,  8,  7,  3, 15, 11,  6,  2, 14, 10}};

	nOrd = malloc(Nn * sizeof *nOrd); // free

	// case 0-7
	for (IndOrd = 0; IndOrd <= 7; IndOrd++) {
		for (i = 0; i < Nn; i++)
			nOrd[i] = 0.0;

		get_facet_ordering(d,IndOrd,FType,Nn,0,NULL,nOrd);

		pass = 0;
		if (array_norm_diff_ui(Nn,nOrd,nOrd3P3Q[IndOrd],"Inf") < EPS)
			pass = 1, TestDB.Npass++;

		//     0         10        20        30        40        50
		printf("                   (d = 3 (QUAD), case %d (P3)):  ",IndOrd);
		test_print(pass);
	}
	free(nOrd);

	/*
	 *	Expected Output (d = 3, TRI):
	 *
	 *		P2:
	 *			nOrd = [ see below ]
	 *
	 *			Possible node positions (with indices as in case 0):
	 *			*** Not identical node ordering the WSH nodes, but same symmetry orbits and hence analogous. ***
	 *
	 *			    2           0           1           1           0           2
	 *			   / \         / \         / \         / \         / \         / \
	 *			  5 - 4       3 - 5       4 - 3       4 - 5       3 - 4       5 - 3
	 *			 / \ / \     / \ / \     / \ / \     / \ / \     / \ / \     / \ / \
	 *			0 - 3 - 1   1 - 4 - 2   2 - 5 - 0   0 - 3 - 2   2 - 5 - 1   1 - 4 - 0
	 *
	 *			case 0      case 1      case 2      case 3      case 4      case 5
	 *
	 *		P3:
	 *			nOrd = [ see below ]
	 *
	 *			Possible node positions (with indices as in case 0):
	 *			*** Not identical node ordering the WSH nodes, but same symmetry orbits and hence analogous. ***
	 *
	 *			      2               0               1               1               0               2
	 *			     / \             / \             / \             / \             / \             / \
	 *			    5 - 7           3 - 8           4 - 6           4 - 8           3 - 7           5 - 6
	 *			   / \ / \         / \ / \         / \ / \         / \ / \         / \ / \         / \ / \
	 *			  8 - 9 - 4       6 - 9 - 5       7 - 9 - 3       7 - 9 - 5       6 - 9 - 4       8 - 9 - 3
	 *			 / \ / \ / \     / \ / \ / \     / \ / \ / \     / \ / \ / \     / \ / \ / \     / \ / \ / \
	 *			0 - 3 - 6 - 1   1 - 4 - 7 - 2   2 - 5 - 8 - 0   0 - 3 - 6 - 2   2 - 5 - 8 - 1   1 - 4 - 7 - 0
	 *
	 *			case 0          case 1          case 2          case 3          case 4          case 5
	 */

	// d == 3 (FType = TRI)
	unsigned int P, NnTRI, Ns, *symms;
	double       *rst, *w;

	d     = 3;
	FType = TRI;

	// P == 2 (WSH)
	P     = 2;
	Nn    = 6;
	unsigned int nOrd3P2T[6][6] = {{ 0, 1, 2, 3, 4, 5},
	                               { 1, 2, 0, 4, 5, 3},
	                               { 2, 0, 1, 5, 3, 4},
	                               { 0, 2, 1, 3, 5, 4},
	                               { 2, 1, 0, 5, 4, 3},
	                               { 1, 0, 2, 4, 3, 5}};

	cubature_TRI(&rst,&w,&symms,&NnTRI,&Ns,0,P,d-1,"WSH");
	free(rst);
	if (NnTRI != Nn)
		printf("Error: Wrong Number of Nodes in test_imp_get_facet_ordering (TRI P2).\n"), exit(1);

	nOrd = malloc(Nn * sizeof *nOrd); // free

	// case 0-5
	for (IndOrd = 0; IndOrd <= 5; IndOrd++) {
		for (i = 0; i < Nn; i++)
			nOrd[i] = 0.0;

		get_facet_ordering(d,IndOrd,FType,Nn,Ns,symms,nOrd);

		pass = 0;
		if (array_norm_diff_ui(Nn,nOrd,nOrd3P2T[IndOrd],"Inf") < EPS)
			pass = 1, TestDB.Npass++;

		//     0         10        20        30        40        50
		printf("                   (d = 3 (TRI), case %d (P3)):   ",IndOrd);
		test_print(pass);
	}

	free(symms);
	free(nOrd);

	// P == 3 (WSH)
	P     = 3;
	Nn    = 10;
	unsigned int nOrd3P3T[6][10] = {{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9},
	                                { 1, 2, 0, 4, 5, 3, 7, 8, 6, 9},
	                                { 2, 0, 1, 5, 3, 4, 8, 6, 7, 9},
	                                { 0, 2, 1, 3, 5, 4, 6, 8, 7, 9},
	                                { 2, 1, 0, 5, 4, 3, 8, 7, 6, 9},
	                                { 1, 0, 2, 4, 3, 5, 7, 6, 8, 9}};

	cubature_TRI(&rst,&w,&symms,&NnTRI,&Ns,0,P,d-1,"WSH");
	free(rst);
	if (NnTRI != Nn)
		printf("Error: Wrong Number of Nodes in test_imp_get_facet_ordering (TRI P2).\n"), exit(1);

	nOrd = malloc(Nn * sizeof *nOrd); // free

	// case 0-5
	for (IndOrd = 0; IndOrd <= 5; IndOrd++) {
		for (i = 0; i < Nn; i++)
			nOrd[i] = 0.0;

		get_facet_ordering(d,IndOrd,FType,Nn,Ns,symms,nOrd);

		pass = 0;
		if (array_norm_diff_ui(Nn,nOrd,nOrd3P3T[IndOrd],"Inf") < EPS)
			pass = 1, TestDB.Npass++;

		//     0         10        20        30        40        50
		printf("                   (d = 3 (TRI), case %d (P3)):   ",IndOrd);
		test_print(pass);
	}
	free(symms);
	free(nOrd);
}
