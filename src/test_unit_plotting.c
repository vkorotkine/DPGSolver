// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "test_unit_plotting.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "mkl.h"

#include "Parameters.h"
#include "Test.h"

#include "test_support.h"
#include "array_norm.h"
#include "plotting_element_info.h"

#include "array_print.h"

/*
 *	Purpose:
 *		Test correctness of implementation of plotting_element_info.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

void test_unit_plotting(void)
{
	unsigned int pass;

	/*	Tensor-Product:
	 *
	 *		Input:
	 *
	 *			P, d
	 *
	 *		Expected Output (d = 1):
	 *
	 *			P = 3:
	 *				rst = [ See below ]
	 *				Nn  = 4
	 *
	 *				Nodes:
	 *				0 - 1 - 2 - 3
	 *
	 *				connect:
	 *				[  0  1  0  0  0  0  0  0;
	 *				   1  2  0  0  0  0  0  0;
	 *				   2  3  0  0  0  0  0  0]
	 *
	 *				types:
	 *				[  3  3  3]
	 */

	unsigned int i, j, k, iMax, jMax, kMax, row, dim, Indr, d, P, Nn, NE;
	unsigned int *connect, *types, *connectE;
	double *rst;

	// TP (d = 1, P = 3)
	d = 1;
	P = 3;

	unsigned int NE13_TP, Nn13_TP;

	Nn13_TP = P+1;
	NE13_TP = pow(P,d);

	unsigned int connect13_TP[24] = {  0,  1,  0,  0,  0,  0,  0,  0,
	                                   1,  2,  0,  0,  0,  0,  0,  0,
	                                   2,  3,  0,  0,  0,  0,  0,  0},
	             types13_TP[3];
	double rst13_TP[4] = { -1.0, -1.0/3.0, 1.0/3.0, 1.0};

	for (i = 0; i < NE13_TP; i++)
		types13_TP[i] = 3;

//	for (i = 0, iMax = P+1; i < iMax; i++)
//		rst13_TP[i] = -1.0 + (2.0*i)/P;
//	rst13_TP[0] = -2.0;

	plotting_element_info(&rst,&connect,&types,&connectE,&Nn,&NE,P,LINE); // free

	pass = 0;
	if (array_norm_diff_d(Nn13_TP*d,rst13_TP,rst,"Inf")          < EPS &&
	    array_norm_diff_ui(NE13_TP*8,connect13_TP,connect,"Inf") < EPS &&
	    array_norm_diff_ui(NE13_TP,types13_TP,types,"Inf")       < EPS &&
	    Nn == Nn13_TP && NE == NE13_TP)
			pass = 1;

	test_print2(pass,"plotting_element_info (d1, P3, TP   ):");

	free(rst);
	free(connect);
	free(types);
	free(connectE);


	/*
	 *		Expected Output (d = 2):
	 *
	 *			P = 3:
	 *				rst = [ See below ]
	 *				Nn  = 16
	 *
	 *				Nodes:
	 *				12 - 13 - 14 - 15
	 *				 |    |    |    |
	 *				 8 -  9 - 10 - 11
	 *				 |    |    |    |
	 *				 4 -  5 -  6 -  7
	 *				 |    |    |    |
	 *				 0 -  1 -  2 -  3
	 *
	 *				connect:
	 *				[  0  1  5  4  0  0  0  0;
	 *				   1  2  6  5  0  0  0  0;
	 *				   2  3  7  6  0  0  0  0;
	 *				   4  5  9  8  0  0  0  0;
	 *				   5  6 10  9  0  0  0  0;
	 *				   6  7 11 10  0  0  0  0;
	 *				   8  9 13 12  0  0  0  0;
	 *				   9 10 14 13  0  0  0  0;
	 *				  10 11 15 14  0  0  0  0]
	 *
	 *				types:
	 *				[  9  9  9  9  9  9  9  9  9]
	 */

	// TP (d = 2, P = 3)
	d = 2;
	P = 3;

	unsigned int NE23_TP, Nn23_TP, Ne23_TP;

	Nn23_TP = pow(P+1,d);
	NE23_TP = pow(P,d);
	Ne23_TP = 4;

	unsigned int connect23_TP[72] = {  0,  1,  5,  4,  0,  0,  0,  0,
	                                   1,  2,  6,  5,  0,  0,  0,  0,
	                                   2,  3,  7,  6,  0,  0,  0,  0,
	                                   4,  5,  9,  8,  0,  0,  0,  0,
	                                   5,  6, 10,  9,  0,  0,  0,  0,
	                                   6,  7, 11, 10,  0,  0,  0,  0,
	                                   8,  9, 13, 12,  0,  0,  0,  0,
	                                   9, 10, 14, 13,  0,  0,  0,  0,
	                                  10, 11, 15, 14,  0,  0,  0,  0},
	             connectE23_TP[16] = {  0,  4,  8, 12,
	                                    3,  7, 11, 15,
	                                    0,  1,  2,  3,
	                                   12, 13, 14, 15},
	             types23_TP[9];
	double rst23_TP[Nn23_TP*d];

	for (i = 0; i < NE23_TP; i++)
		types23_TP[i] = 9;

	row = 0;
	for (j = 0, jMax = P+1; j < jMax; j++) {
	for (i = 0, iMax = P+1; i < iMax; i++) {
		for (dim = 0; dim < d; dim++) {
			if      (dim == 0) Indr = i;
			else if (dim == 1) Indr = j;
			rst23_TP[dim*Nn23_TP+row] = rst13_TP[Indr];
		}
		row++;
	}}

	plotting_element_info(&rst,&connect,&types,&connectE,&Nn,&NE,P,QUAD); // free

	pass = 0;
	if (array_norm_diff_d(Nn23_TP*d,rst23_TP,rst,"Inf")                < EPS &&
	    array_norm_diff_ui(NE23_TP*8,connect23_TP,connect,"Inf")       < EPS &&
	    array_norm_diff_ui(NE23_TP,types23_TP,types,"Inf")             < EPS &&
	    array_norm_diff_ui(Ne23_TP*(P+1),connectE23_TP,connectE,"Inf") < EPS &&
	    Nn == Nn23_TP && NE == NE23_TP)
			pass = 1;

	test_print2(pass,"                      (d2, P3, TP   ):");

	free(rst);
	free(connect);
	free(types);
	free(connectE);

	/*
	 *		Expected Output (d = 3):
	 *
	 *			P = 3:
	 *				rst = [ See below ]
	 *				Nn  = 64
	 *
	 *				Nodes:
	 *				Extension of 2d.
	 *
	 *				connect:
	 *				[ See below ]
	 *
	 *				types:
	 *				[ 12 12 ... (27x) ... 12]
	 */

	// TP (d = 3, P = 3)
	d = 3;
	P = 3;

	unsigned int NE33_TP, Nn33_TP, Ne33_TP;

	Nn33_TP = pow(P+1,d);
	NE33_TP = pow(P,d);
	Ne33_TP = 12;

	unsigned int connect33_TP[216] = {  0,  1,  5,  4, 16, 17, 21, 20,
	                                    1,  2,  6,  5, 17, 18, 22, 21,
	                                    2,  3,  7,  6, 18, 19, 23, 22,
	                                    4,  5,  9,  8, 20, 21, 25, 24,
	                                    5,  6, 10,  9, 21, 22, 26, 25,
	                                    6,  7, 11, 10, 22, 23, 27, 26,
	                                    8,  9, 13, 12, 24, 25, 29, 28,
	                                    9, 10, 14, 13, 25, 26, 30, 29,
	                                   10, 11, 15, 14, 26, 27, 31, 30,
	                                   16, 17, 21, 20, 32, 33, 37, 36,
	                                   17, 18, 22, 21, 33, 34, 38, 37,
	                                   18, 19, 23, 22, 34, 35, 39, 38,
	                                   20, 21, 25, 24, 36, 37, 41, 40,
	                                   21, 22, 26, 25, 37, 38, 42, 41,
	                                   22, 23, 27, 26, 38, 39, 43, 42,
	                                   24, 25, 29, 28, 40, 41, 45, 44,
	                                   25, 26, 30, 29, 41, 42, 46, 45,
	                                   26, 27, 31, 30, 42, 43, 47, 46,
	                                   32, 33, 37, 36, 48, 49, 53, 52,
	                                   33, 34, 38, 37, 49, 50, 54, 53,
	                                   34, 35, 39, 38, 50, 51, 55, 54,
	                                   36, 37, 41, 40, 52, 53, 57, 56,
	                                   37, 38, 42, 41, 53, 54, 58, 57,
	                                   38, 39, 43, 42, 54, 55, 59, 58,
	                                   40, 41, 45, 44, 56, 57, 61, 60,
	                                   41, 42, 46, 45, 57, 58, 62, 61,
	                                   42, 43, 47, 46, 58, 59, 63, 62},
	             connectE33_TP[48] = {  0,  1,  2,  3,
	                                   12, 13, 14, 15,
	                                   48, 49, 50, 51,
	                                   60, 61, 62, 63,
	                                    0,  4,  8, 12,
	                                    3,  7, 11, 15,
	                                   48, 52, 56, 60,
	                                   51, 55, 59, 63,
	                                    0, 16, 32, 48,
	                                    3, 19, 35, 51,
	                                   12, 28, 44, 60,
	                                   15, 31, 47, 63},
	             types33_TP[27];
	double rst33_TP[Nn33_TP*d];

	for (i = 0; i < NE33_TP; i++)
		types33_TP[i] = 12;

	row = 0;
	for (k = 0, kMax = P+1; k < kMax; k++) {
	for (j = 0, jMax = P+1; j < jMax; j++) {
	for (i = 0, iMax = P+1; i < iMax; i++) {
		for (dim = 0; dim < d; dim++) {
			if      (dim == 0) Indr = i;
			else if (dim == 1) Indr = j;
			else if (dim == 2) Indr = k;
			rst33_TP[dim*Nn33_TP+row] = rst13_TP[Indr];
		}
		row++;
	}}}

	plotting_element_info(&rst,&connect,&types,&connectE,&Nn,&NE,P,HEX); // free

	pass = 0;
	if (array_norm_diff_d(Nn33_TP*d,rst33_TP,rst,"Inf")                < EPS &&
	    array_norm_diff_ui(NE33_TP*8,connect33_TP,connect,"Inf")       < EPS &&
	    array_norm_diff_ui(NE33_TP,types33_TP,types,"Inf")             < EPS &&
	    array_norm_diff_ui(Ne33_TP*(P+1),connectE33_TP,connectE,"Inf") < EPS &&
	    Nn == Nn33_TP && NE == NE33_TP)
			pass = 1;

	test_print2(pass,"                      (d3, P3, TP   ):");

	free(rst);
	free(connect);
	free(types);
	free(connectE);

	/*	Simplex
	 *
	 *		Input:
	 *
	 *			P, d
	 *
	 *		Expected Output (d = 2):
	 *
	 *			P = 3:
	 *				rst = [ See below ]
	 *				Nn  = 10
	 *
	 *				Nodes:
	 *				      9
	 *				     / \
	 *				    7 - 8
	 *				   / \ / \
	 *				  4 - 5 - 6
	 *				 / \ / \ / \
	 *				0 - 1 - 2 - 3
	 *
	 *				connect:
	 *				[  0  1  4  0  0  0  0  0;
	 *				   1  2  5  0  0  0  0  0;
	 *				   2  3  6  0  0  0  0  0;
	 *				   4  5  7  0  0  0  0  0;
	 *				   5  6  8  0  0  0  0  0;
	 *				   7  8  9  0  0  0  0  0;
	 *				   1  4  5  0  0  0  0  0;
	 *				   2  5  6  0  0  0  0  0;
	 *				   5  7  8  0  0  0  0  0]
	 *
	 *				types:
	 *				[  5  5  5  5  5  5  5  5  5]
	 */

	// SI (d = 2, P = 3)
	d = 2;
	P = 3;

	unsigned int NE23_SI, Nn23_SI, Ne23_SI;

	Nn23_SI = 1.0/2.0*(P+1)*(P+2);
	NE23_SI = 0;
	for (i = 0; i < P; i++)
		NE23_SI += 2*i+1;
	Ne23_SI = 3;

	unsigned int connect23_SI[72] = {  0,  1,  4,  0,  0,  0,  0,  0,
	                                   1,  2,  5,  0,  0,  0,  0,  0,
	                                   2,  3,  6,  0,  0,  0,  0,  0,
	                                   4,  5,  7,  0,  0,  0,  0,  0,
	                                   5,  6,  8,  0,  0,  0,  0,  0,
	                                   7,  8,  9,  0,  0,  0,  0,  0,
	                                   1,  4,  5,  0,  0,  0,  0,  0,
	                                   2,  5,  6,  0,  0,  0,  0,  0,
	                                   5,  7,  8,  0,  0,  0,  0,  0},
	             connectE23_SI[12] = { 3,  6,  8,  9,
	                                   0,  4,  7,  9,
	                                   0,  1,  2,  3},
	             types23_SI[9];
	double rst23_SI[20] = { -1.0    , -1.0/sqrt(3.0),
	                        -1.0/3.0, -1.0/sqrt(3.0),
	                         1.0/3.0, -1.0/sqrt(3.0),
	                         1.0    , -1.0/sqrt(3.0),
	                        -2.0/3.0,  0.0          ,
	                         0.0    ,  0.0          ,
	                         2.0/3.0,  0.0          ,
	                        -1.0/3.0,  1.0/sqrt(3.0),
	                         1.0/3.0,  1.0/sqrt(3.0),
	                         0.0    ,  2.0/sqrt(3.0)};

	for (i = 0; i < NE23_SI; i++)
		types23_SI[i] = 5;

	// Convert to column-major ordering
	mkl_dimatcopy('R','T',Nn23_SI,d,1.0,rst23_SI,d,Nn23_SI);

	plotting_element_info(&rst,&connect,&types,&connectE,&Nn,&NE,P,TRI); // free

	pass = 0;
	if (array_norm_diff_d(Nn23_SI*d,rst23_SI,rst,"Inf")                < EPS &&
	    array_norm_diff_ui(NE23_SI*8,connect23_SI,connect,"Inf")       < EPS &&
	    array_norm_diff_ui(NE23_SI,types23_SI,types,"Inf")             < EPS &&
	    array_norm_diff_ui(Ne23_SI*(P+1),connectE23_SI,connectE,"Inf") < EPS &&
	    Nn == Nn23_SI && NE == NE23_SI)
			pass = 1;

	test_print2(pass,"                      (d2, P3, SI   ):");

	free(rst);
	free(connect);
	free(types);
	free(connectE);


	/*
	 *		Expected Output (d = 3):
	 *
	 *			P = 3:
	 *				rst = [ See below ]
	 *				Nn  = 20
	 *
	 *				Nodes:
	 *				      9
	 *				     / \
	 *				    7 - 8            15
	 *				   / \ / \          /  \
	 *				  4 - 5 - 6       13 -- 14        18
	 *				 / \ / \ / \     /  \  /  \      /  \
	 *				0 - 1 - 2 - 3   10 - 11 - 12   16 - 17   19
	 *
	 *				Layer 0 (bot)   Layer 1        Layer 2   Layer 3
	 *
	 *
	 *				connect:
	 *				[ See Below ]
	 *
	 *				types:
	 *				[ 10 (6x) 14 (6x) 10 (1x) 10 (3x) 14 (2x) 10 (1x) ]
	 */

	// SI (d = 3, P = 3)
	d = 3;
	P = 3;

	unsigned int NE33_SI, Nn33_SI, Ne33_SI;

	Nn33_SI = 1.0/6.0*(P+1)*(P+2)*(P+3);
	NE33_SI = 19;
	Ne33_SI = 6;

	unsigned int connect33_SI[152] = {  0,  1,  4, 10,  0,  0,  0,  0,
	                                    1,  2,  5, 11,  0,  0,  0,  0,
	                                    2,  3,  6, 12,  0,  0,  0,  0,
	                                    4,  5,  7, 13,  0,  0,  0,  0,
	                                    5,  6,  8, 14,  0,  0,  0,  0,
	                                    7,  8,  9, 15,  0,  0,  0,  0,
	                                    4,  5, 11, 10,  1,  0,  0,  0,
	                                    4,  5, 11, 10, 13,  0,  0,  0,
	                                    5,  6, 12, 11,  2,  0,  0,  0,
	                                    5,  6, 12, 11, 14,  0,  0,  0,
	                                    7,  8, 14, 13,  5,  0,  0,  0,
	                                    7,  8, 14, 13, 15,  0,  0,  0,
	                                    5, 11, 13, 14,  0,  0,  0,  0,
	                                   10, 11, 13, 16,  0,  0,  0,  0,
	                                   11, 12, 14, 17,  0,  0,  0,  0,
	                                   13, 14, 15, 18,  0,  0,  0,  0,
	                                   13, 14, 17, 16, 11,  0,  0,  0,
	                                   13, 14, 17, 16, 18,  0,  0,  0,
	                                   16, 17, 18, 19,  0,  0,  0,  0},
	             connectE33_SI[24] = {  0,  1,  2,  3,
	                                    0,  4,  7,  9,
	                                    0, 10, 16, 19,
	                                    3,  6,  8,  9,
	                                    3, 12, 17, 19,
	                                    9, 15, 18, 19},
	             types33_SI[NE33_SI];
	double rst33_SI[60] = { -1.0    , -1.0/sqrt(3.0)      , -1.0/sqrt(6.0)      ,
	                        -1.0/3.0, -1.0/sqrt(3.0)      , -1.0/sqrt(6.0)      ,
	                         1.0/3.0, -1.0/sqrt(3.0)      , -1.0/sqrt(6.0)      ,
	                         1.0    , -1.0/sqrt(3.0)      , -1.0/sqrt(6.0)      ,
	                        -2.0/3.0,  0.0                , -1.0/sqrt(6.0)      ,
	                         0.0    ,  0.0                , -1.0/sqrt(6.0)      ,
	                         2.0/3.0,  0.0                , -1.0/sqrt(6.0)      ,
	                        -1.0/3.0,  1.0/sqrt(3.0)      , -1.0/sqrt(6.0)      ,
	                         1.0/3.0,  1.0/sqrt(3.0)      , -1.0/sqrt(6.0)      ,
	                         0.0    ,  2.0/sqrt(3.0)      , -1.0/sqrt(6.0)      ,
	                        -2.0/3.0, -2.0/(3.0*sqrt(3.0)),  1.0/(3.0*sqrt(6.0)),
	                         0.0    , -2.0/(3.0*sqrt(3.0)),  1.0/(3.0*sqrt(6.0)),
	                         2.0/3.0, -2.0/(3.0*sqrt(3.0)),  1.0/(3.0*sqrt(6.0)),
	                        -1.0/3.0,  1.0/(3.0*sqrt(3.0)),  1.0/(3.0*sqrt(6.0)),
	                         1.0/3.0,  1.0/(3.0*sqrt(3.0)),  1.0/(3.0*sqrt(6.0)),
	                         0.0    ,  4.0/(3.0*sqrt(3.0)),  1.0/(3.0*sqrt(6.0)),
	                        -1.0/3.0, -1.0/(3.0*sqrt(3.0)),  5.0/(3.0*sqrt(6.0)),
	                         1.0/3.0, -1.0/(3.0*sqrt(3.0)),  5.0/(3.0*sqrt(6.0)),
	                         0.0    ,  2.0/(3.0*sqrt(3.0)),  5.0/(3.0*sqrt(6.0)),
	                         0.0    ,  0.0                ,  3.0/sqrt(6.0)      };

	j = 0;
	for ( i = 0; i < 6; i++) types33_SI[j++] = 10;
	for ( i = 0; i < 6; i++) types33_SI[j++] = 14;
	for ( i = 0; i < 1; i++) types33_SI[j++] = 10;
	for ( i = 0; i < 3; i++) types33_SI[j++] = 10;
	for ( i = 0; i < 2; i++) types33_SI[j++] = 14;
	for ( i = 0; i < 1; i++) types33_SI[j++] = 10;

	// Convert to column-major ordering
	mkl_dimatcopy('R','T',Nn33_SI,d,1.0,rst33_SI,d,Nn33_SI);

	plotting_element_info(&rst,&connect,&types,&connectE,&Nn,&NE,P,TET); // free

	pass = 0;
	if (array_norm_diff_d(Nn33_SI*d,rst33_SI,rst,"Inf")                < EPS &&
	    array_norm_diff_ui(NE33_SI*8,connect33_SI,connect,"Inf")       < EPS &&
	    array_norm_diff_ui(NE33_SI,types33_SI,types,"Inf")             < EPS &&
	    array_norm_diff_ui(Ne33_SI*(P+1),connectE33_SI,connectE,"Inf") < EPS &&
	    Nn == Nn33_SI && NE == NE33_SI)
			pass = 1;

	test_print2(pass,"                      (d3, P3, SI   ):");

	free(rst);
	free(connect);
	free(types);
	free(connectE);

	/*		Wedge
	 *
	 *		Input:
	 *
	 *			P
	 *
	 *		Expected Output:
	 *
	 *			P = 3:
	 *				rst = [ See below ]
	 *				Nn  = 40
	 *
	 *				Nodes:
	 *				         09                     19                     29                     39
	 *				        /  \                   /  \                   /  \                   /  \
	 *				      07 -- 08               17 -- 18               27 -- 28               37 -- 38
	 *				     /  \  /  \             /  \  /  \             /  \  /  \             /  \  /  \
	 *				   04 -- 05 -- 06         14 -- 15 -- 16         24 -- 25 -- 26         34 -- 35 -- 36
	 *				  /  \  /  \  /  \       /  \  /  \  /  \       /  \  /  \  /  \       /  \  /  \  /  \
	 *				00 -- 01 -- 02 -- 03   10 -- 11 -- 12 -- 13   20 -- 21 -- 22 -- 23   30 -- 31 -- 32 -- 33
	 *
	 *				Layer 0 (bot)          Layer 1                Layer 2                Layer 3
	 *
	 *				connect:
	 *				[ See Below ]
	 *
	 *				types:
	 *				[ 13 (27x) ]
	 */

	// WEDGE (P = 3)
	d = 3;
	P = 3;

	unsigned int NE3_WEDGE, Nn3_WEDGE, Ne3_WEDGE;

	Nn3_WEDGE = (P+1)*1.0/2.0*(P+1)*(P+2);
	NE3_WEDGE = 27;
	Ne3_WEDGE = 9;

	unsigned int connect3_WEDGE[216] = {  0,  1,  4, 10, 11, 14,  0,  0,
	                                      1,  2,  5, 11, 12, 15,  0,  0,
	                                      2,  3,  6, 12, 13, 16,  0,  0,
	                                      4,  5,  7, 14, 15, 17,  0,  0,
	                                      5,  6,  8, 15, 16, 18,  0,  0,
	                                      7,  8,  9, 17, 18, 19,  0,  0,
	                                      1,  4,  5, 11, 14, 15,  0,  0,
	                                      2,  5,  6, 12, 15, 16,  0,  0,
	                                      5,  7,  8, 15, 17, 18,  0,  0,
	                                     10, 11, 14, 20, 21, 24,  0,  0,
	                                     11, 12, 15, 21, 22, 25,  0,  0,
	                                     12, 13, 16, 22, 23, 26,  0,  0,
	                                     14, 15, 17, 24, 25, 27,  0,  0,
	                                     15, 16, 18, 25, 26, 28,  0,  0,
	                                     17, 18, 19, 27, 28, 29,  0,  0,
	                                     11, 14, 15, 21, 24, 25,  0,  0,
	                                     12, 15, 16, 22, 25, 26,  0,  0,
	                                     15, 17, 18, 25, 27, 28,  0,  0,
	                                     20, 21, 24, 30, 31, 34,  0,  0,
	                                     21, 22, 25, 31, 32, 35,  0,  0,
	                                     22, 23, 26, 32, 33, 36,  0,  0,
	                                     24, 25, 27, 34, 35, 37,  0,  0,
	                                     25, 26, 28, 35, 36, 38,  0,  0,
	                                     27, 28, 29, 37, 38, 39,  0,  0,
	                                     21, 24, 25, 31, 34, 35,  0,  0,
	                                     22, 25, 26, 32, 35, 36,  0,  0,
	                                     25, 27, 28, 35, 37, 38,  0,  0},
	             connectE3_WEDGE[36] = {  3,  6,  8,  9,
	                                      0,  4,  7,  9,
	                                      0,  1,  2,  3,
	                                     33, 36, 38, 39,
	                                     30, 34, 37, 39,
	                                     30, 31, 32, 33,
	                                      0, 10, 20, 30,
	                                      1, 11, 21, 31,
	                                      2, 12, 22, 32},
	             types3_WEDGE[NE3_WEDGE];
	double rst3_WEDGE[120] = { -1.0    , -1.0/sqrt(3.0)      , -1.0    ,
	                           -1.0/3.0, -1.0/sqrt(3.0)      , -1.0    ,
	                            1.0/3.0, -1.0/sqrt(3.0)      , -1.0    ,
	                            1.0    , -1.0/sqrt(3.0)      , -1.0    ,
	                           -2.0/3.0,  0.0                , -1.0    ,
	                            0.0    ,  0.0                , -1.0    ,
	                            2.0/3.0,  0.0                , -1.0    ,
	                           -1.0/3.0,  1.0/sqrt(3.0)      , -1.0    ,
	                            1.0/3.0,  1.0/sqrt(3.0)      , -1.0    ,
	                            0.0    ,  2.0/sqrt(3.0)      , -1.0    ,
	                           -1.0    , -1.0/sqrt(3.0)      , -1.0/3.0,
	                           -1.0/3.0, -1.0/sqrt(3.0)      , -1.0/3.0,
	                            1.0/3.0, -1.0/sqrt(3.0)      , -1.0/3.0,
	                            1.0    , -1.0/sqrt(3.0)      , -1.0/3.0,
	                           -2.0/3.0,  0.0                , -1.0/3.0,
	                            0.0    ,  0.0                , -1.0/3.0,
	                            2.0/3.0,  0.0                , -1.0/3.0,
	                           -1.0/3.0,  1.0/sqrt(3.0)      , -1.0/3.0,
	                            1.0/3.0,  1.0/sqrt(3.0)      , -1.0/3.0,
	                            0.0    ,  2.0/sqrt(3.0)      , -1.0/3.0,
	                           -1.0    , -1.0/sqrt(3.0)      ,  1.0/3.0,
	                           -1.0/3.0, -1.0/sqrt(3.0)      ,  1.0/3.0,
	                            1.0/3.0, -1.0/sqrt(3.0)      ,  1.0/3.0,
	                            1.0    , -1.0/sqrt(3.0)      ,  1.0/3.0,
	                            -2.0/3.0,  0.0                ,  1.0/3.0,
	                            0.0    ,  0.0                ,  1.0/3.0,
	                            2.0/3.0,  0.0                ,  1.0/3.0,
	                           -1.0/3.0,  1.0/sqrt(3.0)      ,  1.0/3.0,
	                            1.0/3.0,  1.0/sqrt(3.0)      ,  1.0/3.0,
	                            0.0    ,  2.0/sqrt(3.0)      ,  1.0/3.0,
	                           -1.0    , -1.0/sqrt(3.0)      ,  1.0    ,
	                           -1.0/3.0, -1.0/sqrt(3.0)      ,  1.0    ,
	                            1.0/3.0, -1.0/sqrt(3.0)      ,  1.0    ,
	                            1.0    , -1.0/sqrt(3.0)      ,  1.0    ,
	                           -2.0/3.0,  0.0                ,  1.0    ,
	                            0.0    ,  0.0                ,  1.0    ,
	                            2.0/3.0,  0.0                ,  1.0    ,
	                           -1.0/3.0,  1.0/sqrt(3.0)      ,  1.0    ,
	                            1.0/3.0,  1.0/sqrt(3.0)      ,  1.0    ,
	                            0.0    ,  2.0/sqrt(3.0)      ,  1.0    };

	for ( i = 0; i < NE3_WEDGE; i++)
		types3_WEDGE[i] = 13;

	// Convert to column-major ordering
	mkl_dimatcopy('R','T',Nn3_WEDGE,d,1.0,rst3_WEDGE,d,Nn3_WEDGE);

	plotting_element_info(&rst,&connect,&types,&connectE,&Nn,&NE,P,WEDGE); // free

	pass = 0;
	if (array_norm_diff_d(Nn3_WEDGE*d,rst3_WEDGE,rst,"Inf")                < EPS &&
	    array_norm_diff_ui(NE3_WEDGE*8,connect3_WEDGE,connect,"Inf")       < EPS &&
	    array_norm_diff_ui(NE3_WEDGE,types3_WEDGE,types,"Inf")             < EPS &&
	    array_norm_diff_ui(Ne3_WEDGE*(P+1),connectE3_WEDGE,connectE,"Inf") < EPS &&
	    Nn == Nn3_WEDGE && NE == NE3_WEDGE)
			pass = 1;

	test_print2(pass,"                      (    P3, WEDGE):");

	free(rst);
	free(connect);
	free(types);
	free(connectE);

	/*		Pyramid
	 *
	 *		Input:
	 *
	 *			P
	 *
	 *		Expected Output:
	 *
	 *			P = 3:
	 *				rst = [ See below ]
	 *				Nn  = 30
	 *
	 *				Nodes:
	 *
	 *				012 - 013 - 014 - 015
	 *				 |     |     |     |
	 *				008 - 009 - 010 - 011   022 - 023 - 024
	 *				 |     |     |     |     |     |     |
	 *				004 - 005 - 006 - 007   019 - 020 - 021   027 - 028
	 *				 |     |     |     |     |     |     |     |     |
	 *				000 - 001 - 002 - 003   016 - 017 - 018   025 - 026   029
	 *
	 *				Layer 0 (bot)           Layer 1           Layer 2     Layer 3
	 *
	 *				connect:
	 *				[ See Below ]
	 *
	 *				types:
	 *				[ 14 (9x) 10 (2*2*3x) 14 (4x) 14 (4x) 10 (2*1*2x) 14 (1x) 14 (1x)]
	 */

	// PYR (P = 3)
	d = 3;
	P = 3;

	unsigned int NE3_PYR, Nn3_PYR, Ne3_PYR;

	Nn3_PYR = 16+9+4+1;
	NE3_PYR = 35;
	Ne3_PYR = 8;

	unsigned int connect3_PYR[280] = {  0,  1,  5,  4, 16,  0,  0,  0,
	                                    1,  2,  6,  5, 17,  0,  0,  0,
	                                    2,  3,  7,  6, 18,  0,  0,  0,
	                                    4,  5,  9,  8, 19,  0,  0,  0,
	                                    5,  6, 10,  9, 20,  0,  0,  0,
	                                    6,  7, 11, 10, 21,  0,  0,  0,
	                                    8,  9, 13, 12, 22,  0,  0,  0,
	                                    9, 10, 14, 13, 23,  0,  0,  0,
	                                   10, 11, 15, 14, 24,  0,  0,  0,
	                                    1,  5, 16, 17,  0,  0,  0,  0,
	                                    2,  6, 17, 18,  0,  0,  0,  0,
	                                    5,  9, 19, 20,  0,  0,  0,  0,
	                                    6, 10, 20, 21,  0,  0,  0,  0,
	                                    9, 13, 22, 23,  0,  0,  0,  0,
	                                   10, 14, 23, 24,  0,  0,  0,  0,
	                                    4,  5, 16, 19,  0,  0,  0,  0,
	                                    8,  9, 19, 22,  0,  0,  0,  0,
	                                    5,  6, 17, 20,  0,  0,  0,  0,
	                                    9, 10, 20, 23,  0,  0,  0,  0,
	                                    6,  7, 18, 21,  0,  0,  0,  0,
	                                   10, 11, 21, 24,  0,  0,  0,  0,
	                                   16, 17, 20, 19,  5,  0,  0,  0,
	                                   17, 18, 21, 20,  6,  0,  0,  0,
	                                   19, 20, 23, 22,  9,  0,  0,  0,
	                                   20, 21, 24, 23, 10,  0,  0,  0,
	                                   16, 17, 20, 19, 25,  0,  0,  0,
	                                   17, 18, 21, 20, 26,  0,  0,  0,
	                                   19, 20, 23, 22, 27,  0,  0,  0,
	                                   20, 21, 24, 23, 28,  0,  0,  0,
	                                   17, 20, 25, 26,  0,  0,  0,  0,
	                                   20, 23, 27, 28,  0,  0,  0,  0,
	                                   19, 20, 25, 27,  0,  0,  0,  0,
	                                   20, 21, 26, 28,  0,  0,  0,  0,
	                                   25, 26, 28, 27, 20,  0,  0,  0,
	                                   25, 26, 28, 27, 29,  0,  0,  0},
	             connectE3_PYR[32] = {  0,  1,  2,  3,
	                                   12, 13, 14, 15,
	                                    0,  4,  8, 12,
	                                    3,  7, 11, 15,
	                                    0, 16, 25, 29,
	                                    3, 18, 26, 29,
	                                   12, 22, 27, 29,
	                                   15, 24, 28, 29},
	             types3_PYR[NE3_PYR];
	double rst3_PYR[90] = { -1.0    , -1.0    , -1.0/5.0 *sqrt(2.0),
	                        -1.0/3.0, -1.0    , -1.0/5.0 *sqrt(2.0),
	                         1.0/3.0, -1.0    , -1.0/5.0 *sqrt(2.0),
	                         1.0    , -1.0    , -1.0/5.0 *sqrt(2.0),
	                        -1.0    , -1.0/3.0, -1.0/5.0 *sqrt(2.0),
	                        -1.0/3.0, -1.0/3.0, -1.0/5.0 *sqrt(2.0),
	                         1.0/3.0, -1.0/3.0, -1.0/5.0 *sqrt(2.0),
	                         1.0    , -1.0/3.0, -1.0/5.0 *sqrt(2.0),
	                        -1.0    ,  1.0/3.0, -1.0/5.0 *sqrt(2.0),
	                        -1.0/3.0,  1.0/3.0, -1.0/5.0 *sqrt(2.0),
	                         1.0/3.0,  1.0/3.0, -1.0/5.0 *sqrt(2.0),
	                         1.0    ,  1.0/3.0, -1.0/5.0 *sqrt(2.0),
	                        -1.0    ,  1.0    , -1.0/5.0 *sqrt(2.0),
	                        -1.0/3.0,  1.0    , -1.0/5.0 *sqrt(2.0),
	                         1.0/3.0,  1.0    , -1.0/5.0 *sqrt(2.0),
	                         1.0    ,  1.0    , -1.0/5.0 *sqrt(2.0),
	                        -2.0/3.0, -2.0/3.0,  2.0/15.0*sqrt(2.0),
	                         0.0    , -2.0/3.0,  2.0/15.0*sqrt(2.0),
	                         2.0/3.0, -2.0/3.0,  2.0/15.0*sqrt(2.0),
	                        -2.0/3.0,  0.0    ,  2.0/15.0*sqrt(2.0),
	                         0.0    ,  0.0    ,  2.0/15.0*sqrt(2.0),
	                         2.0/3.0,  0.0    ,  2.0/15.0*sqrt(2.0),
	                        -2.0/3.0,  2.0/3.0,  2.0/15.0*sqrt(2.0),
	                         0.0    ,  2.0/3.0,  2.0/15.0*sqrt(2.0),
	                         2.0/3.0,  2.0/3.0,  2.0/15.0*sqrt(2.0),
	                        -1.0/3.0, -1.0/3.0,  7.0/15.0*sqrt(2.0),
	                         1.0/3.0, -1.0/3.0,  7.0/15.0*sqrt(2.0),
	                        -1.0/3.0,  1.0/3.0,  7.0/15.0*sqrt(2.0),
	                         1.0/3.0,  1.0/3.0,  7.0/15.0*sqrt(2.0),
	                         0.0    ,  0.0    ,  4.0/5.0 *sqrt(2.0)};

	j = 0;
	for ( i = 0; i < 9 ; i++) types3_PYR[j++] = 14;
	for ( i = 0; i < 12; i++) types3_PYR[j++] = 10;
	for ( i = 0; i < 4 ; i++) types3_PYR[j++] = 14;
	for ( i = 0; i < 4 ; i++) types3_PYR[j++] = 14;
	for ( i = 0; i < 4 ; i++) types3_PYR[j++] = 10;
	for ( i = 0; i < 1 ; i++) types3_PYR[j++] = 14;
	for ( i = 0; i < 1 ; i++) types3_PYR[j++] = 14;

	// Convert to column-major ordering
	mkl_dimatcopy('R','T',Nn3_PYR,d,1.0,rst3_PYR,d,Nn3_PYR);

	plotting_element_info(&rst,&connect,&types,&connectE,&Nn,&NE,P,PYR); // free

	pass = 0;
	if (array_norm_diff_d(Nn3_PYR*d,rst3_PYR,rst,"Inf")                < EPS &&
	    array_norm_diff_ui(NE3_PYR*8,connect3_PYR,connect,"Inf")       < EPS &&
	    array_norm_diff_ui(NE3_PYR,types3_PYR,types,"Inf")             < EPS &&
	    array_norm_diff_ui(Ne3_PYR*(P+1),connectE3_PYR,connectE,"Inf") < EPS &&
	    Nn == Nn3_PYR && NE == NE3_PYR)
			pass = 1;

	test_print2(pass,"                      (    P3, PYR  ):");

	free(rst);
	free(connect);
	free(types);
	free(connectE);
}
