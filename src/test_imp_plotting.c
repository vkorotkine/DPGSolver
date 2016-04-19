#include <stdlib.h>
#include <stdio.h>
#include <math.h>
//#include <string.h>

#include "test.h"
//#include "database.h"
#include "parameters.h"
#include "functions.h"

//#include "petscsys.h"

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

void test_imp_plotting(void)
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

	unsigned int i, j, k, iMax, jMax, kMax, row, dim, Indr, d, P, Nn;
	unsigned int *connect, *types;
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

	plotting_element_info(&rst,&connect,&types,&Nn,P,LINE); // free

	pass = 0;
	if (array_norm_diff_d(Nn13_TP*d,rst13_TP,rst,"Inf")          < EPS &&
		array_norm_diff_ui(NE13_TP*8,connect13_TP,connect,"Inf") < EPS &&
		array_norm_diff_ui(NE13_TP,types13_TP,types,"Inf")       < EPS &&
		Nn == Nn13_TP)
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("plotting_element_info (d1, P3, TP):              ");
	test_print(pass);

	free(rst);
	free(connect);
	free(types);


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

	unsigned int NE23_TP, Nn23_TP;

	Nn23_TP = pow(P+1,d);
	NE23_TP = pow(P,d);

	unsigned int connect23_TP[72] = {  0,  1,  5,  4,  0,  0,  0,  0,
	                                   1,  2,  6,  5,  0,  0,  0,  0,
	                                   2,  3,  7,  6,  0,  0,  0,  0,
	                                   4,  5,  9,  8,  0,  0,  0,  0,
	                                   5,  6, 10,  9,  0,  0,  0,  0,
	                                   6,  7, 11, 10,  0,  0,  0,  0,
	                                   8,  9, 13, 12,  0,  0,  0,  0,
	                                   9, 10, 14, 13,  0,  0,  0,  0,
	                                  10, 11, 15, 14,  0,  0,  0,  0},
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

	plotting_element_info(&rst,&connect,&types,&Nn,P,QUAD); // free

	pass = 0;
	if (array_norm_diff_d(Nn23_TP,rst23_TP,rst,"Inf")            < EPS &&
		array_norm_diff_ui(NE23_TP*8,connect23_TP,connect,"Inf") < EPS &&
		array_norm_diff_ui(NE23_TP,types23_TP,types,"Inf")       < EPS &&
		Nn == Nn23_TP)
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("                      (d2, P3, TP):              ");
	test_print(pass);

	free(rst);
	free(connect);
	free(types);

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

	unsigned int NE33_TP, Nn33_TP;

	Nn33_TP = pow(P+1,d);
	NE33_TP = pow(P,d);

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
			rst33_TP[dim*Nn23_TP+row] = rst13_TP[Indr];
		}
		row++;
	}}}

	plotting_element_info(&rst,&connect,&types,&Nn,P,HEX); // free

	pass = 0;
	if (array_norm_diff_d(Nn33_TP,rst33_TP,rst,"Inf") < EPS &&
		array_norm_diff_ui(NE33_TP*8,connect33_TP,connect,"Inf")     < EPS &&
		array_norm_diff_ui(NE33_TP,types33_TP,types,"Inf")     < EPS &&
		Nn == Nn33_TP)
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("                      (d3, P3, TP):              ");
	test_print(pass);

	free(rst);
	free(connect);
	free(types);

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

	unsigned int NE23_SI, Nn23_SI;

	Nn23_SI = 1.0/2.0*(P+1)*(P+2);
	NE23_SI = 0;
	for (i = 0; i < P; i++)
		NE23_SI += 2*i+1;

	unsigned int connect23_SI[72] = {  0,  1,  4,  0,  0,  0,  0,  0,
	                                   1,  2,  5,  0,  0,  0,  0,  0,
	                                   2,  3,  6,  0,  0,  0,  0,  0,
	                                   4,  5,  7,  0,  0,  0,  0,  0,
	                                   5,  6,  8,  0,  0,  0,  0,  0,
	                                   7,  8,  9,  0,  0,  0,  0,  0,
	                                   1,  4,  5,  0,  0,  0,  0,  0,
	                                   2,  5,  6,  0,  0,  0,  0,  0,
	                                   5,  7,  8,  0,  0,  0,  0,  0},
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

	plotting_element_info(&rst,&connect,&types,&Nn,P,TRI); // free

	pass = 0;
	if (array_norm_diff_d(Nn23_SI,rst23_SI,rst,"Inf")            < EPS &&
		array_norm_diff_ui(NE23_SI*8,connect23_SI,connect,"Inf") < EPS &&
		array_norm_diff_ui(NE23_SI,types23_SI,types,"Inf")       < EPS &&
		Nn == Nn23_SI)
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("                      (d2, P3, SI):              ");
	test_print(pass);

	free(rst);
	free(connect);
	free(types);


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
	 *				  4 - 5 - 6       13    14        18
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
	 *				[ 10 (6x) 14 (6x) 10 (3x) 14 (2x) 10 (1x) ]
	 */

	// SI (d = 3, P = 3)
	d = 3;
	P = 3;

	unsigned int NE33_SI, Nn33_SI;

	Nn33_SI = 1.0/6.0*(P+1)*(P+2)*(P+3);
	NE33_SI = 18;

	unsigned int connect33_SI[144] = {  0,  1,  4, 10,  0,  0,  0,  0,
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
	                                   10, 11, 13, 16,  0,  0,  0,  0,
	                                   11, 12, 14, 17,  0,  0,  0,  0,
	                                   13, 14, 15, 18,  0,  0,  0,  0,
	                                   13, 14, 17, 16, 11,  0,  0,  0,
	                                   13, 14, 17, 16, 18,  0,  0,  0,
	                                   16, 17, 18, 19,  0,  0,  0,  0},
	             types33_SI[9];
	double rst33_SI[60] = { -1.0    , -1.0/sqrt(3.0), -1.0/sqrt(6.0),
		                    -1.0/3.0, -1.0/sqrt(3.0), -1.0/sqrt(6.0),
		                     1.0/3.0, -1.0/sqrt(3.0), -1.0/sqrt(6.0),
		                     1.0    , -1.0/sqrt(3.0), -1.0/sqrt(6.0),
		                    -2.0/3.0,  0.0          , -1.0/sqrt(6.0),
		                     0.0    ,  0.0          , -1.0/sqrt(6.0),
		                     2.0/3.0,  0.0          , -1.0/sqrt(6.0),
		                    -1.0/3.0,  1.0/sqrt(3.0), -1.0/sqrt(6.0),
		                     1.0/3.0,  1.0/sqrt(3.0), -1.0/sqrt(6.0),
		                     0.0    ,  2.0/sqrt(3.0), -1.0/sqrt(6.0),
							 };
	i = 0;
	for ( ; i < 6; i++) types33_SI[i] = 10;
	for ( ; i < 6; i++) types33_SI[i] = 14;
	for ( ; i < 3; i++) types33_SI[i] = 10;
	for ( ; i < 2; i++) types33_SI[i] = 14;
	for ( ; i < 1; i++) types33_SI[i] = 10;

	plotting_element_info(&rst,&connect,&types,&Nn,P,TET); // free

	pass = 0;
	if (array_norm_diff_d(Nn33_SI,rst33_SI,rst,"Inf")            < EPS &&
		array_norm_diff_ui(NE33_SI*8,connect33_SI,connect,"Inf") < EPS &&
		array_norm_diff_ui(NE33_SI,types33_SI,types,"Inf")       < EPS &&
		Nn == Nn33_SI)
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("                      (d2, P3, SI):              ");
	test_print(pass);

	free(rst);
	free(connect);
	free(types);
}
