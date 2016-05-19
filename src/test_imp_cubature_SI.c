// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>

#include "test.h"
#include "functions.h"
#include "parameters.h"

/*
 *	Purpose:
 *		Test correctness of implementation of cubature_TRI and cubature_TET.
 *
 *	Comments:
 *		Other than the grouping based on 1/3 symmetries, the ordering of the tetrahedral nodes is arbitrary.
 *
 *	Notation:
 *
 *	References:
 */

void test_imp_cubature_SI(void)
{
	unsigned int pass;

	/*
	 *	cubature_TRI:
	 *
	 *		Input:
	 *
	 *			P, NodeType.
	 *
	 *		Expected output:
	 *
	 *			P = 3, AO:
	 *				rst = [ See below ]
	 *				w   = N/A
	 *				symms = [ 3 3 3 1 ]
	 *				Ns  = 4
	 *				Nn  = 10
	 *
	 *			P = 3, WS:
	 *				rst = [ See below ]
	 *				w   = [ See below ]
	 *				symms = [ 3 3 3 1 ]
	 *				Ns  = 4
	 *				Nn  = 10
	 *
	 *			P = 6, WV:
	 *				rst = [ See below ]
	 *				w   = [ See below ]
	 *				symms = [ 3 3 3 3 ]
	 *				Ns  = 4
	 *				Nn  = 12
	 */

	unsigned int Nn, Ns, P, d;
	unsigned int *symms;
	double *rst, *w;

	d = 2;

	// AO (P = 3)
	P = 3;
	unsigned int Nn23_AO = 10, Ns23_AO = 4;
	unsigned int symms23_AO[4] = { 3, 3, 3, 1 };
	double rst23_AO[20] = { -0.447213595499958, -0.577350269189626,
	                         0.723606797749979, -0.098623200025929,
	                        -0.276393202250021,  0.675973469215555,
	                        -0.723606797749979, -0.098623200025929,
	                         0.447213595499958, -0.577350269189626,
	                         0.276393202250021,  0.675973469215555,
	                        -1.000000000000000, -0.577350269189626,
	                         1.000000000000000, -0.577350269189626,
	                        -0.000000000000000,  1.154700538379252,
	                        -0.000000000000000, -0.000000000000000};

	// Convert to column-major ordering
	mkl_dimatcopy('R','T',Nn23_AO,d,1.0,rst23_AO,d,Nn23_AO);

	cubature_TRI(&rst,&w,&symms,&Nn,&Ns,0,P,d,"AO");

	pass = 0;
	if (array_norm_diff_d(Nn23_AO*d,rst23_AO,rst,"Inf")    < EPS &&
	    array_norm_diff_ui(Ns23_AO,symms23_AO,symms,"Inf") < EPS &&
	    Nn23_AO == Nn && Ns23_AO == Ns)
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("cubature_SI (d2, P3, AO):                        ");
	test_print(pass);

	free(rst);
	free(symms);

	// WS (P = 3)
	P = 3;
	unsigned int Nn23_WS = 10, Ns23_WS = 4;
	unsigned int symms23_WS[4] = { 3, 3, 3, 1 };
	double rst23_WS[20] = { -0.338677036009830,  -0.455664103498571,
	                         0.563955207227339,  -0.065470865113645,
	                        -0.225278171217509,   0.521134968612215,
	                        -0.563955207227339,  -0.065470865113645,
	                         0.338677036009830,  -0.455664103498571,
	                         0.225278171217509,   0.521134968612215,
	                        -0.833307841990621,  -0.481110506891111,
	                         0.833307841990621,  -0.481110506891111,
	                        -0.000000000000000,   0.962221013782223,
	                        -0.000000000000000,  -0.000000000000001};

	double w23_WS[10] = { 0.194160145154569,
	                      0.194160145154569,
	                      0.194160145154569,
	                      0.194160145154569,
	                      0.194160145154569,
	                      0.194160145154569,
	                      0.072669080167812,
	                      0.072669080167812,
	                      0.072669080167812,
	                      0.349082696138027};

	// Convert to column-major ordering
	mkl_dimatcopy('R','T',Nn23_WS,d,1.0,rst23_WS,d,Nn23_WS);

	cubature_TRI(&rst,&w,&symms,&Nn,&Ns,1,P,d,"WS");

	pass = 0;
	if (array_norm_diff_d(Nn23_WS*d,rst23_WS,rst,"Inf")    < EPS &&
	    array_norm_diff_d(Nn23_WS,w23_WS,w,"Inf")          < EPS &&
	    array_norm_diff_ui(Ns23_WS,symms23_WS,symms,"Inf") < EPS &&
	    Nn23_WS == Nn && Ns23_WS == Ns)
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("            (d2, P3, WS):                        ");
	test_print(pass);

	free(rst), free(w);
	free(symms);

	// WV (P = 6)
	P = 6;

	unsigned int Nn26_WV = 12, Ns26_WV = 4;
	unsigned int symms26_WV[4] = { 3, 3, 3, 3 };
	double rst26_WV[24] = { -0.326150048087614,  -0.485300342687622,
	                         0.583357449276582,  -0.039804055745579,
	                        -0.257207401188967,   0.525104398433202,
	                        -0.583357449276582,  -0.039804055745579,
	                         0.326150048087614,  -0.485300342687622,
	                         0.257207401188967,   0.525104398433202,
	                        -0.810732956525493,  -0.468076890690895,
	                         0.810732956525494,  -0.468076890690895,
	                        -0.000000000000000,   0.936153781381790,
	                        -0.252139764487269,  -0.145572960900134,
	                         0.252139764487269,  -0.145572960900134,
	                         0.000000000000000,   0.291145921800267};

	double w26_WV[12] = { 0.143502272432754,
	                      0.143502272432754,
	                      0.143502272432754,
	                      0.143502272432754,
	                      0.143502272432754,
	                      0.143502272432754,
	                      0.088065961139281,
	                      0.088065961139281,
	                      0.088065961139281,
	                      0.202279763184837,
	                      0.202279763184837,
	                      0.202279763184837};


	// Convert to column-major ordering
	mkl_dimatcopy('R','T',Nn26_WV,d,1.0,rst26_WV,d,Nn26_WV);

	cubature_TRI(&rst,&w,&symms,&Nn,&Ns,1,P,d,"WV");

	pass = 0;
	if (array_norm_diff_d(Nn26_WV*d,rst26_WV,rst,"Inf")    < EPS &&
	    array_norm_diff_d(Nn26_WV,w26_WV,w,"Inf")          < EPS &&
	    array_norm_diff_ui(Ns26_WV,symms26_WV,symms,"Inf") < EPS &&
	    Nn26_WV == Nn && Ns26_WV == Ns)
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("            (d2, P6, WV):                        ");
	test_print(pass);

	free(rst), free(w);
	free(symms);

	/*
	 *	cubature_TET:
	 *
	 *		Input:
	 *
	 *			P, NodeType.
	 *
	 *		Expected output:
	 *
	 *			P = 2, AO:
	 *				rst = [ See below ]
	 *				w   = N/A
	 *				symms = [ 3 3 3 1 ]
	 *				Ns  = 4
	 *				Nn  = 10
	 *
	 *			P = 2, SH:
	 *				rst = [ See below ]
	 *				w   = [ See below ]
	 *				symms = [ 3 3 3 1 ]
	 *				Ns  = 4
	 *				Nn  = 10
	 *
	 *			P = 4, WV:
	 *				rst = [ See below ]
	 *				w   = [ See below ]
	 *				symms = [ 3 3 3 3 1 1 ]
	 *				Ns  = 6
	 *				Nn  = 14
	 */

	d = 3;

	// AO (P = 3)
	P = 2;
	unsigned int Nn33_AO = 10, Ns33_AO = 4;
	unsigned int symms33_AO[4] = { 3, 3, 3, 1 };
	double rst33_AO[30] = { 0.000000000000000, -0.577350269189626, -0.408248290463863,
	                        0.500000000000000,  0.288675134594813, -0.408248290463863,
	                       -0.500000000000000,  0.288675134594813, -0.408248290463863,
	                       -0.500000000000000, -0.288675134594813,  0.408248290463863,
	                        0.500000000000000, -0.288675134594813,  0.408248290463863,
	                        0.000000000000000,  0.577350269189626,  0.408248290463863,
	                       -1.000000000000000, -0.577350269189626, -0.408248290463863,
	                        1.000000000000000, -0.577350269189626, -0.408248290463863,
	                        0.000000000000000,  1.154700538379252, -0.408248290463863,
	                        0.000000000000000,  0.000000000000000,  1.224744871391589};

	// Convert to column-major ordering
	mkl_dimatcopy('R','T',Nn33_AO,d,1.0,rst33_AO,d,Nn33_AO);

	cubature_TET(&rst,&w,&symms,&Nn,&Ns,0,P,d,"AO");

	pass = 0;
	if (array_norm_diff_d(Nn33_AO*d,rst33_AO,rst,"Inf")    < EPS &&
	    array_norm_diff_ui(Ns33_AO,symms33_AO,symms,"Inf") < EPS &&
	    Nn33_AO == Nn && Ns33_AO == Ns)
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("            (d3, P2, AO):                        ");
	test_print(pass);

	free(rst);
	free(symms);

	// SH (P = 2)
	P = 2;
	unsigned int Nn33_SH = 10, Ns33_SH = 4;
	unsigned int symms33_SH[4] = { 3, 3, 3, 1 };
	double rst33_SH[30] = { 0.000000000000000, -0.360830856003254, -0.255145945141247,
	                        0.312488687768103,  0.180415428001627, -0.255145945141247,
	                       -0.312488687768103,  0.180415428001627, -0.255145945141247,
	                       -0.312488687768103, -0.180415428001627,  0.255145945141247,
	                        0.312488687768103, -0.180415428001627,  0.255145945141247,
	                        0.000000000000000,  0.360830856003254,  0.255145945141247,
	                       -0.704660393095107, -0.406835867640727, -0.287676400838671,
	                        0.704660393095107, -0.406835867640727, -0.287676400838671,
	                        0.000000000000000,  0.813671735281455, -0.287676400838671,
	                        0.000000000000000,  0.000000000000000,  0.863029202516013};

	double w33_SH[10] = { 0.127195540124294,
	                      0.127195540124294,
	                      0.127195540124294,
	                      0.127195540124294,
	                      0.127195540124294,
	                      0.127195540124294,
	                      0.044908950209075,
	                      0.044908950209075,
	                      0.044908950209075,
	                      0.044908950209075};

	// Convert to column-major ordering
	mkl_dimatcopy('R','T',Nn33_SH,d,1.0,rst33_SH,d,Nn33_SH);

	cubature_TET(&rst,&w,&symms,&Nn,&Ns,1,P,d,"SH");

	pass = 0;
	if (array_norm_diff_d(Nn33_SH*d,rst33_SH,rst,"Inf")    < EPS &&
	    array_norm_diff_d(Nn33_SH,w33_SH,w,"Inf")          < EPS &&
	    array_norm_diff_ui(Ns33_SH,symms33_SH,symms,"Inf") < EPS &&
	    Nn33_SH == Nn && Ns33_SH == Ns)
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("            (d3, P2, SH):                        ");
	test_print(pass);

	free(rst), free(w);
	free(symms);

	// WV (P = 4)
	P = 4;

	unsigned int Nn34_WV = 14, Ns34_WV = 6;
	unsigned int symms34_WV[6] = { 3, 3, 3, 3, 1, 1 };
	double rst34_WV[42] = { 0.000000000000000, -0.472263965885350, -0.333941052787584,
	                        0.408992591748701,  0.236131982942675, -0.333941052787584,
	                       -0.408992591748701,  0.236131982942675, -0.333941052787584,
	                       -0.408992591748701, -0.236131982942675,  0.333941052787584,
	                        0.408992591748701, -0.236131982942675,  0.333941052787584,
	                        0.000000000000000,  0.472263965885350,  0.333941052787584,
	                       -0.629058998756435, -0.363187382268184, -0.256812260843224,
	                        0.629058998756435, -0.363187382268184, -0.256812260843224,
	                        0.000000000000000,  0.726374764536369, -0.256812260843224,
	                        0.243543677053203,  0.140610007506098,  0.099426289810253,
	                       -0.243543677053203,  0.140610007506098,  0.099426289810253,
	                        0.000000000000000, -0.281220015012196,  0.099426289810253,
	                        0.000000000000000,  0.000000000000000,  0.770436782529672,
	                        0.000000000000000,  0.000000000000000, -0.298278869430759};

	double w34_WV[14] = { 0.040112773071971,
	                      0.040112773071971,
	                      0.040112773071971,
	                      0.040112773071971,
	                      0.040112773071971,
	                      0.040112773071971,
	                      0.069289905543486,
	                      0.069289905543486,
	                      0.069289905543486,
	                      0.106243195244073,
	                      0.106243195244073,
	                      0.106243195244073,
	                      0.069289905543486,
	                      0.106243195244073};


	// Convert to column-major ordering
	mkl_dimatcopy('R','T',Nn34_WV,d,1.0,rst34_WV,d,Nn34_WV);

	cubature_TET(&rst,&w,&symms,&Nn,&Ns,1,P,d,"WV");

	pass = 0;
	if (array_norm_diff_d(Nn34_WV*d,rst34_WV,rst,"Inf")    < EPS*10 &&
	    array_norm_diff_d(Nn34_WV,w34_WV,w,"Inf")          < EPS    &&
	    array_norm_diff_ui(Ns34_WV,symms34_WV,symms,"Inf") < EPS    &&
	    Nn34_WV == Nn && Ns34_WV == Ns)
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("            (d3, P4, WV):                        ");
	test_print(pass);

	free(rst), free(w);
	free(symms);
}
