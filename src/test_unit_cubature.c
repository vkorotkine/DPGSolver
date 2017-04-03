// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "test_unit_cubature.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "mkl.h"

#include "Parameters.h"
#include "Test.h"

#include "test_support.h"
#include "array_norm.h"
#include "cubature.h"

void test_unit_cubature_TP(void)
{
	/*
	 *	Purpose:
	 *		Test correctness of implementation of cubature_TP.
	 *
	 *	Comments:
	 *		'symms' is only tested for 1d (ToBeModified).
	 *
	 *	Notation:
	 *
	 *	References:
	 *		GL  : http://mathworld.wolfram.com/Legendre-GaussQuadrature.html
	 *		GLL : http://mathworld.wolfram.com/LobattoQuadrature.html
	 */

	unsigned int pass;

	/*
	 *	cubature_TP (dE = 1):
	 *
	 *		Input:
	 *
	 *			P, dE, NodeType.
	 *
	 *		Expected output:
	 *
	 *			P = 3/4, GL:
	 *				rst   = [ See below ]
	 *				w     = [ See below ]
	 *				symms = [ 2 2 ]/[ 2 2 1]
	 *				Nn  = 4/5
	 *				Ns  = 2/3
	 *
	 *			P = 3/4, GLL:
	 *				rst   = [ See below ]
	 *				w     = [ See below ]
	 *				symms = [ 2 2 ]/[ 2 2 1]
	 *				Nn  = 4/5
	 *				Ns  = 2/3
	 */

	unsigned int i, j, k, row;
	unsigned int dE, Nn, Ns, P;
	unsigned int *symms;
	double *rst, *w;

	dE = 1;

	// P = 3
	P = 3;

	unsigned int N13 = P+1, Ns13 = 2;
	unsigned int symms13_GL[Ns13], symms13_GLL[Ns13];
	double rst13_GL[N13], w13_GL[N13], rst13_GLL[N13], w13_GLL[N13];

	rst13_GL[0] = -1.0/35.0*sqrt(525.0+70.0*sqrt(30.0));
	rst13_GL[1] =  1.0/35.0*sqrt(525.0+70.0*sqrt(30.0));
	rst13_GL[2] = -1.0/35.0*sqrt(525.0-70.0*sqrt(30.0));
	rst13_GL[3] =  1.0/35.0*sqrt(525.0-70.0*sqrt(30.0));

	w13_GL[0] = 1.0/36.0*(18.0-sqrt(30.0));
	w13_GL[1] = 1.0/36.0*(18.0-sqrt(30.0));
	w13_GL[2] = 1.0/36.0*(18.0+sqrt(30.0));
	w13_GL[3] = 1.0/36.0*(18.0+sqrt(30.0));

	symms13_GL[0] = 2;
	symms13_GL[1] = 2;

	rst13_GLL[0] = -1.0;
	rst13_GLL[1] =  1.0;
	rst13_GLL[2] = -1.0/5.0*sqrt(5.0);
	rst13_GLL[3] =  1.0/5.0*sqrt(5.0);

	w13_GLL[0] = 1.0/6.0;
	w13_GLL[1] = 1.0/6.0;
	w13_GLL[2] = 5.0/6.0;
	w13_GLL[3] = 5.0/6.0;

	symms13_GLL[0] = 2;
	symms13_GLL[1] = 2;

	// GL
	cubature_TP(&rst,&w,&symms,&Nn,&Ns,1,P,dE,"GL"); // free

	pass = 0;
	if (array_norm_diff_d(N13,rst13_GL,rst,"Inf")       < EPS &&
	    array_norm_diff_d(N13,w13_GL,w,"Inf")           < EPS &&
	    array_norm_diff_ui(Ns13,symms13_GL,symms,"Inf") < EPS &&
	    Nn == N13 && Ns == Ns13)
			pass = 1;

	test_print2(pass,"cubature_TP (d1, P3, GL):");

	free(rst), free(w), free(symms);

	// GLL
	cubature_TP(&rst,&w,&symms,&Nn,&Ns,1,P,dE,"GLL"); // free

	pass = 0;
	if (array_norm_diff_d(N13,rst13_GLL,rst,"Inf")       < EPS &&
	    array_norm_diff_d(N13,w13_GLL,w,"Inf")           < EPS &&
	    array_norm_diff_ui(Ns13,symms13_GLL,symms,"Inf") < EPS &&
	    Nn == N13 && Ns == Ns13)
			pass = 1;

	test_print2(pass,"            (d1, P3, GLL):");

	free(rst), free(w), free(symms);


	// P = 4;
	P = 4;

	unsigned int N14 = P+1, Ns14 = 3;
	unsigned int symms14_GL[Ns14], symms14_GLL[Ns14];
	double rst14_GL[N14], w14_GL[N14], rst14_GLL[N14], w14_GLL[N14];

	rst14_GL[0] = -1.0/21.0*sqrt(245.0+14.0*sqrt(70.0));
	rst14_GL[1] =  1.0/21.0*sqrt(245.0+14.0*sqrt(70.0));
	rst14_GL[2] = -1.0/21.0*sqrt(245.0-14.0*sqrt(70.0));
	rst14_GL[3] =  1.0/21.0*sqrt(245.0-14.0*sqrt(70.0));
	rst14_GL[4] =  0.0;

	w14_GL[0] = 1.0/900.0*(322.0-13.0*sqrt(70.0));
	w14_GL[1] = 1.0/900.0*(322.0-13.0*sqrt(70.0));
	w14_GL[2] = 1.0/900.0*(322.0+13.0*sqrt(70.0));
	w14_GL[3] = 1.0/900.0*(322.0+13.0*sqrt(70.0));
	w14_GL[4] = 128.0/225.0;

	symms14_GL[0] = 2;
	symms14_GL[1] = 2;
	symms14_GL[2] = 1;

	rst14_GLL[0] = -1.0;
	rst14_GLL[1] =  1.0;
	rst14_GLL[2] = -1.0/7.0*sqrt(21.0);
	rst14_GLL[3] =  1.0/7.0*sqrt(21.0);
	rst14_GLL[4] =  0.0;

	w14_GLL[0] = 1.0/10.0;
	w14_GLL[1] = 1.0/10.0;
	w14_GLL[2] = 49.0/90.0;
	w14_GLL[3] = 49.0/90.0;
	w14_GLL[4] = 32.0/45.0;

	symms14_GLL[0] = 2;
	symms14_GLL[1] = 2;
	symms14_GLL[2] = 1;

	// GL
	cubature_TP(&rst,&w,&symms,&Nn,&Ns,1,P,dE,"GL"); // free

	pass = 0;
	if (array_norm_diff_d(N14,rst14_GL,rst,"Inf")       < EPS &&
	    array_norm_diff_d(N14,w14_GL,w,"Inf")           < EPS &&
	    array_norm_diff_ui(Ns14,symms14_GL,symms,"Inf") < EPS &&
	    Nn == N14 && Ns == Ns14)
			pass = 1;

	test_print2(pass,"            (d1, P4, GL):");

	free(rst), free(w), free(symms);

	// GLL
	cubature_TP(&rst,&w,&symms,&Nn,&Ns,1,P,dE,"GLL"); // free

	pass = 0;
	if (array_norm_diff_d(N14,rst14_GLL,rst,"Inf")       < EPS &&
	    array_norm_diff_d(N14,w14_GLL,w,"Inf")           < EPS &&
	    array_norm_diff_ui(Ns14,symms14_GLL,symms,"Inf") < EPS &&
	    Nn == N14 && Ns == Ns14)
			pass = 1;

	test_print2(pass,"            (d1, P4, GLL):");

	free(rst), free(w), free(symms);

	/*
	 *	Cubature_TP (dE = 2):
	 *
	 *		Input:
	 *
	 *			P, dE, NodeType.
	 *
	 *		Expected output:
	 *
	 *			P = 3, GL/GLL:
	 *				rst : Tensor-Product extension of rst13_*
	 *				w   : Tensor-Product extension of w13_*
	 *				Nn = 16
	 */

	dE = 2;

	// P = 3
	P = 3;

	unsigned int N23 = pow(P+1,dE);
	double rst23_GL[dE*N23], w23_GL[dE*N23], rst23_GLL[dE*N23], w23_GLL[dE*N23];

	row = 0;
	for (j = 0; j < N13; j++) {
	for (i = 0; i < N13; i++) {
		w23_GL[row]  = w13_GL[i]*w13_GL[j];
		w23_GLL[row] = w13_GLL[i]*w13_GLL[j];

		rst23_GL[0*N23+row] = rst13_GL[i];
		rst23_GL[1*N23+row] = rst13_GL[j];

		rst23_GLL[0*N23+row] = rst13_GLL[i];
		rst23_GLL[1*N23+row] = rst13_GLL[j];

		row++;
	}}

	// GL
	cubature_TP(&rst,&w,&symms,&Nn,&Ns,1,P,dE,"GL"); // free

	pass = 0;
	if (array_norm_diff_d(N23,rst23_GL,rst,"Inf") < EPS &&
	    array_norm_diff_d(N23,w23_GL,w,"Inf")     < EPS &&
	    Nn == N23)
			pass = 1;

	test_print2(pass,"            (d2, P3, GL):");

	free(rst), free(w), free(symms);

	// GLL
	cubature_TP(&rst,&w,&symms,&Nn,&Ns,1,P,dE,"GLL"); // free

	pass = 0;
	if (array_norm_diff_d(N23,rst23_GLL,rst,"Inf") < EPS &&
	    array_norm_diff_d(N23,w23_GLL,w,"Inf")     < EPS &&
	    Nn == N23)
			pass = 1;

	test_print2(pass,"            (d2, P3, GLL):");

	free(rst), free(w), free(symms);

	/*
	 *	Cubature_TP (dE = 3):
	 *
	 *		Input:
	 *
	 *			P, dE, NodeType.
	 *
	 *		Expected output:
	 *
	 *			P = 3, GL/GLL:
	 *				rst : Tensor-Product extension of rst13_*
	 *				w   : Tensor-Product extension of w13_*
	 *				Nn = 64
	 */

	dE = 3;

	// P = 3
	P = 3;

	unsigned int N33 = pow(P+1,dE);
	double rst33_GL[dE*N33], w33_GL[dE*N33], rst33_GLL[dE*N33], w33_GLL[dE*N33];

	row = 0;
	for (k = 0; k < N13; k++) {
	for (j = 0; j < N13; j++) {
	for (i = 0; i < N13; i++) {
		w33_GL[row]  = w13_GL[i]*w13_GL[j]*w13_GL[k];
		w33_GLL[row] = w13_GLL[i]*w13_GLL[j]*w13_GLL[k];

		rst33_GL[0*N33+row] = rst13_GL[i];
		rst33_GL[1*N33+row] = rst13_GL[j];
		rst33_GL[2*N33+row] = rst13_GL[k];

		rst33_GLL[0*N33+row] = rst13_GLL[i];
		rst33_GLL[1*N33+row] = rst13_GLL[j];
		rst33_GLL[2*N33+row] = rst13_GLL[k];

		row++;
	}}}

	// GL
	cubature_TP(&rst,&w,&symms,&Nn,&Ns,1,P,dE,"GL"); // free

	pass = 0;
	if (array_norm_diff_d(N33,rst33_GL,rst,"Inf") < EPS &&
	    array_norm_diff_d(N33,w33_GL,w,"Inf")     < EPS &&
	    Nn == N33)
			pass = 1;

	test_print2(pass,"            (d3, P3, GL):");

	free(rst), free(w), free(symms);

	// GLL
	cubature_TP(&rst,&w,&symms,&Nn,&Ns,1,P,dE,"GLL"); // free

	pass = 0;
	if (array_norm_diff_d(N33,rst33_GLL,rst,"Inf") < EPS &&
	    array_norm_diff_d(N33,w33_GLL,w,"Inf")     < EPS &&
	    Nn == N33)
			pass = 1;

	test_print2(pass,"            (d3, P3, GLL):");

	free(rst), free(w), free(symms);
}

void test_unit_cubature_SI(void)
{
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
	 *			P = 3, WSH:
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
			pass = 1;

	test_print2(pass,"cubature_SI (d2, P3, AO):");

	free(rst);
	free(symms);

	// WSH (P = 3)
	P = 3;
	unsigned int Nn23_WSH = 10, Ns23_WSH = 4;
	unsigned int symms23_WSH[4] = { 3, 3, 3, 1 };
	double rst23_WSH[20] = { -0.338677036009830,  -0.455664103498571,
	                          0.563955207227339,  -0.065470865113645,
	                         -0.225278171217509,   0.521134968612215,
	                         -0.563955207227339,  -0.065470865113645,
	                          0.338677036009830,  -0.455664103498571,
	                          0.225278171217509,   0.521134968612215,
	                         -0.833307841990621,  -0.481110506891111,
	                          0.833307841990621,  -0.481110506891111,
	                         -0.000000000000000,   0.962221013782223,
	                         -0.000000000000000,  -0.000000000000001};

	double w23_WSH[10] = { 0.194160145154569,
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
	mkl_dimatcopy('R','T',Nn23_WSH,d,1.0,rst23_WSH,d,Nn23_WSH);

	cubature_TRI(&rst,&w,&symms,&Nn,&Ns,1,P,d,"WSH");

	pass = 0;
	if (array_norm_diff_d(Nn23_WSH*d,rst23_WSH,rst,"Inf")    < EPS &&
	    array_norm_diff_d(Nn23_WSH,w23_WSH,w,"Inf")          < EPS &&
	    array_norm_diff_ui(Ns23_WSH,symms23_WSH,symms,"Inf") < EPS &&
	    Nn23_WSH == Nn && Ns23_WSH == Ns)
			pass = 1;

	test_print2(pass,"            (d2, P3, WSH):");

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
			pass = 1;

	test_print2(pass,"            (d2, P6, WV):");

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
	 *			P = 2, WSH:
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
			pass = 1;

	test_print2(pass,"            (d3, P2, AO):");

	free(rst);
	free(symms);

	// WSH (P = 2)
	P = 2;
	unsigned int Nn33_WSH = 10, Ns33_WSH = 4;
	unsigned int symms33_WSH[4] = { 3, 3, 3, 1 };
	double rst33_WSH[30] = { 0.000000000000000, -0.360830856003254, -0.255145945141247,
	                         0.312488687768103,  0.180415428001627, -0.255145945141247,
	                        -0.312488687768103,  0.180415428001627, -0.255145945141247,
	                        -0.312488687768103, -0.180415428001627,  0.255145945141247,
	                         0.312488687768103, -0.180415428001627,  0.255145945141247,
	                         0.000000000000000,  0.360830856003254,  0.255145945141247,
	                        -0.704660393095107, -0.406835867640727, -0.287676400838671,
	                         0.704660393095107, -0.406835867640727, -0.287676400838671,
	                         0.000000000000000,  0.813671735281455, -0.287676400838671,
	                         0.000000000000000,  0.000000000000000,  0.863029202516013};

	double w33_WSH[10] = { 0.127195540124294,
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
	mkl_dimatcopy('R','T',Nn33_WSH,d,1.0,rst33_WSH,d,Nn33_WSH);

	cubature_TET(&rst,&w,&symms,&Nn,&Ns,1,P,d,"WSH");

	pass = 0;
	if (array_norm_diff_d(Nn33_WSH*d,rst33_WSH,rst,"Inf")    < EPS &&
	    array_norm_diff_d(Nn33_WSH,w33_WSH,w,"Inf")          < EPS &&
	    array_norm_diff_ui(Ns33_WSH,symms33_WSH,symms,"Inf") < EPS &&
	    Nn33_WSH == Nn && Ns33_WSH == Ns)
			pass = 1;

	test_print2(pass,"            (d3, P2, WSH):");

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
			pass = 1;

	test_print2(pass,"            (d3, P4, WV):");

	free(rst), free(w);
	free(symms);
}

void test_unit_cubature_PYR(void)
{
	/*
	 *	Purpose:
	 *		Test correctness of implementation of cubature_PYR.
	 *
	 *	Comments:
	 *		Other than the grouping based on 1/4 symmetries, the ordering of the pyramidal nodes is arbitrary.
	 *
	 *	Notation:
	 *
	 *	References:
	 */

	unsigned int pass;

	/*
	 *	Input:
	 *
	 *		P, NodeType.
	 *
	 *	Expected output:
	 *
	 *		P = 2, GL:
	 *			rst = [ See below ]
	 *			w   = N/A
	 *			symms = [ 4 4 4 1 1 ]
	 *			Ns  = 5
	 *			Nn  = 14
	 *
	 *		P = 3, GLL:
	 *			rst = [ See below ]
	 *			w   = N/A
	 *			symms = [ 8 4 4 4 4 1 1 ]
	 *			Ns  = 8
	 *			Nn  = 30
	 *
	 *		P = 4, WV:
	 *			rst = [ See below ]
	 *			w   = [ See below ]
	 *			symms = [ 4 4 1 1 ]
	 *			Ns  = 4
	 *			Nn  = 10
	 */

	unsigned int Nn, Ns, P, d;
	unsigned int *symms;
	double *rst, *w;

	d = 3;

	// GL (P = 2)
	P = 2;
	unsigned int Nn2_GL = 14, Ns2_GL = 5;
	unsigned int symms2_GL[5] = { 4, 4, 4, 1, 1 };
	double rst2_GL[42] = { -0.687298334620742, -0.687298334620742, -0.123458488793238,
	                        0.687298334620742, -0.687298334620742, -0.123458488793238,
	                        0.687298334620742,  0.687298334620742, -0.123458488793238,
	                       -0.687298334620742,  0.687298334620742, -0.123458488793238,
	                        0.000000000000000, -0.687298334620742, -0.123458488793238,
	                        0.687298334620742, -0.000000000000000, -0.123458488793238,
	                        0.000000000000000,  0.687298334620742, -0.123458488793238,
	                       -0.687298334620742,  0.000000000000000, -0.123458488793238,
	                       -0.288675134594813, -0.288675134594813,  0.424264068711929,
	                        0.288675134594813, -0.288675134594813,  0.424264068711929,
	                        0.288675134594813,  0.288675134594813,  0.424264068711929,
	                       -0.288675134594813,  0.288675134594813,  0.424264068711929,
	                        0.000000000000000,  0.000000000000000, -0.123458488793238,
	                        0.000000000000000,  0.000000000000000,  0.971986626217095};

	// Convert to column-major ordering
	mkl_dimatcopy('R','T',Nn2_GL,d,1.0,rst2_GL,d,Nn2_GL);

	cubature_PYR(&rst,&w,&symms,&Nn,&Ns,0,P,d,"GL");

	pass = 0;
	if (array_norm_diff_d(Nn2_GL*d,rst2_GL,rst,"Inf")    < EPS &&
	    array_norm_diff_ui(Ns2_GL,symms2_GL,symms,"Inf") < EPS &&
	    Nn2_GL == Nn && Ns2_GL == Ns)
			pass = 1;

	test_print2(pass,"cubature_PYR (P2, GL) :");

	free(rst);
	free(symms);

	// GLL (P = 3)
	P = 3;
	unsigned int Nn3_GLL = 30, Ns3_GLL = 9;
	unsigned int symms3_GLL[9] = { 4, 4, 4, 4, 4, 4, 4, 1, 1 };
	double rst3_GLL[90] = { -0.447213595499958, -1.000000000000000, -0.282842712474619,
	                         1.000000000000000, -0.447213595499958, -0.282842712474619,
	                         0.447213595499958,  1.000000000000000, -0.282842712474619,
	                        -1.000000000000000,  0.447213595499958, -0.282842712474619,
	                         0.447213595499958, -1.000000000000000, -0.282842712474619,
	                         1.000000000000000,  0.447213595499958, -0.282842712474619,
	                        -0.447213595499958,  1.000000000000000, -0.282842712474619,
	                        -1.000000000000000, -0.447213595499958, -0.282842712474619,
	                        -1.000000000000000, -1.000000000000000, -0.282842712474619,
	                         1.000000000000000, -1.000000000000000, -0.282842712474619,
	                         1.000000000000000,  1.000000000000000, -0.282842712474619,
	                        -1.000000000000000,  1.000000000000000, -0.282842712474619,
	                        -0.723606797749979, -0.723606797749979,  0.108036302695091,
	                         0.723606797749979, -0.723606797749979,  0.108036302695091,
	                         0.723606797749979,  0.723606797749979,  0.108036302695091,
	                        -0.723606797749979,  0.723606797749979,  0.108036302695091,
	                        -0.447213595499958, -0.447213595499958, -0.282842712474619,
	                         0.447213595499958, -0.447213595499958, -0.282842712474619,
	                         0.447213595499958,  0.447213595499958, -0.282842712474619,
	                        -0.447213595499958,  0.447213595499958, -0.282842712474619,
	                         0.000000000000000, -0.723606797749979,  0.108036302695091,
	                         0.723606797749979,  0.000000000000000,  0.108036302695091,
	                         0.000000000000000,  0.723606797749979,  0.108036302695091,
	                        -0.723606797749979,  0.000000000000000,  0.108036302695091,
	                        -0.276393202250021, -0.276393202250021,  0.740491834728766,
	                         0.276393202250021, -0.276393202250021,  0.740491834728766,
	                         0.276393202250021,  0.276393202250021,  0.740491834728766,
	                        -0.276393202250021,  0.276393202250021,  0.740491834728766,
	                        -0.000000000000000, -0.000000000000000,  0.108036302695091,
	                         0.000000000000000,  0.000000000000000,  1.131370849898476};

	// Convert to column-major ordering
	mkl_dimatcopy('R','T',Nn3_GLL,d,1.0,rst3_GLL,d,Nn3_GLL);

	cubature_PYR(&rst,&w,&symms,&Nn,&Ns,0,P,d,"GLL");

	pass = 0;
	if (array_norm_diff_d(Nn3_GLL*d,rst3_GLL,rst,"Inf")    < EPS &&
	    array_norm_diff_ui(Ns3_GLL,symms3_GLL,symms,"Inf") < EPS &&
	    Nn3_GLL == Nn && Ns3_GLL == Ns)
			pass = 1;

	test_print2(pass,"             (P3, GLL):");

	free(rst);
	free(symms);

	// WV (P = 4)
	P = 4;
	unsigned int Nn4_WV = 10, Ns4_WV = 4;
	unsigned int symms4_WV[4] = { 4, 4, 1, 1 };
	double rst4_WV[30] = { -0.657966997121690, -0.657966997121690, -0.227337257085045,
	                        0.657966997121690, -0.657966997121690, -0.227337257085045,
	                        0.657966997121690,  0.657966997121690, -0.227337257085045,
	                       -0.657966997121690,  0.657966997121690, -0.227337257085045,
	                       -0.000000000000000, -0.650581556398233,  0.173077324153007,
	                        0.650581556398233,  0.000000000000000,  0.173077324153007,
	                        0.000000000000000,  0.650581556398233,  0.173077324153007,
	                       -0.650581556398233,  0.000000000000000,  0.173077324153007,
	                        0.000000000000000,  0.000000000000000, -0.105872336234184,
	                        0.000000000000000,  0.000000000000000,  0.674909082451912};

	double w4_WV[10] = { 0.119772745994674,
	                     0.119772745994674,
	                     0.119772745994674,
	                     0.119772745994674,
	                     0.200487565609086,
	                     0.200487565609086,
	                     0.200487565609086,
	                     0.200487565609086,
	                     0.390103085029384,
	                     0.214473751719704};

	// Convert to column-major ordering
	mkl_dimatcopy('R','T',Nn4_WV,d,1.0,rst4_WV,d,Nn4_WV);

	cubature_PYR(&rst,&w,&symms,&Nn,&Ns,1,P,d,"WV");

	pass = 0;
	if (array_norm_diff_d(Nn4_WV*d,rst4_WV,rst,"Inf")    < EPS*10 &&
	    array_norm_diff_d(Nn4_WV,w4_WV,w,"Inf")          < EPS    &&
	    array_norm_diff_ui(Ns4_WV,symms4_WV,symms,"Inf") < EPS    &&
	    Nn4_WV == Nn && Ns4_WV == Ns)
			pass = 1;

	test_print2(pass,"             (P4, WV) :");

	free(rst);
	free(w);
	free(symms);
}
