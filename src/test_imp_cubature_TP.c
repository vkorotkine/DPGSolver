#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "test.h"
#include "functions.h"
#include "parameters.h"

/*
 *	Purpose:
 *		Test correctness of implementation of cubature_TP.
 *
 *	Comments:
 *		'Con' arrays may need to be modified based on plotting software used. (ToBeDeleted)
 *
 *	Notation:
 *
 *	References:
 *		GL  : http://mathworld.wolfram.com/Legendre-GaussQuadrature.html
 *		GLL : http://mathworld.wolfram.com/LobattoQuadrature.html
 */

void test_imp_cubature_TP(void)
{
	unsigned int pass;

	/*
	 *	cubature_TP (dE = 1):
	 *
	 *		Input:
	 *
	 *			ToReturn, P, dE, NodeType.
	 *
	 *		Expected output:
	 *
	 *			P = 3/4, GL:
	 *				rst = [ See below ]
	 *				w   = [ See below ]
	 *				Nn  = 4/5
	 *
	 *			P = 3/4, GLL:
	 *				rst = [ See below ]
	 *				w   = [ See below ]
	 *				Nn  = 4/5
	 *
	 *			P = 3/4, ES:
	 *				rst = [ See below ]
	 *				Nn  = 4/5
	 *
	 *			P = 3:
	 *				Nodes:
	 *				0 - 1 - 2 - 3
	 *
	 *				Con:
	 *				[ 0 1; 1 2; 2 3 ]
	 */

	unsigned int i, j, k, row;
	unsigned int dE, ToReturn[4], Nn, P;
	unsigned int *Con;
	double *rst, *W;

	dE = 1;

	// P = 3
	P = 3;

	unsigned int N13 = P+1;
	unsigned int Con13[6] = { 0, 1, 1, 2, 2, 3 };
	double rst13_GL[N13], W13_GL[N13], rst13_GLL[N13], W13_GLL[N13], rst13_ES[N13];

	rst13_GL[0] = -1.0/35.0*sqrt(525.0-70.0*sqrt(30.0));
	rst13_GL[1] =  1.0/35.0*sqrt(525.0-70.0*sqrt(30.0));
	rst13_GL[2] = -1.0/35.0*sqrt(525.0+70.0*sqrt(30.0));
	rst13_GL[3] =  1.0/35.0*sqrt(525.0+70.0*sqrt(30.0));

	W13_GL[0] = 1.0/36.0*(18.0+sqrt(30.0));
	W13_GL[1] = 1.0/36.0*(18.0+sqrt(30.0));
	W13_GL[2] = 1.0/36.0*(18.0-sqrt(30.0));
	W13_GL[3] = 1.0/36.0*(18.0-sqrt(30.0));

	rst13_GLL[0] = -1.0/5.0*sqrt(5.0);
	rst13_GLL[1] =  1.0/5.0*sqrt(5.0);
	rst13_GLL[2] = -1.0;
	rst13_GLL[3] =  1.0;

	W13_GLL[0] = 5.0/6.0;
	W13_GLL[1] = 5.0/6.0;
	W13_GLL[2] = 1.0/6.0;
	W13_GLL[3] = 1.0/6.0;

	rst13_ES[0] = -1.0;
	rst13_ES[1] = -1.0/3.0;
	rst13_ES[2] =  1.0/3.0;
	rst13_ES[3] =  1.0;


	// GL
	ToReturn[0] = 1; ToReturn[1] = 1; ToReturn[2] = 0; ToReturn[3] = 1;
	cubature_TP(&rst,&W,&Con,&Nn,ToReturn,P,dE,"GL"); // free

	pass = 0;
	if (array_norm_diff_d(N13,rst13_GL,rst,"Inf") < EPS &&
		array_norm_diff_d(N13,W13_GL,W,"Inf")     < EPS &&
		Nn == N13)
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("cubature_TP (d1, P3, GL):                        ");
	test_print(pass);

	free(rst), free(W);

	// GLL
	ToReturn[0] = 1; ToReturn[1] = 1; ToReturn[2] = 0; ToReturn[3] = 1;
	cubature_TP(&rst,&W,&Con,&Nn,ToReturn,P,dE,"GLL"); // free

	pass = 0;
	if (array_norm_diff_d(N13,rst13_GLL,rst,"Inf") < EPS &&
		array_norm_diff_d(N13,W13_GLL,W,"Inf")     < EPS &&
		Nn == N13)
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("cubature_TP (d1, P3, GLL):                       ");
	test_print(pass);

	free(rst), free(W);

	// ES
	ToReturn[0] = 1; ToReturn[1] = 0; ToReturn[2] = 1; ToReturn[3] = 1;
	cubature_TP(&rst,&W,&Con,&Nn,ToReturn,P,dE,"ES"); // free

	pass = 0;
	if (array_norm_diff_d(N13,rst13_ES,rst,"Inf") < EPS &&
		array_norm_diff_ui(2*Nn,Con13,Con,"Inf")  < EPS &&
		Nn == P)
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("cubature_TP (d1, P3, ES):                        ");
	test_print(pass);

	free(rst), free(Con);


	// P = 4;
	P = 4;

	unsigned int N14 = P+1;
	unsigned int Con14[8] = { 0, 1, 1, 2, 2, 3, 3, 4 };
	double rst14_GL[N14], W14_GL[N14], rst14_GLL[N14], W14_GLL[N14], rst14_ES[N14];

	rst14_GL[0] =  0.0;
	rst14_GL[1] = -1.0/21.0*sqrt(245.0-14.0*sqrt(70.0));
	rst14_GL[2] =  1.0/21.0*sqrt(245.0-14.0*sqrt(70.0));
	rst14_GL[3] = -1.0/21.0*sqrt(245.0+14.0*sqrt(70.0));
	rst14_GL[4] =  1.0/21.0*sqrt(245.0+14.0*sqrt(70.0));

	W14_GL[0] = 128.0/225.0;
	W14_GL[1] = 1.0/900.0*(322.0+13.0*sqrt(70.0));
	W14_GL[2] = 1.0/900.0*(322.0+13.0*sqrt(70.0));
	W14_GL[3] = 1.0/900.0*(322.0-13.0*sqrt(70.0));
	W14_GL[4] = 1.0/900.0*(322.0-13.0*sqrt(70.0));

	rst14_GLL[0] =  0.0;
	rst14_GLL[1] = -1.0/7.0*sqrt(21.0);
	rst14_GLL[2] =  1.0/7.0*sqrt(21.0);
	rst14_GLL[3] = -1.0;
	rst14_GLL[4] =  1.0;

	W14_GLL[0] = 32.0/45.0;
	W14_GLL[1] = 49.0/90.0;
	W14_GLL[2] = 49.0/90.0;
	W14_GLL[3] = 1.0/10.0;
	W14_GLL[4] = 1.0/10.0;

	rst14_ES[0] = -1.0;
	rst14_ES[1] = -0.5;
	rst14_ES[2] =  0.0;
	rst14_ES[3] =  0.5;
	rst14_ES[4] =  1.0;


	// GL
	ToReturn[0] = 1; ToReturn[1] = 1; ToReturn[2] = 0; ToReturn[3] = 1;
	cubature_TP(&rst,&W,&Con,&Nn,ToReturn,P,dE,"GL"); // free

	pass = 0;
	if (array_norm_diff_d(N14,rst14_GL,rst,"Inf") < EPS &&
		array_norm_diff_d(N14,W14_GL,W,"Inf")     < EPS &&
		Nn == N14)
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("cubature_TP (d1, P4, GL):                        ");
	test_print(pass);

	free(rst), free(W);

	// GLL
	ToReturn[0] = 1; ToReturn[1] = 1; ToReturn[2] = 0; ToReturn[3] = 1;
	cubature_TP(&rst,&W,&Con,&Nn,ToReturn,P,dE,"GLL"); // free

	pass = 0;
	if (array_norm_diff_d(N14,rst14_GLL,rst,"Inf") < EPS &&
		array_norm_diff_d(N14,W14_GLL,W,"Inf")     < EPS &&
		Nn == N14)
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("cubature_TP (d1, P4, GLL):                       ");
	test_print(pass);

	free(rst), free(W);


	// ES
	ToReturn[0] = 1; ToReturn[1] = 0; ToReturn[2] = 1; ToReturn[3] = 1;
	cubature_TP(&rst,&W,&Con,&Nn,ToReturn,P,dE,"ES"); // free

	pass = 0;
	if (array_norm_diff_d(N14,rst14_ES,rst,"Inf") < EPS &&
		array_norm_diff_ui(2*Nn,Con14,Con,"Inf")   < EPS &&
		Nn == P)
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("cubature_TP (d1, P4, ES):                        ");
	test_print(pass);

	free(rst), free(Con);

	/*
	 *	Cubature_TP (dE = 2):
	 *
	 *		Input:
	 *
	 *			ToReturn, P, dE, NodeType.
	 *
	 *		Expected output:
	 *
	 *			P = 3, GL/GLL/ES:
	 *				rst : Tensor-Product extension of rst13_*
	 *				w   : Tensor-Product extension of w13_*
	 *				Nn = 16
	 *
	 *				Nodes:
	 *				12 - 13 - 14 - 15
	 *				 8 -  9 - 10 - 11
	 *				 4 -  5 -  6 -  7
	 *				 0 -  1 -  2 -  3
	 *
	 *				Con:
	 *				[  0  1  5  4;  4  5  9  7;  8  9 13 12;
	 *				   1  2  6  5;  5  6 10  9;  9 10 14 13;
	 *				   2  3  7  6;  6  7 11 10; 10 11 15 14 ]
	 */

	dE = 2;

	// P = 3
	P = 3;

	unsigned int N23 = pow(P+1,dE);
	unsigned int Con23[9*4] = { 0,  1,  5,  4,   4,  5,  9,  8,   8,  9,  13, 12,
	                            1,  2,  6,  5,   5,  6,  10, 9,   9,  10, 14, 13,
	                            2,  3,  7,  6,   6,  7,  11, 10,  10, 11, 15, 14 };
	double rst23_GL[dE*N23], W23_GL[dE*N23], rst23_GLL[dE*N23], W23_GLL[dE*N23], rst23_ES[dE*N23];

	row = 0;
	for (j = 0; j < N13; j++) {
	for (i = 0; i < N13; i++) {
		W23_GL[row]  = W13_GL[i]*W13_GL[j];
		W23_GLL[row] = W13_GLL[i]*W13_GLL[j];

		rst23_GL[0*N23+row] = rst13_GL[i];
		rst23_GL[1*N23+row] = rst13_GL[j];

		rst23_GLL[0*N23+row] = rst13_GLL[i];
		rst23_GLL[1*N23+row] = rst13_GLL[j];

		rst23_ES[0*N23+row] = rst13_ES[i];
		rst23_ES[1*N23+row] = rst13_ES[j];

		row++;
	}}

	// GL
	ToReturn[0] = 1; ToReturn[1] = 1; ToReturn[2] = 0; ToReturn[3] = 1;
	cubature_TP(&rst,&W,&Con,&Nn,ToReturn,P,dE,"GL"); // free

	pass = 0;
	if (array_norm_diff_d(N23,rst23_GL,rst,"Inf") < EPS &&
		array_norm_diff_d(N23,W23_GL,W,"Inf")     < EPS &&
		Nn == N23)
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("cubature_TP (d2, P3, GL):                        ");
	test_print(pass);

	free(rst), free(W);

	// GLL
	ToReturn[0] = 1; ToReturn[1] = 1; ToReturn[2] = 0; ToReturn[3] = 1;
	cubature_TP(&rst,&W,&Con,&Nn,ToReturn,P,dE,"GLL"); // free

	pass = 0;
	if (array_norm_diff_d(N23,rst23_GLL,rst,"Inf") < EPS &&
		array_norm_diff_d(N23,W23_GLL,W,"Inf")     < EPS &&
		Nn == N23)
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("cubature_TP (d2, P3, GLL):                       ");
	test_print(pass);

	free(rst), free(W);

	// ES
	ToReturn[0] = 1; ToReturn[1] = 0; ToReturn[2] = 1; ToReturn[3] = 1;
	cubature_TP(&rst,&W,&Con,&Nn,ToReturn,P,dE,"ES"); // free

	pass = 0;
	if (array_norm_diff_d(N23,rst23_ES,rst,"Inf") < EPS &&
		array_norm_diff_ui(4*Nn,Con23,Con,"Inf")  < EPS &&
		Nn == (unsigned int)pow(P,dE))
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("cubature_TP (d2, P3, ES):                        ");
	test_print(pass);

	free(rst), free(Con);

	/*
	 *	Cubature_TP (dE = 3):
	 *
	 *		Input:
	 *
	 *			ToReturn, P, dE, NodeType.
	 *
	 *		Expected output:
	 *
	 *			P = 3, GL/GLL/ES:
	 *				rst : Tensor-Product extension of rst13_*
	 *				w   : Tensor-Product extension of w13_*
	 *				Con = [ Extrapolate from 2D, See below ]
	 *				Nn = 64
	 */

	dE = 3;

	// P = 3
	P = 3;

	unsigned int N33 = pow(P+1,dE);
	unsigned int Con33[27*8] = {  0,  1,  5,  4, 16, 17, 21, 20,   16, 17, 21, 20, 32, 33, 37, 36,
	                             32, 33, 37, 36, 48, 49, 53, 52,    4,  5,  9,  8, 20, 21, 25, 24,
	                             20, 21, 25, 24, 36, 37, 41, 40,   36, 37, 41, 40, 52, 53, 57, 56,
	                              8,  9, 13, 12, 24, 25, 29, 28,   24, 25, 29, 28, 40, 41, 45, 44,
	                             40, 41, 45, 44, 56, 57, 61, 60,    1,  2,  6,  5, 17, 18, 22, 21,
	                             17, 18, 22, 21, 33, 34, 38, 37,   33, 34, 38, 37, 49, 50, 54, 53,
	                              5,  6, 10,  9, 21, 22, 26, 25,   21, 22, 26, 25, 37, 38, 42, 41,
	                             37, 38, 42, 41, 53, 54, 58, 57,    9, 10, 14, 13, 25, 26, 30, 29,
	                             25, 26, 30, 29, 41, 42, 46, 45,   41, 42, 46, 45, 57, 58, 62, 61,
	                              2,  3,  7,  6, 18, 19, 23, 22,   18, 19, 23, 22, 34, 35, 39, 38,
	                             34, 35, 39, 38, 50, 51, 55, 54,    6,  7, 11, 10, 22, 23, 27, 26,
	                             22, 23, 27, 26, 38, 39, 43, 42,   38, 39, 43, 42, 54, 55, 59, 58,
	                             10, 11, 15, 14, 26, 27, 31, 30,   26, 27, 31, 30, 42, 43, 47, 46,
	                             42, 43, 47, 46, 58, 59, 63, 62 };
	double rst33_GL[dE*N33], W33_GL[dE*N33], rst33_GLL[dE*N33], W33_GLL[dE*N33], rst33_ES[dE*N33];

	row = 0;
	for (k = 0; k < N13; k++) {
	for (j = 0; j < N13; j++) {
	for (i = 0; i < N13; i++) {
		W33_GL[row]  = W13_GL[i]*W13_GL[j]*W13_GL[k];
		W33_GLL[row] = W13_GLL[i]*W13_GLL[j]*W13_GLL[k];

		rst33_GL[0*N33+row] = rst13_GL[i];
		rst33_GL[1*N33+row] = rst13_GL[j];
		rst33_GL[2*N33+row] = rst13_GL[k];

		rst33_GLL[0*N33+row] = rst13_GLL[i];
		rst33_GLL[1*N33+row] = rst13_GLL[j];
		rst33_GLL[2*N33+row] = rst13_GLL[k];

		rst33_ES[0*N33+row] = rst13_ES[i];
		rst33_ES[1*N33+row] = rst13_ES[j];
		rst33_ES[2*N33+row] = rst13_ES[k];

		row++;
	}}}

	// GL
	ToReturn[0] = 1; ToReturn[1] = 1; ToReturn[2] = 0; ToReturn[3] = 1;
	cubature_TP(&rst,&W,&Con,&Nn,ToReturn,P,dE,"GL"); // free

	pass = 0;
	if (array_norm_diff_d(N33,rst33_GL,rst,"Inf") < EPS &&
		array_norm_diff_d(N33,W33_GL,W,"Inf")     < EPS &&
		Nn == N33)
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("cubature_TP (d3, P3, GL):                        ");
	test_print(pass);

	free(rst), free(W);

	// GLL
	ToReturn[0] = 1; ToReturn[1] = 1; ToReturn[2] = 0; ToReturn[3] = 1;
	cubature_TP(&rst,&W,&Con,&Nn,ToReturn,P,dE,"GLL"); // free

	pass = 0;
	if (array_norm_diff_d(N33,rst33_GLL,rst,"Inf") < EPS &&
		array_norm_diff_d(N33,W33_GLL,W,"Inf")     < EPS &&
		Nn == N33)
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("cubature_TP (d3, P3, GLL):                       ");
	test_print(pass);

	free(rst), free(W);


	// ES
	ToReturn[0] = 1; ToReturn[1] = 0; ToReturn[2] = 1; ToReturn[3] = 1;
	cubature_TP(&rst,&W,&Con,&Nn,ToReturn,P,dE,"ES"); // free

	pass = 0;
	if (array_norm_diff_d(N33,rst33_ES,rst,"Inf") < EPS &&
		array_norm_diff_ui(8*Nn,Con33,Con,"Inf")  < EPS &&
		Nn == (unsigned int)pow(P,dE))
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("cubature_TP (d3, P3, ES):                        ");
	test_print(pass);

	free(rst), free(Con);

}
