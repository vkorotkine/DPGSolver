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
	 *			P, dE, NodeType.
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
	 */

	unsigned int i, j, k, row;
	unsigned int dE, Nn, P;
	double *rst, *w;

	dE = 1;

	// P = 3
	P = 3;

	unsigned int N13 = P+1;
	double rst13_GL[N13], w13_GL[N13], rst13_GLL[N13], w13_GLL[N13];

	rst13_GL[0] = -1.0/35.0*sqrt(525.0+70.0*sqrt(30.0));
	rst13_GL[1] =  1.0/35.0*sqrt(525.0+70.0*sqrt(30.0));
	rst13_GL[2] = -1.0/35.0*sqrt(525.0-70.0*sqrt(30.0));
	rst13_GL[3] =  1.0/35.0*sqrt(525.0-70.0*sqrt(30.0));

	w13_GL[0] = 1.0/36.0*(18.0-sqrt(30.0));
	w13_GL[1] = 1.0/36.0*(18.0-sqrt(30.0));
	w13_GL[2] = 1.0/36.0*(18.0+sqrt(30.0));
	w13_GL[3] = 1.0/36.0*(18.0+sqrt(30.0));

	rst13_GLL[0] = -1.0;
	rst13_GLL[1] =  1.0;
	rst13_GLL[2] = -1.0/5.0*sqrt(5.0);
	rst13_GLL[3] =  1.0/5.0*sqrt(5.0);

	w13_GLL[0] = 1.0/6.0;
	w13_GLL[1] = 1.0/6.0;
	w13_GLL[2] = 5.0/6.0;
	w13_GLL[3] = 5.0/6.0;

	// GL
	cubature_TP(&rst,&w,&Nn,1,P,dE,"GL"); // free

	pass = 0;
	if (array_norm_diff_d(N13,rst13_GL,rst,"Inf") < EPS &&
		array_norm_diff_d(N13,w13_GL,w,"Inf")     < EPS &&
		Nn == N13)
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("cubature_TP (d1, P3, GL):                        ");
	test_print(pass);

	free(rst), free(w);

	// GLL
	cubature_TP(&rst,&w,&Nn,1,P,dE,"GLL"); // free

	pass = 0;
	if (array_norm_diff_d(N13,rst13_GLL,rst,"Inf") < EPS &&
		array_norm_diff_d(N13,w13_GLL,w,"Inf")     < EPS &&
		Nn == N13)
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("            (d1, P3, GLL):                       ");
	test_print(pass);

	free(rst), free(w);


	// P = 4;
	P = 4;

	unsigned int N14 = P+1;
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

	// GL
	cubature_TP(&rst,&w,&Nn,1,P,dE,"GL"); // free

	pass = 0;
	if (array_norm_diff_d(N14,rst14_GL,rst,"Inf") < EPS &&
		array_norm_diff_d(N14,w14_GL,w,"Inf")     < EPS &&
		Nn == N14)
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("            (d1, P4, GL):                        ");
	test_print(pass);

	free(rst), free(w);

	// GLL
	cubature_TP(&rst,&w,&Nn,1,P,dE,"GLL"); // free

	pass = 0;
	if (array_norm_diff_d(N14,rst14_GLL,rst,"Inf") < EPS &&
		array_norm_diff_d(N14,w14_GLL,w,"Inf")     < EPS &&
		Nn == N14)
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("            (d1, P4, GLL):                       ");
	test_print(pass);

	free(rst), free(w);

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
	 *
	 *				Nodes:
	 *				12 - 13 - 14 - 15
	 *				 8 -  9 - 10 - 11
	 *				 4 -  5 -  6 -  7
	 *				 0 -  1 -  2 -  3
	 *
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
	cubature_TP(&rst,&w,&Nn,1,P,dE,"GL"); // free

	pass = 0;
	if (array_norm_diff_d(N23,rst23_GL,rst,"Inf") < EPS &&
		array_norm_diff_d(N23,w23_GL,w,"Inf")     < EPS &&
		Nn == N23)
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("            (d2, P3, GL):                        ");
	test_print(pass);

	free(rst), free(w);

	// GLL
	cubature_TP(&rst,&w,&Nn,1,P,dE,"GLL"); // free

	pass = 0;
	if (array_norm_diff_d(N23,rst23_GLL,rst,"Inf") < EPS &&
		array_norm_diff_d(N23,w23_GLL,w,"Inf")     < EPS &&
		Nn == N23)
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("            (d2, P3, GLL):                       ");
	test_print(pass);

	free(rst), free(w);

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
	cubature_TP(&rst,&w,&Nn,1,P,dE,"GL"); // free

	pass = 0;
	if (array_norm_diff_d(N33,rst33_GL,rst,"Inf") < EPS &&
		array_norm_diff_d(N33,w33_GL,w,"Inf")     < EPS &&
		Nn == N33)
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("            (d3, P3, GL):                        ");
	test_print(pass);

	free(rst), free(w);

	// GLL
	cubature_TP(&rst,&w,&Nn,1,P,dE,"GLL"); // free

	pass = 0;
	if (array_norm_diff_d(N33,rst33_GLL,rst,"Inf") < EPS &&
		array_norm_diff_d(N33,w33_GLL,w,"Inf")     < EPS &&
		Nn == N33)
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("            (d3, P3, GLL):                       ");
	test_print(pass);

	free(rst), free(w);
}
