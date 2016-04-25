#include <stdlib.h>
#include <stdio.h>

#include "test.h"
#include "functions.h"
#include "parameters.h"

/*
 *	Purpose:
 *		Test correctness of implementation of cubature_TRI.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

void test_imp_cubature_TRI(void)
{
	unsigned int pass;

	/*
	 *	Input:
	 *
	 *		P, NodeType.
	 *
	 *	Expected output:
	 *
	 *		P = 3, AO:
	 *			rst = [ See below ]
	 *			w   = N/A
	 *			symms = [ 3 3 3 1 ]
	 *			Ns  = 4
	 *			Nn  = 10
	 *
	 *		P = 3, WS:
	 *			rst = [ See below ]
	 *			w   = [ See below ]
	 *			symms = [ 3 3 3 1 ]
	 *			Ns  = 4
	 *			Nn  = 10
	 *
	 *		P = 6, WV:
	 *			rst = [ See below ]
	 *			w   = [ See below ]
	 *			symms = [ 3 3 3 3 ]
	 *			Ns  = 4
	 *			Nn  = 12
	 */

	unsigned int Nn, Ns, P, d;
	unsigned int *symms;
	double *rst, *w;

	d = 2;

	// AO (P = 3)
	P = 3;
	unsigned int Nn3_AO = 10, Ns3_AO = 4;
	unsigned int symms3_AO[4] = { 3, 3, 3, 1 };
	double rst3_AO[20] = { -0.447213595499958, -0.577350269189626,
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
	mkl_dimatcopy('R','T',Nn3_AO,d,1.0,rst3_AO,d,Nn3_AO);

	cubature_TRI(&rst,&w,&symms,&Nn,&Ns,0,P,d,"AO");

	pass = 0;
	if (array_norm_diff_d(Nn3_AO*d,rst3_AO,rst,"Inf")    < EPS &&
	    array_norm_diff_ui(Ns3_AO,symms3_AO,symms,"Inf") < EPS &&
	    Nn3_AO == Nn && Ns3_AO == Ns)
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("cubature_TRI (P3, AO):                           ");
	test_print(pass);

	free(rst);
	free(symms);

	// WS (P = 3)
	P = 3;
	unsigned int Nn3_WS = 10, Ns3_WS = 4;
	unsigned int symms3_WS[4] = { 3, 3, 3, 1 };
	double rst3_WS[20] = { -0.338677036009830,  -0.455664103498571,
	                        0.563955207227339,  -0.065470865113645,
	                       -0.225278171217509,   0.521134968612215,
	                       -0.563955207227339,  -0.065470865113645,
	                        0.338677036009830,  -0.455664103498571,
	                        0.225278171217509,   0.521134968612215,
	                       -0.833307841990621,  -0.481110506891111,
	                        0.833307841990621,  -0.481110506891111,
	                       -0.000000000000000,   0.962221013782223,
	                       -0.000000000000000,  -0.000000000000001};
	double w3_WS[10] = { 0.194160145154569,
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
	mkl_dimatcopy('R','T',Nn3_WS,d,1.0,rst3_WS,d,Nn3_WS);

	cubature_TRI(&rst,&w,&symms,&Nn,&Ns,1,P,d,"WS");

	pass = 0;
	if (array_norm_diff_d(Nn3_WS*d,rst3_WS,rst,"Inf")    < EPS &&
	    array_norm_diff_d(Nn3_WS,w3_WS,w,"Inf")          < EPS &&
	    array_norm_diff_ui(Ns3_WS,symms3_WS,symms,"Inf") < EPS &&
	    Nn3_WS == Nn && Ns3_WS == Ns)
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("             (P3, WS):                           ");
	test_print(pass);

	free(rst), free(w);
	free(symms);

	// WV (P = 6)
	P = 6;

	unsigned int Nn6_WV = 12, Ns6_WV = 4;
	unsigned int symms6_WV[4] = { 3, 3, 3, 3 };
	double rst6_WV[24] = { -0.326150048087614,  -0.485300342687622,
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

	double w6_WV[12] = { 0.143502272432754,
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
	mkl_dimatcopy('R','T',Nn6_WV,d,1.0,rst6_WV,d,Nn6_WV);

	cubature_TRI(&rst,&w,&symms,&Nn,&Ns,1,P,d,"WV");

	pass = 0;
	if (array_norm_diff_d(Nn6_WV*d,rst6_WV,rst,"Inf")    < EPS &&
	    array_norm_diff_d(Nn6_WV,w6_WV,w,"Inf")          < EPS &&
	    array_norm_diff_ui(Ns6_WV,symms6_WV,symms,"Inf") < EPS &&
	    Nn6_WV == Nn && Ns6_WV == Ns)
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("             (P6, WV):                           ");
	test_print(pass);

	free(rst), free(w);
	free(symms);
}
