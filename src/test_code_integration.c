// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "test_code_integration.h"

#include <string.h>

#include "petscsys.h"

#include "Parameters.h"
#include "Macros.h"
#include "Test.h"
#include "S_DB.h"

#include "initialization.h"
#include "setup_parameters.h"
#include "setup_mesh.h"
#include "setup_operators.h"
#include "setup_structures.h"
#include "setup_geometry.h"
#include "initialize_test_case.h"
#include "adaptation.h"
#include "memory_free.h"

#include "array_print.h"

/*
 *	Purpose:
 *		Provide functions for integration testing.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

static void update_MeshFile(void)
{
	// Standard datatypes
	char *d, *ML;

	d  = malloc(STRLEN_MIN * sizeof *d);  // free
	ML = malloc(STRLEN_MIN * sizeof *ML); // free

	sprintf(d,"%d",DB.d);
	sprintf(ML,"%d",DB.ML);

	strcpy(DB.MeshFile,"");
	strcat(DB.MeshFile,DB.MeshPath);
	strcat(DB.MeshFile,DB.TestCase);
	strcat(DB.MeshFile,"/");
	strcat(DB.MeshFile,DB.TestCase);
	strcat(DB.MeshFile,strcat(d,"D_"));
	strcat(DB.MeshFile,DB.MeshType);
	strcat(DB.MeshFile,strcat(ML,"x.msh"));

	free(d);
	free(ML);
}

void code_startup(int nargc, char **argv, const unsigned int Nref, const unsigned int update_argv)
{
	int  MPIrank, MPIsize;

	// Start MPI and PETSC
	PetscInitialize(&nargc,&argv,PETSC_NULL,PETSC_NULL);
	MPI_Comm_size(MPI_COMM_WORLD,&MPIsize);
	MPI_Comm_rank(MPI_COMM_WORLD,&MPIrank);

	// Test memory leaks only from Petsc and MPI using valgrind
	//PetscFinalize(), exit(1);

	DB.MPIsize = MPIsize;
	DB.MPIrank = MPIrank;

	// Initialization
	initialization(nargc,argv);
	if (update_argv) {
		strcpy(DB.TestCase,TestDB.TestCase);
		DB.PGlobal = TestDB.PGlobal;
		DB.ML      = TestDB.ML;
		update_MeshFile();
	}

	setup_parameters();
	if (update_argv)
		setup_parameters_L2proj();

	setup_mesh();
	setup_operators();
	setup_structures();
	setup_geometry();

	initialize_test_case(Nref);
}

void code_cleanup(void)
{
	mesh_to_level(0);
	memory_free();
}

void check_convergence_orders(const unsigned int MLMin, const unsigned int MLMax, const unsigned int PMin,
                              const unsigned int PMax, unsigned int *pass)
{
	// Initialize DB Parameters
	char         *TestCase = TestDB.TestCase,
	             *MeshType = DB.MeshType;
	unsigned int d         = DB.d;

	// Standard datatypes
	char         f_name[STRLEN_MAX], string[STRLEN_MAX], StringRead[STRLEN_MAX], *data;
	unsigned int i, ML, P, NVars, NML, NP, Indh, *VarsToCheck, *OrderIncrement;
	int          offset;
	double       **L2Errors, **ConvOrders, *h, tmp_d;

	FILE *fID;

	// silence
	NVars = 0;

	if (strstr(TestCase,"Poisson")) {
		NVars = DMAX+1;
	} else if (strstr(TestCase,"SupersonicVortex")) {
		NVars = DMAX+2+1;
	} else {
		printf("Error: Unsupported TestCase.\n"), EXIT_MSG;
	}

	VarsToCheck    = malloc(NVars * sizeof *VarsToCheck);    // free
	OrderIncrement = malloc(NVars * sizeof *OrderIncrement); // free
	if (strstr(TestCase,"Poisson")) {
		for (i = 0; i < NVars; i++) {
			OrderIncrement[i] = 0;
			if (i == 0) {
				OrderIncrement[i] = 1;
				VarsToCheck[i]    = 1;
			} else if (i <= d) {
				VarsToCheck[i]    = 1;
			} else {
				VarsToCheck[i]    = 0;
			}
		}
	} else if (strstr(TestCase,"SupersonicVortex")) {
		for (i = 0; i < NVars; i++) {
			OrderIncrement[i] = 1;
			if (i <= d || i > DMAX) {
				VarsToCheck[i] = 1;
			} else {
				VarsToCheck[i] = 0;
			}
		}
	}


	NML = MLMax-MLMin+1;
	NP  = PMax-PMin+1;

	L2Errors   = malloc(NVars * sizeof *L2Errors);   // free
	ConvOrders = malloc(NVars * sizeof *ConvOrders); // free
	for (i = 0; i < NVars; i++) {
		L2Errors[i]   = calloc(NML*NP , sizeof *L2Errors[i]);   // free
		ConvOrders[i] = calloc(NML*NP , sizeof *ConvOrders[i]); // free
	}
	h = calloc(NML*NP , sizeof *h); // free

	// Read in data and compute convergence orders
	for (ML = MLMin; ML <= MLMax; ML++) {
	for (P = PMin; P <= PMax; P++) {
		Indh = (ML-MLMin)*NP+(P-PMin);

		strcpy(f_name,"../cases/results/");
		strcat(f_name,TestCase); strcat(f_name,"/");
		strcat(f_name,MeshType); strcat(f_name,"/");
		strcat(f_name,"L2errors_");
		sprintf(string,"%dD_",d);   strcat(f_name,string);
		                            strcat(f_name,MeshType);
		sprintf(string,"_ML%d",ML); strcat(f_name,string);
		sprintf(string,"P%d",P);    strcat(f_name,string);
		strcat(f_name,".txt");

		if ((fID = fopen(f_name,"r")) == NULL)
			printf("Error: File: %s, did not open.\n",f_name), exit(1);

		if (fscanf(fID,"%[^\n]\n",StringRead) == 1) { ; }
		if (fscanf(fID,"%[^\n]\n",StringRead) == 1) {
			i = 0;
			data = StringRead;
			if (sscanf(data," %lf%n",&tmp_d,&offset) == 1) {
				data += offset;
				h[Indh] = 1.0/pow(tmp_d,1.0/d);
			}
			while (sscanf(data," %lf%n",&tmp_d,&offset) == 1) {
				L2Errors[i++][Indh] = tmp_d;
				data += offset;
			}
		}
		fclose(fID);
	}}

	for (ML = MLMin+1; ML <= MLMax; ML++) {
	for (P = PMin; P <= PMax; P++) {
		Indh = (ML-MLMin)*NP+(P-PMin);
		for (i = 0; i < NVars; i++) {
			if (fabs(L2Errors[i][Indh]) > EPS && fabs(L2Errors[i][Indh-NP]) > EPS)
				ConvOrders[i][Indh] = log10(L2Errors[i][Indh]/L2Errors[i][Indh-NP])/log10(h[Indh]/h[Indh-NP]);
		}
	}}

	*pass = 1;

	ML = MLMax;
	for (i = 0; i < NVars; i++) {
		if (!VarsToCheck[i])
			continue;

		for (P = PMin; P <= PMax; P++) {
			Indh = (ML-MLMin)*NP+(P-PMin);
			if (fabs(ConvOrders[i][Indh]-(P+OrderIncrement[i])) > 0.125) {
				*pass = 0;

				printf("i = %d, P = %d\n",i,P);
				break;
			}
		}
	}

	free(VarsToCheck);
	free(OrderIncrement);

	if (!(*pass)) {
		for (i = 0; i < NVars; i++)
			array_print_d(NML,NP,ConvOrders[i],'R');
	} else {
		TestDB.Npass++;
	}

	for (i = 0; i < NVars; i++) {
		free(L2Errors[i]);
		free(ConvOrders[i]);
	}
	free(L2Errors);
	free(ConvOrders);
	free(h);
}
