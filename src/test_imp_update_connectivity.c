// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <mpi.h>
#include <petscksp.h>

#include "test.h"
#include "parameters.h"
#include "functions.h"
#include "database.h"

/*
 *	Purpose:
 *		Test correctness of implementation of update_connectivity.
 *
 *	Comments:
 *		Ensure to test with BCs other than periodic as well. (ToBeDeleted)
 *		Fix memory leaks. (ToBeDeleted)
 *
 *	Notation:
 *
 *	References:
 */

struct S_OPERATORS {
	unsigned int NfnS, NvnGs, *nOrdInOut, *nOrdOutIn;
	double       **I_vGs_fS;
};

static void init_ops(struct S_OPERATORS *OPS, const struct S_VOLUME *VOLUME, const struct S_FACET *FACET,
                     const unsigned int IndFType)
{
	unsigned int PV, PF, VType, IndOrdInOut, IndOrdOutIn;

	struct S_ELEMENT *ELEMENT, *ELEMENT_FACET;

	PV    = VOLUME->P;
	VType = VOLUME->type;

	PF = FACET->P;
	IndOrdInOut = FACET->IndOrdInOut;
	IndOrdOutIn = FACET->IndOrdOutIn;

	ELEMENT = get_ELEMENT_type(VType);
	ELEMENT_FACET = get_ELEMENT_FACET(VType,IndFType);

	OPS->NvnGs = ELEMENT->NvnGs[1];
	OPS->NfnS  = ELEMENT_FACET->NvnS[PF];

	OPS->I_vGs_fS = ELEMENT->I_vGs_fS[PV][PF];

	OPS->nOrdInOut = ELEMENT_FACET->nOrd_fS[PF][IndOrdInOut];
	OPS->nOrdOutIn = ELEMENT_FACET->nOrd_fS[PF][IndOrdOutIn];
}

static void code_startup(int nargc, char **argv)
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
	setup_parameters();
	setup_mesh();
	setup_operators();
	setup_structures();
	setup_geometry();
}

static void check_correspondence(unsigned int *pass)
{
	// Initialize DB Parameters
	unsigned int d = DB.d;

	// Standard datatypes
	unsigned int Vf, IndFType, NfnS, *nOrdInOut, *nOrdOutIn,
	             dim, n, Indd, BC, FACET_is_internal;
	double       *XYZ_fSIn, *XYZ_fSOut, *XYZ_fSInOut, *XYZ_fSOutIn;

	struct S_OPERATORS *OPS;
	struct S_VOLUME    *VOLUME;
	struct S_FACET     *FACET;

	OPS = malloc(sizeof *OPS); // free

	*pass = 1;
	for (FACET = DB.FACET; FACET != NULL; FACET = FACET->next) {
		VOLUME = FACET->VIn;
		Vf     = FACET->VfIn;

		IndFType = get_IndFType(VOLUME->Eclass,Vf/NFREFMAX);
		init_ops(OPS,VOLUME,FACET,IndFType);

		NfnS = OPS->NfnS;

		nOrdInOut = OPS->nOrdInOut;
		nOrdOutIn = OPS->nOrdOutIn;

		XYZ_fSIn = mm_Alloc_d(CBCM,CBT,CBNT,NfnS,d,OPS->NvnGs,1.0,OPS->I_vGs_fS[Vf],VOLUME->XYZ_vC);

		VOLUME = FACET->VOut;
		Vf     = FACET->VfOut;

		IndFType = get_IndFType(VOLUME->Eclass,Vf/NFREFMAX);
		init_ops(OPS,VOLUME,FACET,IndFType);
		XYZ_fSOut = mm_Alloc_d(CBCM,CBT,CBNT,NfnS,d,OPS->NvnGs,1.0,OPS->I_vGs_fS[Vf],VOLUME->XYZ_vC);

		XYZ_fSInOut = malloc(NfnS*d * sizeof *XYZ_fSInOut); // free
		XYZ_fSOutIn = malloc(NfnS*d * sizeof *XYZ_fSOutIn); // free

		for (dim = 0; dim < d; dim++) {
			Indd = dim*NfnS;
			for (n = 0; n < NfnS; n++) {
				XYZ_fSInOut[Indd+n] = XYZ_fSIn[Indd+nOrdInOut[n]];
				XYZ_fSOutIn[Indd+n] = XYZ_fSOut[Indd+nOrdOutIn[n]];
			}
		}

		BC = FACET->BC;
		FACET_is_internal = (BC == 0 || (BC % BC_STEP_SC > 50));

		if (FACET_is_internal && (array_norm_diff_d(NfnS*d,XYZ_fSIn,XYZ_fSOutIn,"Inf")  > EPS ||
		                          array_norm_diff_d(NfnS*d,XYZ_fSInOut,XYZ_fSOut,"Inf") > EPS)) {
				*pass = 0;
				printf("Problem in check_correspondence\n");
printf("%d\n",FACET->indexg);
printf("%d %d\n",FACET->VIn->indexg,FACET->VfIn);
printf("%d %d\n",FACET->VOut->indexg,FACET->VfOut);
				printf("Errors: %e %e\n\n",array_norm_diff_d(NfnS*d,XYZ_fSIn,XYZ_fSOutIn,"Inf"),
		                                   array_norm_diff_d(NfnS*d,XYZ_fSInOut,XYZ_fSOut,"Inf"));
				array_print_d(NfnS,d,XYZ_fSIn,'C');
				array_print_d(NfnS,d,XYZ_fSOutIn,'C');
				array_print_d(NfnS,d,XYZ_fSOut,'C');
				array_print_d(NfnS,d,XYZ_fSInOut,'C');
				break;
		}

		free(XYZ_fSInOut);
		free(XYZ_fSOutIn);
	}
	free(OPS);

	if (*pass)
		TestDB.Npass++;
}

static void renumber_VOLUMEs(void)
{
	unsigned int NV = 0;

	struct S_VOLUME *VOLUME;

	for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next)
		VOLUME->indexg = NV++;
	
	DB.NV = NV;
}

void test_imp_update_connectivity(int nargc, char **argv)
{
	unsigned int pass;

	/*
	 *	Input:
	 *
	 *		ToBeModified
	 *
	 *	Expected Output (d = 1, Trivial):
	 *
	 *		ToBeModified
	 *
	 */

	char         **argvNew;
	unsigned int indexg;

	struct S_VOLUME *VOLUME;
	struct S_FACET  *FACET;

	// LINEs


	// TRIs
	argvNew = malloc(2 * sizeof *argvNew);
	argvNew[0] = malloc(STRLEN_MAX * sizeof **argvNew);
	argvNew[1] = malloc(STRLEN_MAX * sizeof **argvNew);

	strcpy(argvNew[0],argv[0]);
	strcpy(argvNew[1],"test/TestTRI");

	code_startup(nargc,argvNew);
	output_to_paraview("ZTest_GeomTRI");

	// Mark VOLUMEs for isotropic h-refinement
	for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next) {
		indexg = VOLUME->indexg;

		if (indexg == 0 || indexg == 3) {
			VOLUME->Vadapt = 1;
			VOLUME->adapt_type = HREFINE;
			VOLUME->hrefine_type = 0;
		}
	}

	update_VOLUME_hp();
	update_FACET_hp();

	update_VOLUME_list();
	renumber_VOLUMEs();

	for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next) {
		printf("VOL: %d %d ",VOLUME->indexg,VOLUME->level);
		if (VOLUME->level > 0)
			printf("%d ",VOLUME->parent->indexg);
		printf("\n");
	}

	printf("\n\n\n");
	for (FACET = DB.FACET; FACET != NULL; FACET = FACET->next) {
		printf("%d %d   ",FACET->indexg,FACET->update);
		if (FACET->level > 0 && FACET->parent)
			printf("%d ",FACET->parent->indexg);
		printf("\n");
	}

	output_to_paraview("ZTest_GeomTRI_up");
	output_to_paraview("ZTest_NormalsTRI_up");

printf("\n\n\n");
//exit(1);

	//     0         10        20        30        40        50
	check_correspondence(&pass);
	printf("update_connectivity (d = 2, TRI):                ");
	test_print(pass);
//exit(1);

	for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next) {
//printf("%d\n",VOLUME->indexg);
//array_print_ui(3,NFREFMAX,VOLUME->neigh_f,'R');
		VOLUME->Vadapt = 1;
		VOLUME->adapt_type = HREFINE;
		VOLUME->hrefine_type = 0;
	}
//exit(1);

	update_VOLUME_hp();
	update_FACET_hp();

	update_VOLUME_list();
	renumber_VOLUMEs();

	output_to_paraview("ZTest_GeomTRI_up");
	output_to_paraview("ZTest_NormalsTRI_up");

/*
printf("\n\n\n");
	for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next) {
		printf("VOL: %d %d ",VOLUME->indexg,VOLUME->level);
		if (VOLUME->level > 0)
			printf("%d ",VOLUME->parent->indexg);
		printf("\n");
	}
*/

	//     0         10        20        30        40        50
	check_correspondence(&pass);
	printf("update_connectivity (d = 2, TRI):                ");
	test_print(pass);


/*
for (VOLUME = DB.VOLUME; VOLUME->indexg != 5; VOLUME = VOLUME->next)
	;
for (int i = 0; i < 3; i++) {
	printf("NsubF %d",VOLUME->NsubF[i]);
	for (int j = 0; j < VOLUME->NsubF[i]; j++)
		printf(" %d",VOLUME->FACET[i*NSUBFMAX+j]->indexg);
	printf("\n");
}
*/

exit(1);

/*
	pass = 0;
	if (array_norm_diff_ui(Nn,nOrd,nOrd10,"Inf") < EPS)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("get_facet_ordering (d = 1, case 0):              ");
	test_print(pass);
*/

}
