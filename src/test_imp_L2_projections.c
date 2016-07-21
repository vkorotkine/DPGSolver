// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include <petscksp.h>

#include "test.h"
#include "parameters.h"
#include "functions.h"
#include "database.h"

/*
 *	Purpose:
 *		Test correctness of implementation of L2 projection operators.
 *
 *	Comments:
 *		Fix memory leaks. (ToBeDeleted)
 *
 *	Notation:
 *
 *	References:
 */

static void   code_startup  (int nargc, char **argv, const unsigned int Nref);
static void   code_cleanup  (const unsigned int final);
static double *get_L2err    (void);
static void   mark_VOLUMEs  (const unsigned int adapt_type);
static void   mesh_update   (void);
static void   mesh_to_level (const unsigned int level);

void test_imp_L2_projections(int nargc, char **argv)
{
	unsigned int pass;
	char         **argvNew;
	double       *L2err[2];

	argvNew    = malloc(2          * sizeof *argvNew);  // free
	argvNew[0] = malloc(STRLEN_MAX * sizeof **argvNew); // free
	argvNew[1] = malloc(STRLEN_MAX * sizeof **argvNew); // free

	/*
	 *	Input:
	 *
	 *		Mesh file to be refined/coarsened.
	 *
	 *	Expected Output:
	 *
	 *		The L2 error of the initial solution after refinement and coarsening should not change.
	 *
	 */

	// LINEs


	// TRIs
	strcpy(argvNew[0],argv[0]);
	strcpy(argvNew[1],"test/Test_L2_proj_TRI");

	code_startup(nargc,argvNew,2);

	L2err[0] = get_L2err();
	mark_VOLUMEs(PREFINE); mesh_update();
	mark_VOLUMEs(PCOARSE); mesh_update();
	L2err[1] = get_L2err();

	pass = 0;
	if (array_norm_diff_d(NVAR3D+1,L2err[0],L2err[1],"Inf") < 10*EPS)
		pass = 1, TestDB.Npass++;
	//     0         10        20        30        40        50
	printf("L2_projections (TRI, ADAPT_P):                   ");
	test_print(pass);

	code_cleanup(0);

/*
	strcpy(argvNew[0],argv[0]);
	strcpy(argvNew[1],"test/Test_L2_proj_TRI");

	code_startup(nargc,argvNew,2);

	L2err[0] = get_L2err();
	mark_VOLUMEs(HREFINE); mesh_update();
	mark_VOLUMEs(HCOARSE); mesh_update();
	L2err[1] = get_L2err();

	pass = 0;
	if (array_norm_diff_d(1,&L2err[0],&L2err[1],"Inf") < EPS)
		pass = 1, TestDB.Npass++;
	//     0         10        20        30        40        50
	printf("L2_projections (TRI, ADAPT_H):                   ");
	test_print(pass);

	code_cleanup(0);
*/

	// QUADs
	strcpy(argvNew[0],argv[0]);
	strcpy(argvNew[1],"test/Test_L2_proj_QUAD");

//	code_startup(nargc,argvNew,3);


	free(argvNew[0]); free(argvNew[1]); free(argvNew);
}

static void code_startup(int nargc, char **argv, const unsigned int Nref)
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

	initialize_test_case(Nref);
}

static void code_cleanup(const unsigned int final)
{
	memory_free();
	if (final)
		PetscFinalize();
}



struct S_OPERATORS {
	unsigned int NfnS, NvnGs, *nOrdInOut, *nOrdOutIn;
	double       **I_vGs_fS;
};

static void init_ops(struct S_OPERATORS *OPS, const struct S_VOLUME *VOLUME, const struct S_FACET *FACET,
                     const unsigned int IndFType)
{
	unsigned int PF, VType, IndOrdInOut, IndOrdOutIn;

	struct S_ELEMENT *ELEMENT, *ELEMENT_FACET;

//	PV    = VOLUME->P;
	VType = VOLUME->type;

	PF = FACET->P;
	IndOrdInOut = FACET->IndOrdInOut;
	IndOrdOutIn = FACET->IndOrdOutIn;

	ELEMENT = get_ELEMENT_type(VType);
	ELEMENT_FACET = get_ELEMENT_FACET(VType,IndFType);

	OPS->NvnGs = ELEMENT->NvnGs[1];
	OPS->NfnS  = ELEMENT_FACET->NvnS[PF];

	OPS->I_vGs_fS = ELEMENT->I_vGs_fS[1][PF];

	OPS->nOrdInOut = ELEMENT_FACET->nOrd_fS[PF][IndOrdInOut];
	OPS->nOrdOutIn = ELEMENT_FACET->nOrd_fS[PF][IndOrdOutIn];
}

static double *get_L2err(void)
{
	unsigned int i, iMax, dummy_ui;
	double       dummy_d, *L2Error2, *L2Error;

	struct S_VOLUME *VOLUME;

	iMax = NVAR3D+1;
	L2Error2 = malloc(iMax * sizeof *L2Error2); // free
	L2Error  = calloc(iMax , sizeof *L2Error);  // keep (requires external free)
	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		compute_errors(VOLUME,L2Error2,&dummy_d,&dummy_ui,0);

//		printf("%d\n",VOLUME->indexg);
//		array_print_d(1,NVAR3D+1,L2Error2,'R');

		for (i = 0; i < iMax; i++)
			L2Error[i] += sqrt(L2Error2[i]);
	}
	array_print_d(1,iMax,L2Error,'R');

	free(L2Error2);
	return L2Error;
}

static void mark_VOLUMEs(const unsigned int adapt_type)
{
	struct S_VOLUME *VOLUME;
	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		VOLUME->Vadapt = 1;
		VOLUME->adapt_type = adapt_type;
	}
}

static void mesh_update(void)
{
	update_VOLUME_hp();
	update_FACET_hp();
	update_VOLUME_list();
	memory_free_children();
	update_VOLUME_finalize();
}

static void mesh_to_level(const unsigned int level)
{
	unsigned int updated = 1;
	struct S_VOLUME *VOLUME;

	while (updated) {
		updated = 0;
		for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
			if (VOLUME->level != level) {
				updated = 1;

				VOLUME->Vadapt = 1;
				if (VOLUME->level > level)
					VOLUME->adapt_type = HCOARSE;
				else
					VOLUME->adapt_type = HREFINE;
			}
		}

		if (updated)
			mesh_update();
	}
}
