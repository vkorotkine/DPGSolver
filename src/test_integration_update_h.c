// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "test_integration_update_h.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "mkl.h"

#include "Parameters.h"
#include "Macros.h"
#include "Test.h"
#include "S_DB.h"
#include "S_ELEMENT.h"
#include "S_VOLUME.h"
#include "S_FACE.h"

#include "test_code_integration.h"
#include "test_support.h"
#include "array_norm.h"
#include "array_print.h"
#include "adaptation.h"
#include "element_functions.h"
#include "matrix_functions.h"

/*
 *	Purpose:
 *		Test correctness of implementation of update_VOLUME_hp and update_FACE_hp for h-adaptation. Also ensures that
 *		Jacobians of refined VOLUMEs are positive.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

struct S_Limits {
	unsigned int index; // 0: < (All dims), 1: > (All dims)
	char         type;  // (a)nd, (d)iagonal, (o)r
	double       XYZ[3];
};

static void mark_VOLUMEs(const unsigned int adapt_type, const struct S_Limits *Lmts);
static void check_correspondence(unsigned int *pass);
static void check_Jacobians(unsigned int *pass);
static void run_test(unsigned int *pass, const char *test_type);

static void test_update_h(int nargc, char **argvNew, const unsigned int Nref, const unsigned int update_argv,
                          const char *EName, struct S_Limits **Lmts)
{
	unsigned int TETrefineType = DB.TETrefineType;

	unsigned int pass = 0, refType, NrefTypes;

	if (strstr(EName,"TET"))
		NrefTypes = 3;
	else
		NrefTypes = 1;

	for (refType = 0; refType < NrefTypes; refType++) {
		code_startup_mod_prmtrs(nargc,argvNew,Nref,update_argv,1);
		if      (refType == 0) DB.TETrefineType = TET8;
		else if (refType == 1) DB.TETrefineType = TET12;
		else if (refType == 2) DB.TETrefineType = TET6;
		code_startup_mod_prmtrs(nargc,argvNew,Nref,update_argv,2);
		if (DB.PGlobal <= 1)
			printf("Please increase PGlobal above 1 in the ctrl file (%s.ctrl)\n",argvNew[1]), TestDB.Nwarnings++;

		run_test(&pass,"FullREFINE");
		if (strstr(EName,"TRI"))
			printf("update_h (%s%d FullREFINE):                   ",EName,refType);
		else
			printf("         (%s%d FullREFINE):                   ",EName,refType);
		test_print(pass);

		//     0         10        20        30        40        50
		run_test(&pass,"FullCOARSE");
		printf("         (        FullCOARSE):                   ");
		test_print(pass);

		mark_VOLUMEs(HREFINE,Lmts[0]);
		mark_VOLUMEs(HCOARSE,Lmts[1]);

		mesh_update();

		mark_VOLUMEs(HREFINE,Lmts[2]);
		mark_VOLUMEs(HCOARSE,Lmts[3]);

		//     0         10        20        30        40        50
		run_test(&pass,"Mixed");
		printf("         (        Mixed):                        ");
		test_print(pass);

		pass = 1;
		check_Jacobians(&pass);
		//     0         10        20        30        40        50
		printf("         (        Jacobians):                    ");
		test_print(pass);

		code_cleanup();
	}
	DB.TETrefineType = TETrefineType;
}

void test_integration_update_h(int nargc, char **argv)
{
//	unsigned int pass = 0;
	char         **argvNew;

	argvNew    = malloc(2          * sizeof *argvNew);  // free
	argvNew[0] = malloc(STRLEN_MAX * sizeof **argvNew); // free
	argvNew[1] = malloc(STRLEN_MAX * sizeof **argvNew); // free

	strcpy(argvNew[0],argv[0]);

	/*
	 *	Input:
	 *
	 *		Mesh file to be h-refined/coarsened.
	 *
	 *	Expected Output:
	 *
	 *		Interpolation of VOLUME geometry nodes to maximally 1-irregular FACEs from both sides should match.
	 *
	 */

	unsigned int    i;
	char            *EName;
	struct S_Limits **Lmts;

	EName   = malloc(STRLEN_MIN * sizeof *EName); // free
	Lmts = malloc(4 * sizeof *Lmts); // free
	for (i = 0; i < 4; i++)
		Lmts[i] = malloc(sizeof *Lmts[i]);

	// **************************************************************************************************** //
	// TRIs
	strcpy(argvNew[1],"test/Test_update_h_TRI");
	strcpy(EName,"TRI   ");
	Lmts[0]->XYZ[0] = -0.75; Lmts[0]->XYZ[1] = -0.75; Lmts[0]->type = 'd'; Lmts[0]->index = 0;
	Lmts[1]->XYZ[0] =  0.00; Lmts[1]->XYZ[1] =  0.00; Lmts[1]->type = 'd'; Lmts[1]->index = 1;
	Lmts[2]->XYZ[0] = -0.50; Lmts[2]->XYZ[1] = -0.50; Lmts[2]->type = 'd'; Lmts[2]->index = 0;
	Lmts[3]->XYZ[0] =  0.00; Lmts[3]->XYZ[1] =  0.00; Lmts[3]->type = 'o'; Lmts[3]->index = 1;
	test_update_h(nargc,argvNew,2,0,EName,Lmts);

	// **************************************************************************************************** //
	// QUADs
	strcpy(argvNew[1],"test/Test_update_h_QUAD");
	strcpy(EName,"QUAD  ");
	Lmts[0]->XYZ[0] = -0.50; Lmts[0]->XYZ[1] = -0.50; Lmts[0]->type = 'a'; Lmts[0]->index = 0;
	Lmts[1]->XYZ[0] =  0.00; Lmts[1]->XYZ[1] =  0.00; Lmts[1]->type = 'a'; Lmts[1]->index = 1;
	Lmts[2]->XYZ[0] = -0.25; Lmts[2]->XYZ[1] = -0.25; Lmts[2]->type = 'a'; Lmts[2]->index = 0;
	Lmts[3]->XYZ[0] =  0.00; Lmts[3]->XYZ[1] =  0.00; Lmts[3]->type = 'o'; Lmts[3]->index = 1;
	test_update_h(nargc,argvNew,3,0,EName,Lmts);

	// **************************************************************************************************** //
	// TETs
	strcpy(argvNew[1],"test/Test_update_h_TET");
	strcpy(EName,"TET   ");
	Lmts[0]->XYZ[0] =  1.00; Lmts[0]->XYZ[1] =  1.00; Lmts[0]->XYZ[2] =  0.00; Lmts[0]->type = 'd'; Lmts[0]->index = 1;
	Lmts[1]->XYZ[0] = -1.00; Lmts[1]->XYZ[1] = -1.00; Lmts[1]->XYZ[2] =  1.00; Lmts[1]->type = 'd'; Lmts[1]->index = 0;
	Lmts[2]->XYZ[0] =  1.00; Lmts[2]->XYZ[1] =  1.00; Lmts[2]->XYZ[2] = -0.50; Lmts[2]->type = 'd'; Lmts[2]->index = 1;
	Lmts[3]->XYZ[0] =  1.00; Lmts[3]->XYZ[1] =  1.00; Lmts[3]->XYZ[2] = -1.00; Lmts[3]->type = 'd'; Lmts[3]->index = 0;
	test_update_h(nargc,argvNew,2,0,EName,Lmts);

	// **************************************************************************************************** //
	// HEXs
	strcpy(argvNew[1],"test/Test_update_h_HEX");
	strcpy(EName,"HEX   ");
	Lmts[0]->XYZ[0] = -0.50; Lmts[0]->XYZ[1] = -0.50; Lmts[0]->XYZ[2] = -0.50; Lmts[0]->type = 'a'; Lmts[0]->index = 0;
	Lmts[1]->XYZ[0] =  0.00; Lmts[1]->XYZ[1] =  0.00; Lmts[1]->XYZ[2] =  0.00; Lmts[1]->type = 'a'; Lmts[1]->index = 1;
	Lmts[2]->XYZ[0] = -0.25; Lmts[2]->XYZ[1] = -0.25; Lmts[2]->XYZ[2] = -0.25; Lmts[2]->type = 'a'; Lmts[2]->index = 0;
	Lmts[3]->XYZ[0] =  0.00; Lmts[3]->XYZ[1] =  0.00; Lmts[3]->XYZ[2] =  0.00; Lmts[3]->type = 'o'; Lmts[3]->index = 1;
	test_update_h(nargc,argvNew,3,0,EName,Lmts);

	// **************************************************************************************************** //
	// WEDGEs
	strcpy(argvNew[1],"test/Test_update_h_WEDGE");
	strcpy(EName,"WEDGE ");
	Lmts[0]->XYZ[0] = -0.50; Lmts[0]->XYZ[1] = -0.50; Lmts[0]->XYZ[2] = -0.50; Lmts[0]->type = 'a'; Lmts[0]->index = 0;
	Lmts[1]->XYZ[0] =  0.00; Lmts[1]->XYZ[1] =  0.00; Lmts[1]->XYZ[2] =  0.00; Lmts[1]->type = 'a'; Lmts[1]->index = 1;
	Lmts[2]->XYZ[0] = -0.25; Lmts[2]->XYZ[1] = -0.25; Lmts[2]->XYZ[2] = -0.25; Lmts[2]->type = 'a'; Lmts[2]->index = 0;
	Lmts[3]->XYZ[0] =  0.00; Lmts[3]->XYZ[1] =  0.00; Lmts[3]->XYZ[2] =  0.00; Lmts[3]->type = 'o'; Lmts[3]->index = 1;
	test_update_h(nargc,argvNew,3,0,EName,Lmts);

	// **************************************************************************************************** //
	// PYRs
	strcpy(argvNew[1],"test/Test_update_h_PYR");
	strcpy(EName,"PYR   ");
	Lmts[0]->XYZ[0] = -0.50; Lmts[0]->XYZ[1] = -0.50; Lmts[0]->XYZ[2] = -0.50; Lmts[0]->type = 'a'; Lmts[0]->index = 0;
	Lmts[1]->XYZ[0] =  0.00; Lmts[1]->XYZ[1] =  0.00; Lmts[1]->XYZ[2] =  0.00; Lmts[1]->type = 'a'; Lmts[1]->index = 1;
	Lmts[2]->XYZ[0] = -0.25; Lmts[2]->XYZ[1] = -0.25; Lmts[2]->XYZ[2] = -0.25; Lmts[2]->type = 'a'; Lmts[2]->index = 0;
	Lmts[3]->XYZ[0] =  0.00; Lmts[3]->XYZ[1] =  0.00; Lmts[3]->XYZ[2] =  0.00; Lmts[3]->type = 'o'; Lmts[3]->index = 1;
	test_update_h(nargc,argvNew,2,0,EName,Lmts);

//	output_to_paraview("ZTest_Geomadapt");
//	output_to_paraview("ZTest_Normals");
//	EXIT_MSG;

	free(argvNew[0]); free(argvNew[1]); free(argvNew);
	free(EName);
	for (i = 0; i < 4; i++)
		free(Lmts[i]);
	free(Lmts);
}

static void mark_VOLUMEs(const unsigned int adapt_type, const struct S_Limits *Lmts)
{
	// Initialize DB Parameters
	unsigned int d = DB.d;

	// Standard datatypes
	unsigned int i, dim, Nve, IndXYZ, update;
	double       XYZ_cent[3], *XYZ_vV, Lmts_sum, XYZ_cent_sum;

	struct S_ELEMENT *ELEMENT;
	struct S_VOLUME  *VOLUME;

	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		// Calculate centroid
		ELEMENT = get_ELEMENT_type(VOLUME->type);
		Nve = ELEMENT->Nve;

		XYZ_vV = VOLUME->XYZ_vV;

		for (i = 0; i < 3; i++)
			XYZ_cent[i] = 0.0;

		for (dim = 0; dim < d; dim++) {
			IndXYZ = dim*Nve;
			for (i = 0; i < Nve; i++)
				XYZ_cent[dim] += XYZ_vV[IndXYZ+i];
			XYZ_cent[dim] /= Nve;
		}

		// Mark VOLUME if centroid is within limits
		switch (Lmts->type) {
		case 'a':
			update = 1;
			for (dim = 0; dim < d; dim++) {
				switch (Lmts->index) {
				case 0: if (!(XYZ_cent[dim] < Lmts->XYZ[dim])) update = 0; break;
				case 1: if (!(XYZ_cent[dim] > Lmts->XYZ[dim])) update = 0; break;
				default:
					printf("Error: Unsupported index in mark_VOLUMEs for type (a).\n"), EXIT_MSG;
					break;
				}
			}
			break;
		case 'o':
			update = 0;
			for (dim = 0; dim < d; dim++) {
				switch (Lmts->index) {
				case 0: if (XYZ_cent[dim] < Lmts->XYZ[dim]) update = 1; break;
				case 1: if (XYZ_cent[dim] > Lmts->XYZ[dim]) update = 1; break;
				default:
					printf("Error: Unsupported index in mark_VOLUMEs for type (o).\n"), EXIT_MSG;
					break;
				}
			}
			break;
		case 'd':
			update = 1;
			Lmts_sum  = 0.0;
			XYZ_cent_sum = 0.0;
			for (dim = 0; dim < d; dim++) {
				Lmts_sum  += Lmts->XYZ[dim];
				XYZ_cent_sum += XYZ_cent[dim];
			}
			switch (Lmts->index) {
			case 0: if (!(XYZ_cent_sum < Lmts_sum)) update = 0; break;
			case 1: if (!(XYZ_cent_sum > Lmts_sum)) update = 0; break;
			default:
				printf("Error: Unsupported index in mark_VOLUMEs for type (d).\n"), EXIT_MSG;
				break;
			}
			break;
		default:
			printf("Error: Unsupported type in mark_VOLUMEs.\n"), EXIT_MSG;
			break;
		}

		if (update) {
			VOLUME->Vadapt = 1;
			VOLUME->adapt_type = adapt_type;
		}
	}
}


struct S_OPERATORS {
	unsigned int NfnS, NvnGs, *nOrdInOut, *nOrdOutIn;
	double       **I_vGs_fS;
};

static void init_ops(struct S_OPERATORS *OPS, const struct S_VOLUME *VOLUME, const struct S_FACE *FACE,
                     const unsigned int IndFType)
{
	unsigned int PF, VType, IndOrdInOut, IndOrdOutIn;

	struct S_ELEMENT *ELEMENT, *ELEMENT_FACE;

//	PV    = VOLUME->P;
	VType = VOLUME->type;

	PF = FACE->P;
	IndOrdInOut = FACE->IndOrdInOut;
	IndOrdOutIn = FACE->IndOrdOutIn;

	ELEMENT = get_ELEMENT_type(VType);
	ELEMENT_FACE = get_ELEMENT_FACE(VType,IndFType);

	OPS->NvnGs = ELEMENT->NvnGs[1];
	OPS->NfnS  = ELEMENT_FACE->NvnS[PF];

	OPS->I_vGs_fS = ELEMENT->I_vGs_fS[1][PF];

	OPS->nOrdInOut = ELEMENT_FACE->nOrd_fS[PF][IndOrdInOut];
	OPS->nOrdOutIn = ELEMENT_FACE->nOrd_fS[PF][IndOrdOutIn];
}

static void check_correspondence(unsigned int *pass)
{
	// Initialize DB Parameters
	unsigned int d = DB.d;

	// Standard datatypes
	unsigned int Vf, IndFType, NfnS, *nOrdInOut, *nOrdOutIn, vhIn, vhOut,
	             dim, n, Indd, BC, FACE_is_internal;
	double       *XYZ_fSIn, *XYZ_fSOut, *XYZ_fSInOut, *XYZ_fSOutIn;

	struct S_OPERATORS *OPS;
	struct S_VOLUME    *VOLUME, *VOLUMEc;
	struct S_FACE     *FACE;

	OPS = malloc(sizeof *OPS); // free

	*pass = 1;
	for (FACE = DB.FACE; FACE; FACE = FACE->next) {
		VOLUME = FACE->VIn;
		Vf     = FACE->VfIn;

		IndFType = get_IndFType(VOLUME->Eclass,Vf/NFREFMAX);
		init_ops(OPS,VOLUME,FACE,IndFType);

		NfnS = OPS->NfnS;

		nOrdInOut = OPS->nOrdInOut;
		nOrdOutIn = OPS->nOrdOutIn;

		XYZ_fSIn = mm_Alloc_d(CBCM,CBT,CBNT,NfnS,d,OPS->NvnGs,1.0,OPS->I_vGs_fS[Vf],VOLUME->XYZ_vV); // free

		VOLUME = FACE->VOut;
		Vf     = FACE->VfOut;

		IndFType = get_IndFType(VOLUME->Eclass,Vf/NFREFMAX);
		init_ops(OPS,VOLUME,FACE,IndFType);

		XYZ_fSOut = mm_Alloc_d(CBCM,CBT,CBNT,NfnS,d,OPS->NvnGs,1.0,OPS->I_vGs_fS[Vf],VOLUME->XYZ_vV); // free

		XYZ_fSInOut = malloc(NfnS*d * sizeof *XYZ_fSInOut); // free
		XYZ_fSOutIn = malloc(NfnS*d * sizeof *XYZ_fSOutIn); // free

		for (dim = 0; dim < d; dim++) {
			Indd = dim*NfnS;
			for (n = 0; n < NfnS; n++) {
				XYZ_fSInOut[Indd+n] = XYZ_fSIn[Indd+nOrdInOut[n]];
				XYZ_fSOutIn[Indd+n] = XYZ_fSOut[Indd+nOrdOutIn[n]];
			}
		}

		BC = FACE->BC;
		FACE_is_internal = (BC == 0 || (BC % BC_STEP_SC > 50));

		if (FACE_is_internal && (array_norm_diff_d(NfnS*d,XYZ_fSIn,XYZ_fSOutIn,"Inf")  > 10*EPS ||
		                          array_norm_diff_d(NfnS*d,XYZ_fSInOut,XYZ_fSOut,"Inf") > 10*EPS)) {
				*pass = 0;
				printf("Problem in check_correspondence\n");

				vhIn = 0;
				for (VOLUMEc = FACE->VIn->parent->child0; VOLUMEc != FACE->VIn; VOLUMEc = VOLUMEc->next)
					vhIn++;

				vhOut = 0;
				for (VOLUMEc = FACE->VOut->parent->child0; VOLUMEc != FACE->VOut; VOLUMEc = VOLUMEc->next)
					vhOut++;

printf("%d %d %d %d %d\n",FACE->indexg,FACE->IndOrdInOut,FACE->IndOrdOutIn,vhIn,vhOut);
printf("%d %d %d %d\n",FACE->VIn->type,FACE->VIn->indexg,FACE->VfIn,FACE->VIn->level);
printf("%d %d %d %d\n",FACE->VOut->type,FACE->VOut->indexg,FACE->VfOut,FACE->VOut->level);
				printf("Errors: %e %e\n\n",array_norm_diff_d(NfnS*d,XYZ_fSIn,XYZ_fSOutIn,"Inf"),
		                                   array_norm_diff_d(NfnS*d,XYZ_fSInOut,XYZ_fSOut,"Inf"));
				array_print_d(NfnS,d,XYZ_fSIn,'C');
				array_print_d(NfnS,d,XYZ_fSOutIn,'C');
				array_print_d(NfnS,d,XYZ_fSOut,'C');
				array_print_d(NfnS,d,XYZ_fSInOut,'C');
EXIT_MSG;
				break;
		}

		free(XYZ_fSIn);
		free(XYZ_fSOut);
		free(XYZ_fSInOut);
		free(XYZ_fSOutIn);
	}
	free(OPS);

	if (*pass)
		TestDB.Npass++;
}

static void check_Jacobians(unsigned int *pass)
{
	unsigned int i, P, NvnI;
	double       *detJV_vI;

	struct S_ELEMENT *ELEMENT;
	struct S_VOLUME  *VOLUME;

	for (VOLUME = DB.VOLUME; *pass && VOLUME; VOLUME = VOLUME->next) {
		P = VOLUME->P;

		ELEMENT = get_ELEMENT_type(VOLUME->type);

		if (!VOLUME->curved)
			NvnI = ELEMENT->NvnIs[P];
		else
			NvnI = ELEMENT->NvnIc[P];

		detJV_vI = VOLUME->detJV_vI;
		for (i = 0; i < NvnI; i++) {
			if (detJV_vI[i] <= 0.0) {
				*pass = 0;
				break;
			}
		}
	}

	if (*pass)
		TestDB.Npass++;
}

static void run_test(unsigned int *pass, const char *test_type)
{
	struct S_VOLUME *VOLUME;

	if (strstr(test_type,"FullREFINE")) {
		for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
			VOLUME->Vadapt = 1;
			VOLUME->adapt_type = HREFINE;
		}
	} else if (strstr(test_type,"FullCOARSE")) {
		for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
			VOLUME->Vadapt = 1;
			VOLUME->adapt_type = HCOARSE;
		}
	} else {
		// VOLUME processing done outside of this function.
	}
	mesh_update();
	check_correspondence(pass);
}
