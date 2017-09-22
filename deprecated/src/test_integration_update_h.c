// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

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
#include "output_to_paraview.h"

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

	char *PrintName = malloc(STRLEN_MAX * sizeof *PrintName); // free

	if (strstr(EName,"TET"))
		NrefTypes = 3;
	else
		NrefTypes = 1;

	for (refType = 0; refType < NrefTypes; refType++) {
		code_startup_mod_prmtrs(nargc,(char const *const *const) argvNew,Nref,update_argv,1);
		if      (refType == 0) DB.TETrefineType = TET8;
		else if (refType == 1) DB.TETrefineType = TET12;
		else if (refType == 2) DB.TETrefineType = TET6;
		code_startup_mod_prmtrs(nargc,(char const *const *const) argvNew,Nref,update_argv,2);
		if (DB.PGlobal <= 1) {
			test_print_warning("Please increase PGlobal above 1 in the ctrl file");
			printf("Ctrl file: %s.ctrl\n",argvNew[1]);
		}

		run_test(&pass,"FullREFINE");
		if (strstr(EName,"TRI"))
			sprintf(PrintName,"update_h (%s%d FullREFINE):",EName,refType);
		else
			sprintf(PrintName,"         (%s%d FullREFINE):",EName,refType);
		test_print2(pass,PrintName);

		run_test(&pass,"FullCOARSE");
		test_print2(pass,"         (        FullCOARSE):");

		mark_VOLUMEs(HREFINE,Lmts[0]);
		mark_VOLUMEs(HCOARSE,Lmts[1]);

		mesh_update();

		mark_VOLUMEs(HREFINE,Lmts[2]);
		mark_VOLUMEs(HCOARSE,Lmts[3]);

		run_test(&pass,"Mixed");
		test_print2(pass,"         (        Mixed):");

		pass = 1;
		check_Jacobians(&pass);
		test_print2(pass,"         (        Jacobians):");

		code_cleanup();
	}
	DB.TETrefineType = TETrefineType;

	free(PrintName);
}

void test_integration_update_h(int nargc, char **argv)
{
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
	strcpy(argvNew[1],"test/update_h/Test_update_h_TRI");
	strcpy(EName,"TRI   ");
	Lmts[0]->XYZ[0] = -0.75; Lmts[0]->XYZ[1] = -0.75; Lmts[0]->type = 'd'; Lmts[0]->index = 0;
	Lmts[1]->XYZ[0] =  0.00; Lmts[1]->XYZ[1] =  0.00; Lmts[1]->type = 'd'; Lmts[1]->index = 1;
	Lmts[2]->XYZ[0] = -0.50; Lmts[2]->XYZ[1] = -0.50; Lmts[2]->type = 'd'; Lmts[2]->index = 0;
	Lmts[3]->XYZ[0] =  0.00; Lmts[3]->XYZ[1] =  0.00; Lmts[3]->type = 'o'; Lmts[3]->index = 1;
	test_update_h(nargc,argvNew,2,0,EName,Lmts);

	// **************************************************************************************************** //
	// QUADs
	strcpy(argvNew[1],"test/update_h/Test_update_h_QUAD");
	strcpy(EName,"QUAD  ");
	Lmts[0]->XYZ[0] = -0.50; Lmts[0]->XYZ[1] = -0.50; Lmts[0]->type = 'a'; Lmts[0]->index = 0;
	Lmts[1]->XYZ[0] =  0.00; Lmts[1]->XYZ[1] =  0.00; Lmts[1]->type = 'a'; Lmts[1]->index = 1;
	Lmts[2]->XYZ[0] = -0.25; Lmts[2]->XYZ[1] = -0.25; Lmts[2]->type = 'a'; Lmts[2]->index = 0;
	Lmts[3]->XYZ[0] =  0.00; Lmts[3]->XYZ[1] =  0.00; Lmts[3]->type = 'o'; Lmts[3]->index = 1;
	test_update_h(nargc,argvNew,3,0,EName,Lmts);

	// **************************************************************************************************** //
	// TETs
	strcpy(argvNew[1],"test/update_h/Test_update_h_TET");
	strcpy(EName,"TET   ");
	Lmts[0]->XYZ[0] =  1.00; Lmts[0]->XYZ[1] =  1.00; Lmts[0]->XYZ[2] =  0.00; Lmts[0]->type = 'd'; Lmts[0]->index = 1;
	Lmts[1]->XYZ[0] = -1.00; Lmts[1]->XYZ[1] = -1.00; Lmts[1]->XYZ[2] =  1.00; Lmts[1]->type = 'd'; Lmts[1]->index = 0;
	Lmts[2]->XYZ[0] =  1.00; Lmts[2]->XYZ[1] =  1.00; Lmts[2]->XYZ[2] = -0.50; Lmts[2]->type = 'd'; Lmts[2]->index = 1;
	Lmts[3]->XYZ[0] =  1.00; Lmts[3]->XYZ[1] =  1.00; Lmts[3]->XYZ[2] = -1.00; Lmts[3]->type = 'd'; Lmts[3]->index = 0;
	test_update_h(nargc,argvNew,2,0,EName,Lmts);

	// **************************************************************************************************** //
	// HEXs
	strcpy(argvNew[1],"test/update_h/Test_update_h_HEX");
	strcpy(EName,"HEX   ");
	Lmts[0]->XYZ[0] = -0.50; Lmts[0]->XYZ[1] = -0.50; Lmts[0]->XYZ[2] = -0.50; Lmts[0]->type = 'a'; Lmts[0]->index = 0;
	Lmts[1]->XYZ[0] =  0.00; Lmts[1]->XYZ[1] =  0.00; Lmts[1]->XYZ[2] =  0.00; Lmts[1]->type = 'a'; Lmts[1]->index = 1;
	Lmts[2]->XYZ[0] = -0.25; Lmts[2]->XYZ[1] = -0.25; Lmts[2]->XYZ[2] = -0.25; Lmts[2]->type = 'a'; Lmts[2]->index = 0;
	Lmts[3]->XYZ[0] =  0.00; Lmts[3]->XYZ[1] =  0.00; Lmts[3]->XYZ[2] =  0.00; Lmts[3]->type = 'o'; Lmts[3]->index = 1;
	test_update_h(nargc,argvNew,3,0,EName,Lmts);

	// **************************************************************************************************** //
	// WEDGEs
	strcpy(argvNew[1],"test/update_h/Test_update_h_WEDGE");
	strcpy(EName,"WEDGE ");
	Lmts[0]->XYZ[0] = -0.50; Lmts[0]->XYZ[1] = -0.50; Lmts[0]->XYZ[2] = -0.50; Lmts[0]->type = 'a'; Lmts[0]->index = 0;
	Lmts[1]->XYZ[0] =  0.00; Lmts[1]->XYZ[1] =  0.00; Lmts[1]->XYZ[2] =  0.00; Lmts[1]->type = 'a'; Lmts[1]->index = 1;
	Lmts[2]->XYZ[0] = -0.25; Lmts[2]->XYZ[1] = -0.25; Lmts[2]->XYZ[2] = -0.25; Lmts[2]->type = 'a'; Lmts[2]->index = 0;
	Lmts[3]->XYZ[0] =  0.00; Lmts[3]->XYZ[1] =  0.00; Lmts[3]->XYZ[2] =  0.00; Lmts[3]->type = 'o'; Lmts[3]->index = 1;
	test_update_h(nargc,argvNew,3,0,EName,Lmts);

	// **************************************************************************************************** //
	// PYRs
	strcpy(argvNew[1],"test/update_h/Test_update_h_PYR");
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
					EXIT_UNSUPPORTED;
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
					EXIT_UNSUPPORTED;
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
				EXIT_UNSUPPORTED;
				break;
			}
			break;
		default:
			EXIT_UNSUPPORTED;
			break;
		}

		if (update) {
			VOLUME->Vadapt = 1;
			VOLUME->adapt_type = adapt_type;
		}
	}
}


struct S_OPERATORS {
	unsigned int NfnS, NvnGs, *nOrdLR, *nOrdRL;
	double       **I_vGs_fS;
};

static void init_ops(struct S_OPERATORS *OPS, const struct S_VOLUME *VOLUME, const struct S_FACE *FACE,
                     const unsigned int IndFType)
{
	unsigned int PF, VType, IndOrdLR, IndOrdRL;

	struct S_ELEMENT *ELEMENT, *ELEMENT_FACE;

//	PV    = VOLUME->P;
	VType = VOLUME->type;

	PF = FACE->P;
	IndOrdLR = FACE->IndOrdLR;
	IndOrdRL = FACE->IndOrdRL;

	ELEMENT = get_ELEMENT_type(VType);
	ELEMENT_FACE = get_ELEMENT_FACE(VType,IndFType);

	OPS->NvnGs = ELEMENT->NvnGs[1];
	OPS->NfnS  = ELEMENT_FACE->NvnS[PF];

	OPS->I_vGs_fS = ELEMENT->I_vGs_fS[1][PF];

	OPS->nOrdLR = ELEMENT_FACE->nOrd_fS[PF][IndOrdLR];
	OPS->nOrdRL = ELEMENT_FACE->nOrd_fS[PF][IndOrdRL];
}

static void check_correspondence(unsigned int *pass)
{
	// Initialize DB Parameters
	unsigned int d = DB.d;

	// Standard datatypes
	unsigned int Vf, IndFType, NfnS, *nOrdLR, *nOrdRL, vhIn, vhOut,
	             dim, n, Indd, BC, FACE_is_internal;
	double       *XYZ_fSIn, *XYZ_fSOut, *XYZ_fSLR, *XYZ_fSRL;

	struct S_OPERATORS *OPS;
	struct S_VOLUME    *VOLUME, *VOLUMEc;
	struct S_FACE     *FACE;

	OPS = malloc(sizeof *OPS); // free

	*pass = 1;
	for (FACE = DB.FACE; FACE; FACE = FACE->next) {
		VOLUME = FACE->VL;
		Vf     = FACE->VfL;

		IndFType = get_IndFType(VOLUME->Eclass,Vf/NFREFMAX);
		init_ops(OPS,VOLUME,FACE,IndFType);

		NfnS = OPS->NfnS;

		nOrdLR = OPS->nOrdLR;
		nOrdRL = OPS->nOrdRL;

		XYZ_fSIn = mm_Alloc_d(CBCM,CBT,CBNT,NfnS,d,OPS->NvnGs,1.0,OPS->I_vGs_fS[Vf],VOLUME->XYZ_vV); // free

		VOLUME = FACE->VR;
		Vf     = FACE->VfR;

		IndFType = get_IndFType(VOLUME->Eclass,Vf/NFREFMAX);
		init_ops(OPS,VOLUME,FACE,IndFType);

		XYZ_fSOut = mm_Alloc_d(CBCM,CBT,CBNT,NfnS,d,OPS->NvnGs,1.0,OPS->I_vGs_fS[Vf],VOLUME->XYZ_vV); // free

		XYZ_fSLR = malloc(NfnS*d * sizeof *XYZ_fSLR); // free
		XYZ_fSRL = malloc(NfnS*d * sizeof *XYZ_fSRL); // free

		for (dim = 0; dim < d; dim++) {
			Indd = dim*NfnS;
			for (n = 0; n < NfnS; n++) {
				XYZ_fSLR[Indd+n] = XYZ_fSIn[Indd+nOrdLR[n]];
				XYZ_fSRL[Indd+n] = XYZ_fSOut[Indd+nOrdRL[n]];
			}
		}

		BC = FACE->BC;
//		FACE_is_internal = (BC == 0 || (BC % BC_STEP_SC > 50));
		FACE_is_internal = (BC == 0); // Not including periodic faces

		if (FACE_is_internal && (array_norm_diff_d(NfnS*d,XYZ_fSIn,XYZ_fSRL,"Inf")  > 10*EPS ||
		                          array_norm_diff_d(NfnS*d,XYZ_fSLR,XYZ_fSOut,"Inf") > 10*EPS)) {
				*pass = 0;
				printf("Problem in check_correspondence\n");

				vhIn = 0;
				for (VOLUMEc = FACE->VL->parent->child0; VOLUMEc != FACE->VL; VOLUMEc = VOLUMEc->next)
					vhIn++;

				vhOut = 0;
				for (VOLUMEc = FACE->VR->parent->child0; VOLUMEc != FACE->VR ; VOLUMEc = VOLUMEc->next)
					vhOut++;

				printf("%d %d %d %d %d\n",FACE->indexg,FACE->IndOrdLR,FACE->IndOrdRL,vhIn,vhOut);
				printf("%d %d %d %d\n",FACE->VL->type,FACE->VL->indexg,FACE->VfL,FACE->VL->level);
				printf("%d %d %d %d\n",FACE->VR->type,FACE->VR->indexg,FACE->VfR,FACE->VR->level);
				printf("Errors: %e %e\n\n",array_norm_diff_d(NfnS*d,XYZ_fSIn,XYZ_fSRL,"Inf"),
		                                   array_norm_diff_d(NfnS*d,XYZ_fSLR,XYZ_fSOut,"Inf"));
				array_print_d(NfnS,d,XYZ_fSIn,'C');
				array_print_d(NfnS,d,XYZ_fSRL,'C');
				array_print_d(NfnS,d,XYZ_fSOut,'C');
				array_print_d(NfnS,d,XYZ_fSLR,'C');
				EXIT_MSG;
				break;
		}

		free(XYZ_fSIn);
		free(XYZ_fSOut);
		free(XYZ_fSLR);
		free(XYZ_fSRL);
	}
	free(OPS);
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
