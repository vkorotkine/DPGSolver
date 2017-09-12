// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/// \file

#include "cubature.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "jacobi.h"
#include "gsl/gsl_errno.h"

#include "macros.h"
#include "definitions_cubature.h"
#include "definitions_math.h"
#include "definitions_tol.h"
#include "definitions_alloc.h"

#include "matrix.h"
#include "vector.h"

#include "math_functions.h"

// Static function declarations ************************************************************************************* //

/** \brief Constructor for a \ref Cubature container of simplex (tri) type.
 *  \return Standard. */
static const struct const_Cubature* constructor_const_Cubature_tri
	(const int p,        ///< Defined for \ref constructor_const_Cubature_si.
	 const int node_type ///< Defined for \ref constructor_const_Cubature_si.
	);

/** \brief Constructor for a \ref Cubature container of simplex (tet) type.
 *  \return Standard. */
static const struct const_Cubature* constructor_const_Cubature_tet
	(const int p,        ///< Defined for \ref constructor_const_Cubature_si.
	 const int node_type ///< Defined for \ref constructor_const_Cubature_si.
	);

// Interface functions ********************************************************************************************** //

const struct const_Cubature* constructor_const_Cubature_tp (const int d, const int p, const int node_type)
{
	if ((d <= 0) || (p < 0))
		EXIT_ERROR("Invalid dimension (%d) or order (%d).\n",d,p);

	bool has_weights = true;

	const int pp1 = p+1;

	double r[pp1],
	       w_1d[pp1];
	if (node_type == CUB_GL) {
		if (jac_zeros_gj(r,pp1,0.0,0.0) != GSL_SUCCESS)
			EXIT_ERROR("Problem computing nodes.\n");

		has_weights = true;
		if (jac_weights_gj(r,w_1d,pp1,0.0,0.0,NULL) != GSL_SUCCESS)
			EXIT_ERROR("Problem computing weights.\n");
	} else if (node_type == CUB_GLL) {
		if (p == 0)
			EXIT_ERROR("The order of the GLL nodes must be greater than 0.\n");

		if (jac_zeros_glj(r,pp1,0.0,0.0) != GSL_SUCCESS)
			EXIT_ERROR("Problem computing nodes.\n");

		has_weights = true;
		if (jac_weights_glj(r,w_1d,pp1,0.0,0.0,NULL) != GSL_SUCCESS)
			EXIT_ERROR("Problem computing weights.\n");
	} else if (node_type == CUB_EQ) {
		has_weights = false;
		for (int i = 0; i < pp1; ++i)
			r[i] = -1.0 + (2.0/p)*i;
	} else {
		EXIT_UNSUPPORTED;
	}

	const int ext_0 = pow(pp1,d);
	double *rst = malloc(ext_0*d * sizeof *rst); // keep

	int row = 0;
	for (int k = 0, k_max = GSL_MIN(GSL_MAX((d-2)*pp1,1),pp1); k < k_max; k++) {
	for (int j = 0, j_max = GSL_MIN(GSL_MAX((d-1)*pp1,1),pp1); j < j_max; j++) {
	for (int i = 0, i_max = GSL_MIN(GSL_MAX((d-0)*pp1,1),pp1); i < i_max; i++) {
		for (int dim = 0; dim < d; dim++) {
			switch (dim) {
			case 0:
				rst[dim*ext_0+row] = r[i];
				break;
			case 1:
				rst[dim*ext_0+row] = r[j];
				break;
			case 2:
				rst[dim*ext_0+row] = r[k];
				break;
			default:
				EXIT_UNSUPPORTED;
				break;
			}
		}
		row++;
	}}}

	double* w = NULL;
	if (has_weights) {
		w = malloc(ext_0 * sizeof *w); // keep

		int row = 0;
		for (int k = 0, k_max = GSL_MIN(GSL_MAX((d-2)*pp1,1),pp1); k < k_max; k++) {
		for (int j = 0, j_max = GSL_MIN(GSL_MAX((d-1)*pp1,1),pp1); j < j_max; j++) {
		for (int i = 0, i_max = GSL_MIN(GSL_MAX((d-0)*pp1,1),pp1); i < i_max; i++) {
			w[row] = w_1d[i];
			if (d == 2) w[row] *= w_1d[j];
			if (d == 3) w[row] *= w_1d[j]*w_1d[k];
			row++;
		}}}
	}

	struct Cubature* cubature = calloc(1,sizeof *cubature); // returned;

	cubature->p         = p;
	cubature->node_type = node_type;

	cubature->rst = constructor_move_Matrix_d_d('C',ext_0,d,true,rst); // keep

	cubature->has_weights = has_weights;
	if (has_weights)
		cubature->w = constructor_move_Vector_d_d(ext_0,true,w); // keep
	else
		cubature->w = NULL;

	return (const struct const_Cubature*) cubature;
}

const struct const_Cubature* constructor_const_Cubature_si (const int d, const int p, const int node_type)
{
	if ((d <= 1) || (p < 0))
		EXIT_ERROR("Invalid dimension (%d) or order (%d).\n",d,p);

	if (d == 2)
		return constructor_const_Cubature_tri(p,node_type);
	else if (d == 3)
		return constructor_const_Cubature_tet(p,node_type);
	else
		EXIT_UNSUPPORTED;
}

const struct const_Cubature* constructor_const_Cubature_pyr (const int d, const int p, const int node_type)
{
	if (d != 3)
		EXIT_UNSUPPORTED;

	const unsigned int P = p;
	bool has_weights = true;

	static const unsigned int perms_QUAD[] = { 0,1,2,3, 3,0,1,2, 2,3,0,1, 1,2,3,0 };
	static const double rst_V[] = { -1.0,        1.0,        1.0,       -1.0,       0.0,
	                                -1.0,       -1.0,        1.0,        1.0,       0.0,
	                                -SQRT2/5.0, -SQRT2/5.0, -SQRT2/5.0, -SQRT2/5.0, 4.0/5.0*SQRT2};

	unsigned int symms41[2] = {0, 0}, NQUADsymms = 0;
	double BCoords_tmp[5];
	unsigned int i, iMax, j, k, QUADsymm,
	             IndB, IndBC, IndGroup, Ind1, GroupCount, Nc, N1,
	             PMax, NnOut, Ngroups, Nsymms, Nperms;
	unsigned int *symms_Nperms, *symms_count;
	char         *StringRead, *strings, *stringe, *CubFile, *Pc;
	double       *rstOut, *wOut, *BCoords, *BCoords_complete, *w_read;

	FILE *fID;

	// Silence compiler warnings
	PMax = 0; Nsymms = 0; Ngroups = 0;
	w_read = malloc(0 * sizeof *w_read); // silence
	wOut = NULL;

	const char* NodeType = NULL;
	if (node_type == CUB_GL) {
		has_weights = false;
		NodeType = "GL";
		PMax = 6;
	} else if (node_type == CUB_GLL) {
		has_weights = false;
		NodeType = "GLL";
		PMax = 6;
	} else if (node_type == CUB_GLW) {
		NodeType = "GLW";
		PMax = 6;
	} else if (node_type == CUB_GLLW) {
		NodeType = "GLLW";
		PMax = 6;
	} else if (node_type == CUB_GJW) {
		NodeType = "GJW";
		PMax = 10;
	} else if (node_type == CUB_WV) {
		NodeType = "WV";
		PMax = 10;
	} else if (node_type == CUB_WVHToP) {
		NodeType = "WV";
		PMax = 11;
	} else {
		EXIT_UNSUPPORTED;
	}

	if (P > PMax)
		printf("Error: %s PYR nodes of order %d are not available.\n",NodeType,P), EXIT_MSG;

	Nc = 5;

	CubFile = malloc(STRLEN_MAX * sizeof *CubFile); // free
	Pc      = malloc(STRLEN_MIN * sizeof *Pc);      // free
	sprintf(Pc,"%d",P);

	strcpy(CubFile,"../cubature/pyr/");
	strcat(CubFile,NodeType);
	strcat(CubFile,Pc);
	strcat(CubFile,".txt");
	free(Pc);

	if ((fID = fopen(CubFile,"r")) == NULL)
		printf("Error: Cubature file %s not found.\n",CubFile), exit(1);
	free(CubFile);

	StringRead = malloc(STRLEN_MAX * sizeof *StringRead); // free

	// Read header
	if (fscanf(fID,"%[^\n]\n",StringRead) == 1) {
		strings = StringRead;
		Ngroups = strtol(strings,&stringe,10); strings = stringe;
		Nsymms  = strtol(strings,&stringe,10); strings = stringe;
	}

	symms_count  = malloc(Nsymms * sizeof *symms_count);  // free
	symms_Nperms = malloc(Nsymms * sizeof *symms_Nperms); // free

	for (i = 0; i < Nsymms; i++) {
		if (fscanf(fID,"%[^\n]\n",StringRead) == 1) {
			strings = StringRead;
			symms_count[i]  = strtol(strings,&stringe,10); strings = stringe;
			symms_Nperms[i] = strtol(strings,&stringe,10); strings = stringe;
		}
	}

	for (iMax = 1; iMax--; ) {
		if (fscanf(fID,"%[^\n]\n",StringRead) == 1) {
			; // skip 1 line(s)
		}
	}

	BCoords = malloc(Ngroups*Nc * sizeof * BCoords); // free

	for (i = 0; i < Ngroups; i++) {
		if (fscanf(fID,"%[^\n]\n",StringRead) == 1) {
			strings = StringRead;
			for (j = 0; j < Nc; j++) {
				BCoords[i*Nc+j] = strtod(strings,&stringe); strings = stringe;
			}
		}
	}

	// Read weights
	if (has_weights) {
		for (iMax = 1; iMax--; ) {
			if (fscanf(fID,"%[^\n]\n",StringRead) == 1) {
				; // skip 1 line(s)
			}
		}

		free(w_read);
		w_read = malloc(Ngroups * sizeof *w_read); // free

		for (i = 0; i < Ngroups; i++) {
			if (fscanf(fID,"%[^\n]\n",StringRead) == 1) {
				strings = StringRead;
				w_read[i] = strtod(strings,&stringe); strings = stringe;
			}
		}
	}

	fclose(fID);
	free(StringRead);

	// Convert barycentric coordinates to rst nodal coordinates
	// Also set weights in array of size NnOut x 1 if applicable

	// Find number of 1 symmetries and total number of nodes
	for (i = 0, N1 = 0, NnOut = 0; i < Nsymms; i++) {
		if (i >= Nsymms-1)
			N1 += symms_count[i];
		NnOut += symms_count[i]*symms_Nperms[i];
	}

	if (has_weights)
		wOut = malloc(NnOut * sizeof *wOut); // free/keep (conditional return_w)

	BCoords_complete = malloc(NnOut*Nc * sizeof *BCoords_complete); // free

	IndB = 0; IndBC = 0;
	GroupCount = 0;

	IndGroup = 0;
	for (i = 0; i < Nsymms; i++) {
		if (symms_count[i] == 0)
			IndGroup++;
		else
			break;
	}

	Ind1 = 0;
	for (i = 0; i < Ngroups; i++) {
		GroupCount++;

		Nperms = symms_Nperms[IndGroup];
		// Count 4/1 symmetries and establish QUAD symmetries to loop over
		if (Nperms == 8) {
			symms41[0] += 2;
			NQUADsymms = 1;
		} else if (Nperms == 4) {
			symms41[0] += 1;
			NQUADsymms = 0;
		} else if (Nperms == 1) {
			symms41[1] += 1;
			NQUADsymms = 0;
		}

		if (Nperms == 1) {
			for (j = 0; j < Nc; j++)
				BCoords_tmp[j] = BCoords[IndB*Nc+j];

			for (k = 0; k < Nc; k++)
				BCoords_complete[(NnOut-N1+Ind1)*Nc+k] = BCoords_tmp[k];

			if (has_weights)
				wOut[NnOut-N1+Ind1] = w_read[IndB];

			Ind1++;
		} else {
			for (QUADsymm = 0; QUADsymm <= NQUADsymms; QUADsymm++) {
				if (QUADsymm == 0) {
					for (j = 0; j < Nc; j++)
						BCoords_tmp[j] = BCoords[IndB*Nc+j];
				} else {
					BCoords_tmp[0] = BCoords[IndB*Nc+1];
					BCoords_tmp[1] = BCoords[IndB*Nc+0];
					BCoords_tmp[2] = BCoords[IndB*Nc+3];
					BCoords_tmp[3] = BCoords[IndB*Nc+2];
					BCoords_tmp[4] = BCoords[IndB*Nc+4];
				}

				for (j = 0; j < 4; j++) {
					if (has_weights)
						wOut[IndBC] = w_read[IndB];

					for (k = 0; k < Nc-1; k++)
						BCoords_complete[IndBC*Nc+k] = BCoords_tmp[perms_QUAD[j*(Nc-1)+k]];

					k = Nc-1;
					BCoords_complete[IndBC*Nc+k] = BCoords_tmp[4];

					IndBC++;
				}
			}
		}
		IndB++;

		if (symms_count[IndGroup] == GroupCount) {
			GroupCount = 0;
			IndGroup += 1;
		}
	}

	free(symms_count);
	free(symms_Nperms);
	free(BCoords);
	free(w_read);

	struct Matrix_d* BCoords_M = constructor_move_Matrix_d_d('R',NnOut,Nc,false,BCoords_complete); // destructed
	struct Matrix_d* rst_V_M   = constructor_move_Matrix_d_d('C',Nc,d,false,(double*)rst_V); // destructed
	struct Matrix_d* rstOut_M  =
		constructor_mm_Matrix_d('N','N',1.0,0.0,(struct const_Matrix_d*)BCoords_M,
		                        (struct const_Matrix_d*)rst_V_M,'C'); // destructed
	free(BCoords_complete);

	rstOut = rstOut_M->data;
	rstOut_M->owns_data = false;

	destructor_Matrix_d(BCoords_M);
	destructor_Matrix_d(rst_V_M);
	destructor_Matrix_d(rstOut_M);


	struct Cubature* cubature = calloc(1,sizeof *cubature); // returned;

	cubature->p         = p;
	cubature->node_type = node_type;

	const ptrdiff_t ext_0 = NnOut;
	cubature->rst = constructor_move_Matrix_d_d('C',ext_0,d,true,rstOut); // keep

	cubature->has_weights = has_weights;
	if (has_weights)
		cubature->w = constructor_move_Vector_d_d(ext_0,true,wOut); // keep
	else
		cubature->w = NULL;

	return (const struct const_Cubature*) cubature;
}

void destructor_Cubature (struct Cubature* cub)
{
	destructor_Matrix_d(cub->rst);
	if (cub->has_weights)
		destructor_Vector_d(cub->w);
	free(cub);
}

void destructor_const_Cubature (const struct const_Cubature*const cub)
{
	destructor_Cubature((struct Cubature*)cub);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static const struct const_Cubature* constructor_const_Cubature_tri (const int p, const int node_type)
{
	const unsigned int d = 2;
	const unsigned int P = p;
	bool has_weights = true;

	static const unsigned int perms[] = { 0,1,2, 2,0,1, 1,2,0, 0,2,1, 1,0,2, 2,1,0 };
	static const double rst_V[] = { -1.0,        1.0,       0.0,
	                                -1.0/SQRT3, -1.0/SQRT3, 2.0/SQRT3};

	// Standard datatypes
	unsigned int symms31[2] = {0, 0};
	unsigned int i, iMax, j, jMax, k,
	             IndB, IndBC, IndGroup, GroupCount, Nc,
	             PMax, NnOut, Ngroups, Nsymms;
	unsigned int  *symms_Nperms, *symms_count;
	char         *StringRead, *strings, *stringe, *CubFile, *Pc;
	double       *rstOut, *wOut, *BCoords, *BCoords_complete, *w_read;

	// Silence compiler warnings
	PMax = 0; Nsymms = 0; Ngroups = 0;
	w_read = malloc(0 * sizeof *w_read); // silence (ToBeModified: change this to NULL and remove free below)
	wOut = NULL;

	const char* NodeType = NULL;
	if (node_type == CUB_AO) {
		has_weights = false;
		NodeType = "AO";
		PMax = 15;
	} else if (node_type == CUB_WSH) {
		NodeType = "WSH";
		PMax = 8;
	} else if (node_type == CUB_WV) {
		NodeType = "WV";
		PMax = 20;
	} else if (node_type == CUB_EQ) {
		has_weights = false;
		NodeType = "EQ";
		PMax = 8;
	} else {
		EXIT_UNSUPPORTED;
	}

	if (P > PMax)
		EXIT_ERROR("%s TRI nodes of order %d are not available.\n",NodeType,P);

	FILE *fID;
	Nc = 3;

	CubFile = malloc(STRLEN_MAX * sizeof *CubFile); // free
	Pc      = malloc(STRLEN_MIN * sizeof *Pc);      // free
	sprintf(Pc,"%d",P);

	strcpy(CubFile,"../cubature/tri/");
	strcat(CubFile,NodeType);
	strcat(CubFile,Pc);
	strcat(CubFile,".txt");
	free(Pc);

	if ((fID = fopen(CubFile,"r")) == NULL)
		printf("Error: Cubature file %s not found.\n",CubFile), exit(1);
	free(CubFile);

	StringRead = malloc(STRLEN_MAX * sizeof *StringRead); // free

	// Read header
	if (fscanf(fID,"%[^\n]\n",StringRead) == 1) {
		strings = StringRead;
		Ngroups = strtol(strings,&stringe,10); strings = stringe;
		Nsymms  = strtol(strings,&stringe,10); strings = stringe;
	}

	symms_count  = malloc(Nsymms * sizeof *symms_count);  // free
	symms_Nperms = malloc(Nsymms * sizeof *symms_Nperms); // free

	for (i = 0; i < Nsymms; i++) {
		if (fscanf(fID,"%[^\n]\n",StringRead) == 1) {
			strings = StringRead;
			symms_count[i]  = strtol(strings,&stringe,10); strings = stringe;
			symms_Nperms[i] = strtol(strings,&stringe,10); strings = stringe;
		}
	}

	for (iMax = 1; iMax--; ) {
		if (fscanf(fID,"%[^\n]\n",StringRead) == 1) {
			; // skip 1 line(s)
		}
	}

	BCoords = malloc(Ngroups*Nc * sizeof * BCoords); // free

	for (i = 0; i < Ngroups; i++) {
		if (fscanf(fID,"%[^\n]\n",StringRead) == 1) {
			strings = StringRead;
			for (j = 0; j < Nc; j++) {
				BCoords[i*Nc+j] = strtod(strings,&stringe); strings = stringe;
			}
		}
	}

	// Read weights
	if (has_weights) {
		for (iMax = 1; iMax--; ) {
			if (fscanf(fID,"%[^\n]\n",StringRead) == 1) {
				; // skip 1 line(s)
			}
		}

		free(w_read);
		w_read = malloc(Ngroups * sizeof *w_read); // free

		for (i = 0; i < Ngroups; i++) {
			if (fscanf(fID,"%[^\n]\n",StringRead) == 1) {
				strings = StringRead;
				w_read[i] = strtod(strings,&stringe); strings = stringe;
			}
		}
	}

	fclose(fID);
	free(StringRead);

	// Convert barycentric coordinates to rst nodal coordinates
	// Also set weights in array of size NnOut x 1 if applicable
	NnOut = 0;
	for (i = 0; i < Nsymms; i++)
		NnOut += symms_count[i]*symms_Nperms[i];

	if (has_weights)
		wOut = malloc(NnOut * sizeof *wOut); // free/keep (conditional return_w)

	BCoords_complete = malloc(NnOut*Nc * sizeof *BCoords_complete); // free

	IndB = 0; IndBC = 0;
	GroupCount = 0;

	IndGroup = 0;
	for (i = 0; i < Nsymms; i++) {
		if (symms_count[i] == 0)
			IndGroup++;
		else
			break;
	}

	for (i = 0; i < Ngroups; i++) {
		GroupCount++;

		jMax = symms_Nperms[IndGroup];
		// Count 3/1 symmetries
		if      (jMax == 6) symms31[0] += 2;
		else if (jMax == 3) symms31[0] += 1;
		else if (jMax == 1) symms31[1] += 1;

		for (j = 0; j < jMax; j++) {
			if (has_weights)
				wOut[IndBC] = w_read[IndB];

			for (k = 0; k < Nc; k++) {
				BCoords_complete[IndBC*Nc+k] = BCoords[IndB*Nc+perms[j*Nc+k]];
			}
			IndBC++;
		}
		IndB++;

		if (symms_count[IndGroup] == GroupCount) {
			GroupCount = 0;
			IndGroup += 1;
		}
	}
	free(symms_count);
	free(symms_Nperms);
	free(BCoords);
	free(w_read);

	struct Matrix_d* BCoords_M = constructor_move_Matrix_d_d('R',NnOut,Nc,false,BCoords_complete); // destructed
	struct Matrix_d* rst_V_M   = constructor_move_Matrix_d_d('C',Nc,d,false,(double*)rst_V); // destructed
	struct Matrix_d* rstOut_M  =
		constructor_mm_Matrix_d('N','N',1.0,0.0,(struct const_Matrix_d*)BCoords_M,
		                        (struct const_Matrix_d*)rst_V_M,'C'); // destructed
	free(BCoords_complete);

	rstOut = rstOut_M->data;
	rstOut_M->owns_data = false;

	destructor_Matrix_d(BCoords_M);
	destructor_Matrix_d(rst_V_M);
	destructor_Matrix_d(rstOut_M);


	struct Cubature* cubature = calloc(1,sizeof *cubature); // returned;

	cubature->p         = p;
	cubature->node_type = node_type;

	const ptrdiff_t ext_0 = NnOut;
	cubature->rst = constructor_move_Matrix_d_d('C',ext_0,d,true,rstOut); // keep

	cubature->has_weights = has_weights;
	if (has_weights)
		cubature->w = constructor_move_Vector_d_d(ext_0,true,wOut); // keep
	else
		cubature->w = NULL;

	return (const struct const_Cubature*) cubature;
}

static const struct const_Cubature* constructor_const_Cubature_tet (const int p, const int node_type)
{
	const unsigned int d = 3;
	const unsigned int P = p;
	bool has_weights = true;

	static const unsigned int perms_TRI[] = { 0,1,2, 2,0,1, 1,2,0, 0,2,1, 1,0,2, 2,1,0 };
	static const unsigned int perms_TET[] = { 0,1,2,3, 3,0,1,2, 2,3,0,1, 1,2,3,0 };
	static const double rst_V[] = { -1.0,        1.0,        0.0,       0.0,
	                                -1.0/SQRT3, -1.0/SQRT3,  2.0/SQRT3, 0.0,
	                                -1.0/SQRT6, -1.0/SQRT6, -1.0/SQRT6, 3.0/SQRT6};

	unsigned int symms31[2] = {0, 0}, NTRIsymms[3], TETperms[4];
	double BCoords_tmp[4];

	unsigned int i, iMax, j, jMax, k, kMax, l, lMax, TRIsymm,
	             IndB, IndBC, IndGroup, Ind1, Indperm, GroupCount, Nc, N1,
	             PMax, NnOut, Ngroups, Nsymms, Nperms;
	unsigned int *symms_Nperms, *symms_count;
	char         *StringRead, *strings, *stringe, *CubFile, *Pc;
	double       *rstOut, *wOut, *BCoords, *BCoords_complete, *w_read;

	FILE *fID;

	// Silence compiler warnings
	PMax = 0; Nsymms = 0; Ngroups = 0;
	w_read = wOut = NULL;

	const char* NodeType = NULL;
	if (node_type == CUB_AO) {
		has_weights = false;
		NodeType = "AO";
		PMax = 15;
	} else if (node_type == CUB_WSH) {
		NodeType = "WSH";
		PMax = 6;
	} else if (node_type == CUB_WV) {
		NodeType = "WV";
		PMax = 10;
	} else {
		EXIT_UNSUPPORTED;
	}

	if (P > PMax)
		printf("Error: %s TET nodes of order %d are not available.\n",NodeType,P), exit(1);

	Nc = 4;

	CubFile = malloc(STRLEN_MAX * sizeof *CubFile); // free
	Pc      = malloc(STRLEN_MIN * sizeof *Pc);      // free
	sprintf(Pc,"%d",P);

	strcpy(CubFile,"../cubature/tet/");
	strcat(CubFile,NodeType);
	strcat(CubFile,Pc);
	strcat(CubFile,".txt");
	free(Pc);

	if ((fID = fopen(CubFile,"r")) == NULL)
		printf("Error: Cubature file %s not found.\n",CubFile), exit(1);
	free(CubFile);

	StringRead = malloc(STRLEN_MAX * sizeof *StringRead); // free

	// Read header
	if (fscanf(fID,"%[^\n]\n",StringRead) == 1) {
		strings = StringRead;
		Ngroups = strtol(strings,&stringe,10); strings = stringe;
		Nsymms  = strtol(strings,&stringe,10); strings = stringe;
	}

	symms_count  = malloc(Nsymms * sizeof *symms_count);  // free
	symms_Nperms = malloc(Nsymms * sizeof *symms_Nperms); // free

	for (i = 0; i < Nsymms; i++) {
		if (fscanf(fID,"%[^\n]\n",StringRead) == 1) {
			strings = StringRead;
			symms_count[i]  = strtol(strings,&stringe,10); strings = stringe;
			symms_Nperms[i] = strtol(strings,&stringe,10); strings = stringe;
		}
	}

	for (iMax = 1; iMax--; ) {
		if (fscanf(fID,"%[^\n]\n",StringRead) == 1) {
			; // skip 1 line(s)
		}
	}

	BCoords = malloc(Ngroups*Nc * sizeof * BCoords); // free

	for (i = 0; i < Ngroups; i++) {
		if (fscanf(fID,"%[^\n]\n",StringRead) == 1) {
			strings = StringRead;
			for (j = 0; j < Nc; j++) {
				BCoords[i*Nc+j] = strtod(strings,&stringe); strings = stringe;
			}
		}
	}

	// Read weights
	if (has_weights) {
		for (iMax = 1; iMax--; ) {
			if (fscanf(fID,"%[^\n]\n",StringRead) == 1) {
				; // skip 1 line(s)
			}
		}

		w_read = malloc(Ngroups * sizeof *w_read); // free

		for (i = 0; i < Ngroups; i++) {
			if (fscanf(fID,"%[^\n]\n",StringRead) == 1) {
				strings = StringRead;
				w_read[i] = strtod(strings,&stringe); strings = stringe;
			}
		}
	}

	fclose(fID);
	free(StringRead);

	// Convert barycentric coordinates to rst nodal coordinates
	// Also set weights in array of size NnOut x 1 if applicable

	// Find number of 1 symmetries and total number of nodes
	for (i = 0, N1 = 0, NnOut = 0; i < Nsymms; i++) {
		if (i >= Nsymms-2)
			N1 += symms_count[i];
		NnOut += symms_count[i]*symms_Nperms[i];
	}

	if (has_weights)
		wOut = malloc(NnOut * sizeof *wOut); // free/keep (conditional return_w)

	BCoords_complete = malloc(NnOut*Nc * sizeof *BCoords_complete); // free

	IndB = 0; IndBC = 0;
	GroupCount = 0;

	IndGroup = 0;
	for (i = 0; i < Nsymms; i++) {
		if (symms_count[i] == 0)
			IndGroup++;
		else
			break;
	}

	Ind1 = 0;
	for (i = 0; i < Ngroups; i++) {
		GroupCount++;

		Nperms = symms_Nperms[IndGroup];
		// Count 3/1 symmetries and establish TRI symmetries to loop over
		if (Nperms == 24) {
			symms31[0] += 8;

			NTRIsymms[0] = 0; NTRIsymms[1] = 0; NTRIsymms[2] = 4;
			TETperms[0]  = 0; TETperms[1]  = 1; TETperms[2]  = 2; TETperms[3]  = 3;
		} else if (Nperms == 12) {
			symms31[0] += 4;

			NTRIsymms[0] = 0; NTRIsymms[1] = 2; NTRIsymms[2] = 1;
			TETperms[0]  = 3; TETperms[1]  = 2; TETperms[2]  = 0; TETperms[3]  = 9;
		} else if (Nperms == 6) {
			symms31[0] += 2;

			NTRIsymms[0] = 0; NTRIsymms[1] = 2; NTRIsymms[2] = 0;
			TETperms[0]  = 0; TETperms[1]  = 3; TETperms[2]  = 9; TETperms[3]  = 9;
		} else if (Nperms == 4) {
			symms31[0] += 1, symms31[1] += 1;

			NTRIsymms[0] = 1; NTRIsymms[1] = 1; NTRIsymms[2] = 0;
			TETperms[0]  = 3; TETperms[1]  = 0; TETperms[2]  = 9; TETperms[3]  = 9;
		} else if (Nperms == 1) {
			symms31[1] += 1;

			NTRIsymms[0] = 1; NTRIsymms[1] = 0; NTRIsymms[2] = 0;
			TETperms[0]  = 0; TETperms[1]  = 9; TETperms[2]  = 9; TETperms[3]  = 9;
		}

		Indperm = 0;
		if (NTRIsymms[0]) {
			for (j = 0; j < Nc; j++)
				BCoords_tmp[j] = BCoords[IndB*Nc+perms_TET[TETperms[Indperm]*Nc+j]];
			Indperm++;

			for (k = 0; k < Nc; k++)
				BCoords_complete[(NnOut-N1+Ind1)*Nc+k] = BCoords_tmp[k];

			if (has_weights)
				wOut[NnOut-N1+Ind1] = w_read[IndB];

			Ind1++;
		}

		for (TRIsymm = 1; TRIsymm <= 2; TRIsymm++) {
			for (l = 0, lMax = NTRIsymms[TRIsymm]; l < lMax; l++) {

				for (j = 0; j < Nc; j++)
					BCoords_tmp[j] = BCoords[IndB*Nc+perms_TET[TETperms[Indperm]*Nc+j]];
				Indperm++;

				jMax = 0;
				for (k = 1, kMax = TRIsymm+1; k <= kMax; k++)
					jMax += k;

				for (j = 0; j < jMax; j++) {
					if (has_weights)
						wOut[IndBC] = w_read[IndB];

					for (k = 0; k < Nc-1; k++)
						BCoords_complete[IndBC*Nc+k] = BCoords_tmp[perms_TRI[j*(Nc-1)+k]];

					k = Nc-1;
					BCoords_complete[IndBC*Nc+k] = BCoords_tmp[3];

					IndBC++;
				}
			}
		}
		IndB++;

		if (symms_count[IndGroup] == GroupCount) {
			GroupCount = 0;
			IndGroup++;
			while (IndGroup < Nsymms && symms_count[IndGroup] == 0)
				IndGroup++;
		}
	}

	free(symms_count);
	free(symms_Nperms);
	free(BCoords);
	free(w_read);

	struct Matrix_d* BCoords_M = constructor_move_Matrix_d_d('R',NnOut,Nc,false,BCoords_complete); // destructed
	struct Matrix_d* rst_V_M   = constructor_move_Matrix_d_d('C',Nc,d,false,(double*)rst_V); // destructed
	struct Matrix_d* rstOut_M  =
		constructor_mm_Matrix_d('N','N',1.0,0.0,(struct const_Matrix_d*)BCoords_M,
		                        (struct const_Matrix_d*)rst_V_M,'C'); // destructed
	free(BCoords_complete);

	rstOut = rstOut_M->data;
	rstOut_M->owns_data = false;

	destructor_Matrix_d(BCoords_M);
	destructor_Matrix_d(rst_V_M);
	destructor_Matrix_d(rstOut_M);


	struct Cubature* cubature = calloc(1,sizeof *cubature); // returned;

	cubature->p         = p;
	cubature->node_type = node_type;

	const ptrdiff_t ext_0 = NnOut;
	cubature->rst = constructor_move_Matrix_d_d('C',ext_0,d,true,rstOut); // keep

	cubature->has_weights = has_weights;
	if (has_weights)
		cubature->w = constructor_move_Vector_d_d(ext_0,true,wOut); // keep
	else
		cubature->w = NULL;

	return (const struct const_Cubature*) cubature;
}
