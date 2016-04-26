#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "parameters.h"
#include "functions.h"

/*
 *	Purpose:
 *		Return nodes, weights and symmetries for triangular cubature depending on the nodetype.
 *
 *	Comments:
 *		Check that AO nodes correspond to those in pyfr after writing the converter script. (ToBeDeleted)
 *		The WS and WV nodes were determined based off of those from the pyfr code (pyfr/quadrules/tri) after being
 *		transfered to the equilateral reference TRI used in this code.
 *		The order of r, w is important for minimizing memory stride while computing the length 3 Discrete Fourier
 *		Transform (ToBeModified).
 *		Ordering convention:
 *			3-blocks of symmetric nodes going from farthest from center towards the center, followed by 1-block of
 *			center node if present.
 *		rst is stored in memory as r, s, then t (See cubature_TP for the motivation).
 *
 *		WS nodes have the following [order, cubature strength]:
 *			[0,1], [1,2], [2,4], [3,5], [4,7], [5,8], [6,10], [7,12], [8,14]
 *
 *		Input P for WV nodes is the cubature strength desired.
 *
 *	Notation:
 *		rst   : Nodes array of dimension Nn*d (column-major storage)
 *		w     : Weights array of dimension Nn*1
 *		symms : Symmetries array
 *		Nn    : (N)umber of (n)odes
 *		Ns    : (N)umber of (s)ymmetries
 *
 *	References:
 *		pyfr code : http://www.pyfr.org
 *
 *		AO : Hesthaven(2008)-Nodal_Discontinuous_Galerkin_Methods (ToBeModified: Likely use same conversion)
 *		WS : pyfr/quadrules/tri + conversion to barycentric coordinates (ToBeModified: See python script)
 *		WV : pyfr/quadrules/tri + conversion to barycentric coordinates (ToBeModified: See python script)
 */

void cubature_TRI(double **rst, double **w, unsigned int **symms, unsigned int *Nn, unsigned int *Ns,
                  const unsigned int return_w, const unsigned int P, const unsigned int d, const char *NodeType)
{
	// Standard datatypes
	static unsigned int perms[18] = { 0, 1, 2,
	                                  2, 0, 1,
	                                  1, 2, 0,
	                                  0, 2, 1,
	                                  1, 0, 2,
	                                  2, 1, 0};
	unsigned int symms31[2] = { 0, 0};
	double rst_c[6] = { -1.0,            1.0,           0.0,
	                    -1.0/sqrt(3.0), -1.0/sqrt(3.0), 2.0/sqrt(3.0)};
	unsigned int i, iMax, j, jMax, k,
	             IndB, IndBC, IndGroup, GroupCount, Nc,
	             PMax, NnOut, NsOut, Ngroups, Nsymms;
	unsigned int *symmsOut, *symms_Nperms, *symms_count;
	char         *StringRead, *strings, *stringe, *CubFile, *Pc;
	double       *rstOut, *wOut, *BCoords, *BCoords_complete, *w_read;

	FILE *fID;

	// Silence compiler warnings
	PMax = 0; Nsymms = 0; Ngroups = 0;
	w_read = malloc(0 * sizeof *w_read); // silence

	if (strstr(NodeType,"AO") != NULL) {
		if (return_w)
			printf("Error: Invalid value for return_w in cubature_TRI.\n"), exit(1);

		PMax = 16;
	} else if (strstr(NodeType,"WS") != NULL) {
		PMax = 8;
	} else if (strstr(NodeType,"WV") != NULL) {
		PMax = 20;
	}

	if (P > PMax)
		printf("Error: %s TRI nodes of order %d are not available.\n",NodeType,P), exit(1);

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
	if (return_w) {
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
//array_print_d(Ngroups,1,w_read,'R');
	}

	fclose(fID);
	free(StringRead);

	// Convert barycentric coordinates to rst nodal coordinates
	// Also set weights in array of size NnOut x 1 if applicable
	NnOut = 0;
	for (i = 0; i < Nsymms; i++)
		NnOut += symms_count[i]*symms_Nperms[i];

	if (return_w)
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
			if (return_w)
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

	rstOut = mm_Alloc_d(CblasColMajor,CblasTrans,CblasNoTrans,NnOut,d,Nc,1.0,BCoords_complete,rst_c);
	// keep (requires external free)
	free(BCoords_complete);

//array_print_d(NnOut,d,rstOut,'C');

	NsOut = 0;
	for (i = 0; i < 2; i++)
		NsOut += symms31[i];

	symmsOut = malloc(NsOut * sizeof *symmsOut); // keep (requires external free)
	k = 0;
	for (i = 0; i < 2; i++) {
	for (j = 0; jMax = symms31[i], j < jMax; j++) {
		if (i == 0) symmsOut[k] = 3;
		else        symmsOut[k] = 1;
		k++;
	}}

	*rst   = rstOut;
	*symms = symmsOut;

	*Nn = NnOut;
	*Ns = NsOut;

	if (return_w) {
		free(w_read);
		*w = wOut;
	} else {
		*w = NULL;
	}
}
