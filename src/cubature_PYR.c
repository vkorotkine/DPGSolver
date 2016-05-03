#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "parameters.h"
#include "functions.h"

/*
 *	Purpose:
 *		Return nodes, weights and symmetries for pyramidal cubature depending on the nodetype.
 *
 *	Comments:
 *		Barycentric coordinates for the PYR element may be negative, unlike for the SI elements. This does not pose any
 *		challenges or result in any limitations, but is simply a result of the last additional condition (chosen for
 *		good conditioning of the matrix) when computing the barycentric coordinates of the nodes in the pyfr code.
 *
 *		Note:
 *
 *			[ones]_{Nn x 1}  = [BCoords]_{Nn x Nc} * [ones]_{Nc x 1}          => Partition of unity (1 condition)
 *			[rst]_{Nn x d}   = [BCoords]_{Nn x Nc} * [rst_(c)orners]_{Nc x d} => Linear precision   (d conditions)
 *			[r*s*t]_{Nn x 1} = [BCoords]_{Nn x Nc} * [{r*s*t}_c]_{Nc x 1}     => Arbitrary          (1 condition)
 *
 *			Then
 *				[BCoords] = [ones(Nn,1) rst r*s*t]*inv([ones(Nc,1) rst_c {r*s*t}_c])
 *
 *			where
 *				cond([ones(Nc,1) rst_c {r*s*t}_c]) = 2.0
 *
 *		All nodes were determined based off of those from the pyfr code (pyfr/quadrules/pyr) after being transfered to
 *		the regular PYR used in this code. (ToBeModified if other nodes are added)
 *			Options:
 *				GL     : GL nodes,                       no weights
 *				GLL    : GLL nodes,                      no weights
 *				GLW    : GL nodes,                       with weights
 *				GLLW   : GLL nodes,                      with weights
 *				WV     : WV PYR nodes,                   with weights
 *					Exact integration to lower order than expected.
 *				WVHToP : WV HEX nodes transfered to PYR, with weights
 *					Exact integration to lower order than expected.
 *
 *			After implementing the traditional PYR orthogonal basis
 *			(Chan(2015)-Orthogonal_Bases_for_Vertex_Mapped_Pyramids, eq. 2.1) and the orthogonal basis from the pyfr
 *			code (Witherden(2015,Thesis), eq. 3.20), the following conclusions were drawn:
 *				Despite being lower-order in the "c" term, the traditional basis mass matrix is not integrated exactly
 *				using the WV nodes, while the basis in the pyfr code is (eventually); this is odd. Further, given the
 *				maximum strength WV rule (WV10), even the P3 PYR mass matrix is not given exactly with constant
 *				Jacobian. Further testing revealed that the WV nodes integrate exactly (with the expected order) only
 *				when there is variation either in a OR in b, but not in both; the nodes cannot exactly integrate
 *				polynomials on QUAD cross-sections of the pyramid, unlike the GL nodes.
 *
 *				Using GL nodes transfered to the PYR element, the cubature strength for exact mass matrix is as
 *				expected, after accounting for the added contribution to the "c" term from the weight (w =
 *				w_HEX*pow(1-c,2); GL nodes of order P+1 are required for exact integration of all terms). Using the WV
 *				HEX nodes transfered to the PYR element does not result in exact integration (similarly to the PYR
 *				nodes) either for the traditional or pyfr PYR basis.
 *
 *		The order of rst, w is important for minimizing memory stride while computing the length 4 Discrete Fourier
 *		Transform (Note: as w is a 1d matrix, its ordering is actually not relevant) (ToBeModified).
 *		Ordering convention:
 *			4-blocks of symmetric nodes in the same rotational order, followed by 1-block(s) of center node(s) if
 *			present.
 *		rst is stored in memory as r, s, then t (See cubature_TP for the motivation).
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
 *		ToBeModified: Add reference to Witherden(2015)-On_the_Identification_of_Symmetric_...
 *		pyfr code : http://www.pyfr.org
 *
 *		GL  : pyfr/quadrules/pyr + conversion to barycentric coordinates (ToBeModified: See python script)
 *		GLL : pyfr/quadrules/pyr + conversion to barycentric coordinates (ToBeModified: See python script)
 *		WV  : pyfr/quadrules/pyr + conversion to barycentric coordinates (ToBeModified: See python script)
 */

void cubature_PYR(double **rst, double **w, unsigned int **symms, unsigned int *Nn, unsigned int *Ns,
                  const unsigned int return_w, const unsigned int P, const unsigned int d, const char *NodeType)
{
	// Standard datatypes
	static unsigned int perms_QUAD[16] = { 0, 1, 2, 3,
	                                       3, 0, 1, 2,
	                                       2, 3, 0, 1,
	                                       1, 2, 3, 0};
	unsigned int symms41[2] = { 0, 0}, NQUADsymms = 0;
	double rst_c[15] = { -1.0,            1.0,            1.0,           -1.0,           0.0,
	                     -1.0,           -1.0,            1.0,            1.0,           0.0,
	                     -sqrt(2.0)/5.0, -sqrt(2.0)/5.0, -sqrt(2.0)/5.0, -sqrt(2.0)/5.0, 4.0/5.0*sqrt(2.0)},
	       BCoords_tmp[5];
	unsigned int i, iMax, j, jMax, k, QUADsymm,
	             IndB, IndBC, IndGroup, Ind1, GroupCount, Nc, N1,
	             PMax, NnOut, NsOut, Ngroups, Nsymms, Nperms;
	unsigned int *symmsOut, *symms_Nperms, *symms_count;
	char         *StringRead, *strings, *stringe, *CubFile, *Pc;
	double       *rstOut, *wOut, *BCoords, *BCoords_complete, *w_read;

	FILE *fID;

	// Silence compiler warnings
	PMax = 0; Nsymms = 0; Ngroups = 0;
	w_read = malloc(0 * sizeof *w_read); // silence
	wOut = NULL;

	if (strstr(NodeType,"GL") != NULL) {
		if (strstr(NodeType,"W") == NULL) {
			if (return_w)
				printf("Error: Invalid value for return_w in cubature_PYR.\n"), exit(1);

			PMax = 6;
		} else {
			PMax = 6;
		}
	} else if (strstr(NodeType,"WV") != NULL) {
		if (strstr(NodeType,"HToP") == NULL)
			PMax = 10;
		else
			PMax = 11;
	}

	if (P > PMax)
		printf("Error: %s PYR nodes of order %d are not available.\n",NodeType,P), exit(1);

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

			if (return_w)
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
					if (return_w)
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
//array_print_d(NnOut,Nc,BCoords_complete,'R');

	free(symms_count);
	free(symms_Nperms);
	free(BCoords);
	free(w_read);

	rstOut = mm_Alloc_d(CblasColMajor,CblasTrans,CblasNoTrans,NnOut,d,Nc,1.0,BCoords_complete,rst_c);
	// keep (requires external free)
	free(BCoords_complete);

//array_print_d(NnOut,d,rstOut,'C');

	NsOut = 0;
	for (i = 0; i < 2; i++)
		NsOut += symms41[i];

	symmsOut = malloc(NsOut * sizeof *symmsOut); // keep (requires external free)
	k = 0;
	for (i = 0; i < 2; i++) {
	for (j = 0; jMax = symms41[i], j < jMax; j++) {
		if (i == 0) symmsOut[k] = 4;
		else        symmsOut[k] = 1;
		k++;
	}}

	*rst   = rstOut;
	*symms = symmsOut;

	*Nn = NnOut;
	*Ns = NsOut;

	if (return_w) {
		*w = wOut;
	} else {
		free(wOut);
		*w = NULL;
	}
}
