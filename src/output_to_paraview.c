// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "output_to_paraview.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_ELEMENT.h"
#include "S_VOLUME.h"
#include "S_FACE.h"

#include "element_functions.h"
#include "matrix_functions.h"
#include "variable_functions.h"

#include "exact_solutions.h"
#include "array_print.h"

/*
 *	Purpose:
 *		Output data to paraview for visualization.
 *
 *	Comments:
 *		Getting an error from paraview when using less than 3 components for the Points DataArray. (2016-04-21)
 *		Potentially place all of the elements in a single piece to reduce file sizes. (ToBeDeleted)
 *		Possibly convert output to binary if this is found to be slow when profiling. (ToBeDeleted)
 *		The following outputs are supported: (ToBeModified)
 *			G : Geometry only
 *
 *		Mesh visualization in paraview:
 *			0) Load .vtk file
 *			1) Change 'Representation' to 'wireframe'
 *			2) Split the screen and select 'SpreadSheet View' in the second window. Change 'Attribute' to 'Cell Data'.
 *			3) From the 'View' menu, select 'Selection Display Inspector' and enable Cell/Point Label IDs.
 *
 *			Any selection of ELEMENTs in the spreadsheet should not show Cell and Point Labels in the first screen.
 *			Note: Use P1 ELEMENTs so that the cells in paraview are the ELEMENTs themselves.
 *
 *	Notation:
 *
 *	References:
 */

static unsigned int get_ELEMENT_Ncorners(const unsigned int VTK_type)
{
	switch (VTK_type) {
	case 3: // LINE
		return 2;
		break;
	case 5: // TRI
		return 3;
		break;
	case 9: // QUAD
		return 4;
		break;
	case 10: // TET
		return 4;
		break;
	case 12: // HEX
		return 8;
		break;
	case 13: // WEDGE
		return 6;
		break;
	case 14: // PYR
		return 5;
		break;
	default:
		printf("Error: Unsupported VTK_type.\n"), exit(1);
		break;
	}
}

static void fprintf_tn(FILE *fID, unsigned int Ntabs, const char *String)
{
	unsigned int i;

	for (i = 0; i < Ntabs; i++)
		fprintf(fID,"\t");

	fprintf(fID,"%s\n",String);
}

static void output_geom(const char *geom_type)
{
	// Initialize DB Parameters
//	char         *TestCase = DB.TestCase;
	unsigned int d         = DB.d;
	int          MPIrank   = DB.MPIrank,
	             MPIsize   = DB.MPIsize;

	// standard datatypes
	char MPIrank_c[STRLEN_MIN], f_name[STRLEN_MAX], f_parallel[STRLEN_MAX], f_serial[STRLEN_MAX];
	unsigned int i, iMax, j, jMax, dim, sum,
	             P, level, NE, NvnP, NvnG,
	             *connectivity, *types, *VTK_Ncorners;
	double *I_vG_vP, *XYZ_vP, *Input;
	FILE *fID;

	struct S_ELEMENT *ELEMENT;
	struct S_VOLUME *VOLUME;

	// silence
	NvnP = 0;
	XYZ_vP = NULL;

	sprintf(MPIrank_c,"%d",MPIrank);
//	strcpy(f_name,TestCase);
//	strcat(f_name,"_geom");
	strcpy(f_name,geom_type);

	if (!DB.MPIrank) {
		strcpy(f_parallel,"paraview/");
		strcat(f_parallel,f_name);
		strcat(f_parallel,".pvtu");

		if ((fID = fopen(f_parallel,"w")) == NULL)
			printf("Error: File: %s, did not open.\n",f_parallel), exit(1);

		fprintf_tn(fID,0,"<?xml version=\"1.0\"?>");
		fprintf_tn(fID,0,"<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">");
		fprintf_tn(fID,1,"<PUnstructuredGrid GhostLevel=\"0\">\n");

		fprintf_tn(fID,2,"<PPointData Scalars=\"Scalars\" Vectors=\"Vectors\">");
		fprintf_tn(fID,3,"<PDataArray type=\"UInt8\" Name=\"P\" format=\"ascii\"/>");
		fprintf_tn(fID,3,"<PDataArray type=\"UInt8\" Name=\"level\" format=\"ascii\"/>");
		if (strstr(geom_type,"adapt")) {
			fprintf_tn(fID,3,"<PDataArray type=\"UInt8\" Name=\"Vadapt\" format=\"ascii\"/>");
			fprintf_tn(fID,3,"<PDataArray type=\"Int8\" Name=\"adapt_type\" format=\"ascii\"/>");
			fprintf_tn(fID,3,"<PDataArray type=\"UInt32\" Name=\"indexg\" format=\"ascii\"/>");
		}
		fprintf_tn(fID,2,"</PPointData>\n");

		fprintf_tn(fID,2,"<PPoints>");
		fprintf_tn(fID,3,"<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\"/>");
		fprintf_tn(fID,2,"</PPoints>\n");

		fprintf_tn(fID,2,"<PCells>");
		fprintf_tn(fID,3,"<PDataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\"/>");
		fprintf_tn(fID,3,"<PDataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\"/>");
		fprintf_tn(fID,3,"<PDataArray type=\"UInt8\" Name=\"types\" format=\"ascii\"/>");
		fprintf_tn(fID,2,"</PCells>\n");

		for (i = 0, iMax = (unsigned int) MPIsize; i < iMax; i++)
			fprintf(fID,"\t\t<Piece Source=\"%s%d.vtu\"/>\n",f_name,i);

		fprintf_tn(fID,1,"</PUnstructuredGrid>");
		fprintf_tn(fID,0,"</VTKFile>");

		fclose(fID);
	}

	strcpy(f_serial,"paraview/");
	strcat(f_serial,f_name);
	strcat(f_serial,MPIrank_c);
	strcat(f_serial,".vtu");

	if ((fID = fopen(f_serial,"w")) == NULL)
		printf("Error: File f_serial did not open.\n"), exit(1);

	fprintf_tn(fID,0,"<?xml version=\"1.0\"?>");
	fprintf_tn(fID,0,"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">");
	fprintf_tn(fID,1,"<UnstructuredGrid>\n");

	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		ELEMENT = get_ELEMENT_type(VOLUME->type);

		P     = VOLUME->P;
		level = VOLUME->level;

		connectivity = ELEMENT->connectivity[P];
		types        = ELEMENT->connect_types[P];
		NE           = ELEMENT->connect_NE[P];

		NvnP = ELEMENT->NvnP[P];
		NvnG = VOLUME->NvnG;

		if (!VOLUME->curved)
			I_vG_vP = ELEMENT->I_vGs_vP[1][P][0];
		else
			I_vG_vP = ELEMENT->I_vGc_vP[P][P][0];

		if (strstr(geom_type,"straight"))
			Input = VOLUME->XYZ_S;
		else
			Input = VOLUME->XYZ;

		XYZ_vP = mm_Alloc_d(CBCM,CBT,CBNT,NvnP,d,NvnG,1.0,I_vG_vP,Input); // free

		fprintf(fID,"\t\t<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n",NvnP,NE);

			fprintf_tn(fID,3,"<Points>");
				fprintf(fID,"\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"%d\" format=\"ascii\">\n",3);
				for (i = 0; i < NvnP; i++) {
					fprintf(fID,"\t\t\t\t");
					for (dim = 0; dim < d; dim++)
						fprintf(fID,"% .3f ",XYZ_vP[dim*NvnP+i]);
					for (dim = d; dim < DMAX; dim++)
						fprintf(fID,"% .3f ",0.0);
					fprintf(fID,"\n");
				}
				fprintf_tn(fID,4,"</DataArray>");
			fprintf_tn(fID,3,"</Points>");

			fprintf_tn(fID,3,"<PointData Scalars=\"Scalars\" Vectors=\"Vectors\">");
				fprintf_tn(fID,4,"<DataArray type=\"UInt8\" Name=\"P\" format=\"ascii\">");
				fprintf(fID,"\t\t\t\t");
				for (i = 0; i < NvnP; i++) {
					fprintf(fID,"%d ",P);
					if ((i+1) % 5 == 0 && i != NvnP-1)
						fprintf(fID,"\n\t\t\t\t");
					else if (i == NvnP-1)
						fprintf(fID,"\n");
				}
				fprintf_tn(fID,4,"</DataArray>");
				fprintf_tn(fID,4,"<DataArray type=\"UInt8\" Name=\"level\" format=\"ascii\">");
				fprintf(fID,"\t\t\t\t");
				for (i = 0; i < NvnP; i++) {
					fprintf(fID,"%d ",level);
					if ((i+1) % 5 == 0 && i != NvnP-1)
						fprintf(fID,"\n\t\t\t\t");
					else if (i == NvnP-1)
						fprintf(fID,"\n");
				}
				fprintf_tn(fID,4,"</DataArray>");

				if (strstr(geom_type,"adapt")) {
					fprintf_tn(fID,4,"<DataArray type=\"UInt8\" Name=\"Vadapt\" format=\"ascii\">");
					fprintf(fID,"\t\t\t\t");
					for (i = 0; i < NvnP; i++) {
						fprintf(fID,"%d ",VOLUME->Vadapt);
						if ((i+1) % 5 == 0 && i != NvnP-1)
							fprintf(fID,"\n\t\t\t\t");
						else if (i == NvnP-1)
							fprintf(fID,"\n");
					}
					fprintf_tn(fID,4,"</DataArray>");
					fprintf_tn(fID,4,"<DataArray type=\"Int8\" Name=\"adapt_type\" format=\"ascii\">");
					fprintf(fID,"\t\t\t\t");
					for (i = 0; i < NvnP; i++) {
						fprintf(fID,"%d ",VOLUME->adapt_type);
						if ((i+1) % 5 == 0 && i != NvnP-1)
							fprintf(fID,"\n\t\t\t\t");
						else if (i == NvnP-1)
							fprintf(fID,"\n");
					}
					fprintf_tn(fID,4,"</DataArray>");
					fprintf_tn(fID,4,"<DataArray type=\"UInt32\" Name=\"indexg\" format=\"ascii\">");
					fprintf(fID,"\t\t\t\t");
					for (i = 0; i < NvnP; i++) {
						fprintf(fID,"%d ",VOLUME->indexg);
						if ((i+1) % 5 == 0 && i != NvnP-1)
							fprintf(fID,"\n\t\t\t\t");
						else if (i == NvnP-1)
							fprintf(fID,"\n");
					}
					fprintf_tn(fID,4,"</DataArray>");
				}
			fprintf_tn(fID,3,"</PointData>");

			fprintf_tn(fID,3,"<Cells>");

			VTK_Ncorners = malloc(NE * sizeof * VTK_Ncorners); // free
			for (i = 0; i < NE; i++) {
				VTK_Ncorners[i] = get_ELEMENT_Ncorners(types[i]);
			}

				fprintf_tn(fID,4,"<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">");
				for (i = 0; i < NE; i++) {
					fprintf(fID,"\t\t\t\t");
					for (j = 0, jMax = VTK_Ncorners[i]; j < jMax; j++)
						fprintf(fID,"%6d ",connectivity[i*8+j]);
					fprintf(fID,"\n");
				}
				fprintf_tn(fID,4,"</DataArray>");

				fprintf_tn(fID,4,"<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">");
				fprintf(fID,"\t\t\t\t");
				for (i = 0, sum = 0; i < NE; i++) {
					sum += VTK_Ncorners[i];
					fprintf(fID,"%6d ",sum);
					if ((i+1) % 8 == 0 && i != NE-1)
						fprintf(fID,"\n\t\t\t\t");
					else if (i == NE-1)
						fprintf(fID,"\n");
				}
				fprintf_tn(fID,4,"</DataArray>");

				fprintf_tn(fID,4,"<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">");
				fprintf(fID,"\t\t\t\t");
				for (i = 0; i < NE; i++) {
					fprintf(fID,"%6d ",types[i]);
					if ((i+1) % 8 == 0 && i != NE-1)
						fprintf(fID,"\n\t\t\t\t");
					else if (i == NE-1)
						fprintf(fID,"\n");
				}
				fprintf_tn(fID,4,"</DataArray>");

			fprintf_tn(fID,3,"</Cells>");

		fprintf_tn(fID,2,"</Piece>\n");

		free(XYZ_vP);
		free(VTK_Ncorners);
	}

	fprintf_tn(fID,1,"</UnstructuredGrid>");
	fprintf(fID,"</VTKFile>");

	fclose(fID);
}

static void output_normals(const char *normals_type)
{
	// Initialize DB Parameters
//	char         *TestCase = DB.TestCase;
	unsigned int d         = DB.d,
	             NfrefMax  = DB.NfrefMax;
	int          MPIrank   = DB.MPIrank,
	             MPIsize   = DB.MPIsize;

	// Standard datatypes
	char         MPIrank_c[STRLEN_MIN], f_name[STRLEN_MAX], f_parallel[STRLEN_MAX], f_serial[STRLEN_MAX];
	unsigned int i, iMax, dim, nInd, curved, PV, PF, Nfn, NvnG, IndFType, Eclass, VfL;
	double       *Input, *I_vG_f, *XYZ_f, *n;
	FILE         *fID;

	struct S_ELEMENT *ELEMENT;
	struct S_VOLUME  *VL;
	struct S_FACE   *FACE;

	// silence
	XYZ_f = NULL;
	Nfn = 0;

	sprintf(MPIrank_c,"%d",MPIrank);
//	strcpy(f_name,TestCase);
//	strcat(f_name,"_normals");
	strcpy(f_name,normals_type);

	if (!DB.MPIrank) {
		strcpy(f_parallel,"paraview/");
		strcat(f_parallel,f_name);
		strcat(f_parallel,".pvtp");

		if ((fID = fopen(f_parallel,"w")) == NULL)
			printf("Error: File: %s, did not open.\n",f_parallel), exit(1);

		fprintf_tn(fID,0,"<?xml version=\"1.0\"?>");
		fprintf_tn(fID,0,"<VTKFile type=\"PPolyData\" version=\"0.1\" byte_order=\"LittleEndian\">");
		fprintf_tn(fID,1,"<PPolyData GhostLevel=\"0\">\n");

		fprintf_tn(fID,2,"<PPointData Vectors=\"Normals\">");
		fprintf_tn(fID,3,"<PDataArray type=\"Float32\" Name=\"Normals\" NumberOfComponents=\"3\" format=\"ascii\"/>");
		fprintf_tn(fID,2,"</PPointData>\n");

		fprintf_tn(fID,2,"<PPoints>");
		fprintf_tn(fID,3,"<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\"/>");
		fprintf_tn(fID,2,"</PPoints>\n");

		for (i = 0, iMax = (unsigned int) MPIsize; i < iMax; i++)
			fprintf(fID,"\t\t<Piece Source=\"%s%d.vtp\"/>\n",f_name,i);

		fprintf_tn(fID,1,"</PPolyData>");
		fprintf_tn(fID,0,"</VTKFile>");

		fclose(fID);
	}

	strcpy(f_serial,"paraview/");
	strcat(f_serial,f_name);
	strcat(f_serial,MPIrank_c);
	strcat(f_serial,".vtp");

	if ((fID = fopen(f_serial,"w")) == NULL)
		printf("Error: File f_serial did not open.\n"), exit(1);

	fprintf_tn(fID,0,"<?xml version=\"1.0\"?>");
	fprintf_tn(fID,0,"<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">");
	fprintf_tn(fID,1,"<PolyData>\n");

	for (FACE = DB.FACE; FACE; FACE = FACE->next) {
		curved = FACE->curved;

		VL  = FACE->VL;
		VfL = FACE->VfL;

		if (VfL % NFREFMAX != 0)
			printf("Error: VfL should be h-conforming in output_to_paraview (output_normals).\n"), exit(1);

		PV = VL->P;
		PF = FACE->P;

		ELEMENT = get_ELEMENT_type(VL->type);
		Eclass  = get_Eclass(ELEMENT->type);

		IndFType = get_IndFType(Eclass,VfL/NfrefMax);

		NvnG = VL->NvnG;

		Input = VL->XYZ;

		if (FACE->typeInt == 's') {
			Nfn = ELEMENT->NfnIs[PF][IndFType];
			if (!VL->curved) I_vG_f = ELEMENT->I_vGs_fIs[1][PF][VfL];
			else              I_vG_f = ELEMENT->I_vGc_fIs[PV][PF][VfL];
		} else {
			Nfn = ELEMENT->NfnIc[PF][IndFType];
			if (!VL->curved) I_vG_f = ELEMENT->I_vGs_fIc[1][PF][VfL];
			else              I_vG_f = ELEMENT->I_vGc_fIc[PV][PF][VfL];
		}
		n = FACE->n_fI;
		XYZ_f = mm_Alloc_d(CBCM,CBT,CBNT,Nfn,d,NvnG,1.0,I_vG_f,Input); // free

		fprintf(fID,"\t\t<Piece NumberOfPoints=\"%d\" NumberOfVerts=\"%d\" NumberOfLines=\"%d\" NumberOfStrips=\"%d\" "
		            "NumberOfPolys=\"%d\">\n",Nfn,0,0,0,0);

			fprintf_tn(fID,3,"<Points>");
				fprintf(fID,"\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"%d\" format=\"ascii\">\n",3);
				for (i = 0; i < Nfn; i++) {
					fprintf(fID,"\t\t\t\t");
					for (dim = 0; dim < d; dim++)
						fprintf(fID,"% .3f ",XYZ_f[dim*Nfn+i]);
					for (dim = d; dim < 3; dim++)
						fprintf(fID,"% .3f ",0.0);
					fprintf(fID,"\n");
				}
				fprintf_tn(fID,4,"</DataArray>");
			fprintf_tn(fID,3,"</Points>");

			fprintf_tn(fID,3,"<PointData Vectors=\"Normals\">");
				fprintf(fID,"\t\t\t\t<DataArray type=\"Float32\" Name=\"Normals\" NumberOfComponents=\"%d\" "
				                 "format=\"ascii\">\n",3);
				for (i = 0; i < Nfn; i++) {
					fprintf(fID,"\t\t\t\t");
					for (dim = 0; dim < d; dim++) {
						if (!curved) nInd = dim;
						else         nInd = i*d+dim;
						fprintf(fID,"% .3f ",n[nInd]);
					}
					for (dim = d; dim < 3; dim++)
						fprintf(fID,"% .3f ",0.0);
					fprintf(fID,"\n");
				}
				fprintf_tn(fID,4,"</DataArray>");
			fprintf_tn(fID,3,"</PointData>");

		fprintf_tn(fID,2,"</Piece>\n");

		free(XYZ_f);
	}

	fprintf_tn(fID,1,"</PolyData>");
	fprintf(fID,"</VTKFile>");

	fclose(fID);
}

static void output_solution(const char *sol_type)
{
	// Initialize DB Parameters
	char         *TestCase = DB.TestCase,
	             *MeshType = DB.MeshType;
	unsigned int d         = DB.d,
	             Nvar      = DB.Nvar;
	int          MPIrank   = DB.MPIrank,
	             MPIsize   = DB.MPIsize;

	// standard datatypes
	char         MPIrank_c[STRLEN_MIN], f_name[STRLEN_MAX], f_name_source[STRLEN_MAX],
	             f_parallel[STRLEN_MAX], f_serial[STRLEN_MAX];
	unsigned int i, iMax, j, jMax, dim, sum, varMax,
	             P, PP, NE, NvnP, NvnG, NvnS,
	             *connectivity, *types, *VTK_Ncorners;
	double       *I_vG_vP, *ChiS_vP, *XYZ_vP, *W_vP, *U_vP, *rho, *u, *v, *w, *p, *E, *s, *Mach, V2, c2, *q;
	FILE *fID;

	struct S_ELEMENT *ELEMENT;
	struct S_VOLUME *VOLUME;

	// silence
	W_vP = U_vP = NULL;
	rho = u = v = w = p = E = s = Mach = q = NULL;

	sprintf(MPIrank_c,"%d",MPIrank);
	strcpy(f_name,TestCase); strcat(f_name,"/");
	strcat(f_name,MeshType); strcat(f_name,"/");
	strcat(f_name,sol_type);
	strcpy(f_name_source,sol_type);

	if (!DB.MPIrank) {
		strcpy(f_parallel,"paraview/");
		strcat(f_parallel,f_name);
		strcat(f_parallel,".pvtu");

		if ((fID = fopen(f_parallel,"w")) == NULL)
			printf("Error: File: %s, did not open.\n",f_parallel), exit(1);

		fprintf_tn(fID,0,"<?xml version=\"1.0\"?>");
		fprintf_tn(fID,0,"<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">");
		fprintf_tn(fID,1,"<PUnstructuredGrid GhostLevel=\"0\">\n");

		fprintf_tn(fID,2,"<PPointData Scalars=\"Scalars\" Vectors=\"Vectors\">");
		if (strstr(TestCase,"Advection")) {
			fprintf_tn(fID,3,"<PDataArray type=\"Float32\" Name=\"u\" format=\"ascii\"/>");
		} else if (strstr(TestCase,"Poisson")) {
			fprintf_tn(fID,3,"<PDataArray type=\"Float32\" Name=\"u\" format=\"ascii\"/>");
			fprintf_tn(fID,3,"<PDataArray type=\"Float32\" Name=\"q\" NumberOfComponents=\"3\" format=\"ascii\"/>");
		} else if (strstr(TestCase,"Euler") || strstr(TestCase,"NavierStokes")) {
			fprintf_tn(fID,3,"<PDataArray type=\"Float32\" Name=\"rho\" format=\"ascii\"/>");
			fprintf_tn(fID,3,"<PDataArray type=\"Float32\" Name=\"V\" NumberOfComponents=\"3\" format=\"ascii\"/>");
			fprintf_tn(fID,3,"<PDataArray type=\"Float32\" Name=\"p\" format=\"ascii\"/>");
			fprintf_tn(fID,3,"<PDataArray type=\"Float32\" Name=\"E\" format=\"ascii\"/>");
			fprintf_tn(fID,3,"<PDataArray type=\"Float32\" Name=\"s\" format=\"ascii\"/>");
			fprintf_tn(fID,3,"<PDataArray type=\"Float32\" Name=\"Mach\" format=\"ascii\"/>");
		} else {
			EXIT_UNSUPPORTED;
		}
		fprintf_tn(fID,2,"</PPointData>\n");

		fprintf_tn(fID,2,"<PPoints>");
		fprintf_tn(fID,3,"<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\"/>");
		fprintf_tn(fID,2,"</PPoints>\n");

		fprintf_tn(fID,2,"<PCells>");
		fprintf_tn(fID,3,"<PDataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\"/>");
		fprintf_tn(fID,3,"<PDataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\"/>");
		fprintf_tn(fID,3,"<PDataArray type=\"UInt8\" Name=\"types\" format=\"ascii\"/>");
		fprintf_tn(fID,2,"</PCells>\n");

		for (i = 0, iMax = (unsigned int) MPIsize; i < iMax; i++)
			fprintf(fID,"\t\t<Piece Source=\"%s%d.vtu\"/>\n",f_name_source,i);

		fprintf_tn(fID,1,"</PUnstructuredGrid>");
		fprintf_tn(fID,0,"</VTKFile>");

		fclose(fID);
	}

	strcpy(f_serial,"paraview/");
	strcat(f_serial,f_name);
	strcat(f_serial,MPIrank_c);
	strcat(f_serial,".vtu");

	if ((fID = fopen(f_serial,"w")) == NULL)
		printf("Error: File f_serial did not open.\n"), exit(1);

	fprintf_tn(fID,0,"<?xml version=\"1.0\"?>");
	fprintf_tn(fID,0,"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">");
	fprintf_tn(fID,1,"<UnstructuredGrid>\n");

	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		ELEMENT = get_ELEMENT_type(VOLUME->type);

		P = VOLUME->P;

		// Visualize with higher order to see curved geometry if applicable (P+1 is the maximum currently supported)
		if (!VOLUME->curved || DB.Adapt == ADAPT_0 || DB.Adapt == ADAPT_H)
			PP = P;
		else
			PP = P;
//			PP = min(PMax,P+1);

		connectivity = ELEMENT->connectivity[PP];
		types        = ELEMENT->connect_types[PP];
		NE           = ELEMENT->connect_NE[PP];

		NvnP = ELEMENT->NvnP[PP];
		NvnG = VOLUME->NvnG;
		NvnS = VOLUME->NvnS;

		if (!VOLUME->curved)
			I_vG_vP = ELEMENT->I_vGs_vP[1][P][0];
		else
			I_vG_vP = ELEMENT->I_vGc_vP[P][PP][0];
		ChiS_vP = ELEMENT->ChiS_vP[P][PP][0];

		XYZ_vP = mm_Alloc_d(CBCM,CBT,CBNT,NvnP,d,NvnG,1.0,I_vG_vP,VOLUME->XYZ);     // free

		if (strstr(TestCase,"Advection")) {
			u = mm_Alloc_d(CBCM,CBT,CBNT,NvnP,Nvar,NvnS,1.0,ChiS_vP,VOLUME->What); // free
		} else if (strstr(TestCase,"Poisson")) {
			u = mm_Alloc_d(CBCM,CBT,CBNT,NvnP,Nvar,NvnS,1.0,ChiS_vP,VOLUME->What); // free
			q = calloc(DMAX*NvnP , sizeof *q); // free
			for (dim = 0; dim < d; dim++)
				mm_d(CBCM,CBT,CBNT,NvnP,Nvar,NvnS,1.0,0.0,ChiS_vP,VOLUME->Qhat[dim],&q[dim*NvnP]);

			// Store solution error in q3 for d = 2
			if (d == 2) {
				unsigned int n;
				double       *uEx;

				uEx = malloc(NvnP*1 * sizeof *uEx); // free
				compute_exact_solution(NvnP,XYZ_vP,uEx,1);
				for (n = 0; n < NvnP; n++)
					q[d*NvnP+n] = fabs(u[n]-uEx[n]);

				free(uEx);
			}
		} else if (strstr(TestCase,"Euler") || strstr(TestCase,"NavierStokes")) {
			W_vP = mm_Alloc_d(CBCM,CBT,CBNT,NvnP,Nvar,NvnS,1.0,ChiS_vP,VOLUME->What); // free

			U_vP = malloc(NvnP*5 * sizeof *U_vP); // free

			convert_variables(W_vP,U_vP,d,3,NvnP,1,'c','p');
			varMax = Nvar-1;

			rho = &U_vP[NvnP*0];
			u   = &U_vP[NvnP*1];
			v   = &U_vP[NvnP*2];
			w   = &U_vP[NvnP*3];
			p   = &U_vP[NvnP*4];
			E   = &W_vP[NvnP*varMax];

			s    = malloc(NvnP * sizeof *s);    // free
			Mach = malloc(NvnP * sizeof *Mach); // free

			for (i = 0; i < NvnP; i++) {
				V2 = u[i]*u[i] + v[i]*v[i] + w[i]*w[i];
				c2 = GAMMA*p[i]/rho[i];
				if (c2 < 0.0)
					c2 = EPS;

				s[i]    = p[i]/pow(rho[i],GAMMA);
				Mach[i] = sqrt(V2/c2);
			}
		} else {
			EXIT_UNSUPPORTED;
		}

		fprintf(fID,"\t\t<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n",NvnP,NE);

			fprintf_tn(fID,3,"<Points>");
				fprintf(fID,"\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"%d\" format=\"ascii\">\n",3);
				for (i = 0; i < NvnP; i++) {
					fprintf(fID,"\t\t\t\t");
					for (dim = 0; dim < d; dim++)
						fprintf(fID,"% .4e ",XYZ_vP[dim*NvnP+i]);
					for (dim = d; dim < DMAX; dim++)
						fprintf(fID,"% .4e ",0.0);
					fprintf(fID,"\n");
				}
				fprintf_tn(fID,4,"</DataArray>");
			fprintf_tn(fID,3,"</Points>");

			fprintf_tn(fID,3,"<PointData Scalars=\"Scalars\" Vectors=\"Vectors\">");
			if (strstr(TestCase,"Advection")) {
				fprintf_tn(fID,4,"<DataArray type=\"Float32\" Name=\"u\" format=\"ascii\">");
				fprintf(fID,"\t\t\t\t");
				for (i = 0; i < NvnP; i++) {
					fprintf(fID,"% .8e ",u[i]);
					if ((i+1) % 5 == 0 && i != NvnP-1)
						fprintf(fID,"\n\t\t\t\t");
					else if (i == NvnP-1)
						fprintf(fID,"\n");
				}
				fprintf_tn(fID,4,"</DataArray>");
			} else if (strstr(TestCase,"Poisson")) {
				fprintf_tn(fID,4,"<DataArray type=\"Float32\" Name=\"u\" format=\"ascii\">");
				fprintf(fID,"\t\t\t\t");
				for (i = 0; i < NvnP; i++) {
					fprintf(fID,"% .8e ",u[i]);
					if ((i+1) % 5 == 0 && i != NvnP-1)
						fprintf(fID,"\n\t\t\t\t");
					else if (i == NvnP-1)
						fprintf(fID,"\n");
				}
				fprintf_tn(fID,4,"</DataArray>");

				fprintf(fID,"\t\t\t\t<DataArray type=\"Float32\" Name=\"q\" NumberOfComponents=\"%d\" format=\"ascii\">\n",3);
				for (i = 0; i < NvnP; i++)
					fprintf(fID,"\t\t\t\t % .8e % .8e % .8e\n",q[NvnP*0+i],q[NvnP*1+i],q[NvnP*2+i]);
				fprintf_tn(fID,4,"</DataArray>");
			} else if (strstr(TestCase,"Euler") || strstr(TestCase,"NavierStokes")) {
				fprintf_tn(fID,4,"<DataArray type=\"Float32\" Name=\"rho\" format=\"ascii\">");
				fprintf(fID,"\t\t\t\t");
				for (i = 0; i < NvnP; i++) {
					fprintf(fID,"% .8e ",rho[i]);
					if ((i+1) % 5 == 0 && i != NvnP-1)
						fprintf(fID,"\n\t\t\t\t");
					else if (i == NvnP-1)
						fprintf(fID,"\n");
				}
				fprintf_tn(fID,4,"</DataArray>");

				fprintf(fID,"\t\t\t\t<DataArray type=\"Float32\" Name=\"V\" NumberOfComponents=\"%d\" format=\"ascii\">\n",3);
				for (i = 0; i < NvnP; i++)
					fprintf(fID,"\t\t\t\t % .8e % .8e % .8e\n",u[i],v[i],w[i]);
				fprintf_tn(fID,4,"</DataArray>");

				fprintf_tn(fID,4,"<DataArray type=\"Float32\" Name=\"p\" format=\"ascii\">");
				fprintf(fID,"\t\t\t\t");
				for (i = 0; i < NvnP; i++) {
					fprintf(fID,"% .8e ",p[i]);
					if ((i+1) % 5 == 0 && i != NvnP-1)
						fprintf(fID,"\n\t\t\t\t");
					else if (i == NvnP-1)
						fprintf(fID,"\n");
				}
				fprintf_tn(fID,4,"</DataArray>");

				fprintf_tn(fID,4,"<DataArray type=\"Float32\" Name=\"E\" format=\"ascii\">");
				fprintf(fID,"\t\t\t\t");
				for (i = 0; i < NvnP; i++) {
					fprintf(fID,"% .8e ",E[i]);
					if ((i+1) % 5 == 0 && i != NvnP-1)
						fprintf(fID,"\n\t\t\t\t");
					else if (i == NvnP-1)
						fprintf(fID,"\n");
				}
				fprintf_tn(fID,4,"</DataArray>");

				fprintf_tn(fID,4,"<DataArray type=\"Float32\" Name=\"s\" format=\"ascii\">");
				fprintf(fID,"\t\t\t\t");
				for (i = 0; i < NvnP; i++) {
					fprintf(fID,"% .8e ",s[i]);
					if ((i+1) % 5 == 0 && i != NvnP-1)
						fprintf(fID,"\n\t\t\t\t");
					else if (i == NvnP-1)
						fprintf(fID,"\n");
				}
				fprintf_tn(fID,4,"</DataArray>");

				fprintf_tn(fID,4,"<DataArray type=\"Float32\" Name=\"Mach\" format=\"ascii\">");
				fprintf(fID,"\t\t\t\t");
				for (i = 0; i < NvnP; i++) {
					fprintf(fID,"% .8e ",Mach[i]);
					if ((i+1) % 5 == 0 && i != NvnP-1)
						fprintf(fID,"\n\t\t\t\t");
					else if (i == NvnP-1)
						fprintf(fID,"\n");
				}
				fprintf_tn(fID,4,"</DataArray>");
			} else {
				EXIT_UNSUPPORTED;
			}

			fprintf_tn(fID,3,"</PointData>");


			fprintf_tn(fID,3,"<Cells>");

			VTK_Ncorners = malloc(NE * sizeof * VTK_Ncorners); // free
			for (i = 0; i < NE; i++) {
				VTK_Ncorners[i] = get_ELEMENT_Ncorners(types[i]);
			}

				fprintf_tn(fID,4,"<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">");
				for (i = 0; i < NE; i++) {
					fprintf(fID,"\t\t\t\t");
					for (j = 0, jMax = VTK_Ncorners[i]; j < jMax; j++)
						fprintf(fID,"%6d ",connectivity[i*8+j]);
					fprintf(fID,"\n");
				}
				fprintf_tn(fID,4,"</DataArray>");

				fprintf_tn(fID,4,"<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">");
				fprintf(fID,"\t\t\t\t");
				for (i = 0, sum = 0; i < NE; i++) {
					sum += VTK_Ncorners[i];
					fprintf(fID,"%6d ",sum);
					if ((i+1) % 8 == 0 && i != NE-1)
						fprintf(fID,"\n\t\t\t\t");
					else if (i == NE-1)
						fprintf(fID,"\n");
				}
				fprintf_tn(fID,4,"</DataArray>");

				fprintf_tn(fID,4,"<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">");
				fprintf(fID,"\t\t\t\t");
				for (i = 0; i < NE; i++) {
					fprintf(fID,"%6d ",types[i]);
					if ((i+1) % 8 == 0 && i != NE-1)
						fprintf(fID,"\n\t\t\t\t");
					else if (i == NE-1)
						fprintf(fID,"\n");
				}
				fprintf_tn(fID,4,"</DataArray>");

			fprintf_tn(fID,3,"</Cells>");

		fprintf_tn(fID,2,"</Piece>\n");

		free(XYZ_vP);
		if (strstr(TestCase,"Advection")) {
			free(u);
		} else if (strstr(TestCase,"Poisson")) {
			free(u);
			free(q);
		} else if (strstr(TestCase,"Euler") || strstr(TestCase,"NavierStokes")) {
			free(W_vP);
			free(U_vP);
			free(s);
			free(Mach);
		} else {
			EXIT_UNSUPPORTED;
		}
		free(VTK_Ncorners);
	}

	fprintf_tn(fID,1,"</UnstructuredGrid>");
	fprintf(fID,"</VTKFile>");

	fclose(fID);
}

static void output_mesh_edges(const char *mesh_type)
{
	// Initialize DB Parameters
	char         *TestCase = DB.TestCase,
	             *MeshType = DB.MeshType;
	unsigned int d         = DB.d;
	int          MPIrank   = DB.MPIrank,
	             MPIsize   = DB.MPIsize;

	// standard datatypes
	char         MPIrank_c[STRLEN_MIN], f_name[STRLEN_MAX], f_name_source[STRLEN_MAX],
	             f_parallel[STRLEN_MAX], f_serial[STRLEN_MAX];
	unsigned int i, iMax, j, jMax, dim, sum, P, PP, Ne, NvnP, NvnG, NenP, *connectivityE;
	double       *I_vG_vP, *XYZ_vP;
	FILE         *fID;

	struct S_ELEMENT *ELEMENT;
	struct S_VOLUME *VOLUME;

	sprintf(MPIrank_c,"%d",MPIrank);
	strcpy(f_name,TestCase); strcat(f_name,"/");
	strcat(f_name,MeshType); strcat(f_name,"/");
	strcat(f_name,mesh_type);
	strcpy(f_name_source,mesh_type);

	if (!DB.MPIrank) {
		strcpy(f_parallel,"paraview/");
		strcat(f_parallel,f_name);
		strcat(f_parallel,".pvtu");

		if ((fID = fopen(f_parallel,"w")) == NULL)
			printf("Error: File: %s, did not open.\n",f_parallel), EXIT_MSG;

		fprintf_tn(fID,0,"<?xml version=\"1.0\"?>");
		fprintf_tn(fID,0,"<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">");
		fprintf_tn(fID,1,"<PUnstructuredGrid GhostLevel=\"0\">\n");

		fprintf_tn(fID,2,"<PPoints>");
		fprintf_tn(fID,3,"<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\"/>");
		fprintf_tn(fID,2,"</PPoints>\n");

		fprintf_tn(fID,2,"<PCells>");
		fprintf_tn(fID,3,"<PDataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\"/>");
		fprintf_tn(fID,3,"<PDataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\"/>");
		fprintf_tn(fID,3,"<PDataArray type=\"UInt8\" Name=\"types\" format=\"ascii\"/>");
		fprintf_tn(fID,2,"</PCells>\n");

		for (i = 0, iMax = (unsigned int) MPIsize; i < iMax; i++)
			fprintf(fID,"\t\t<Piece Source=\"%s%d.vtu\"/>\n",f_name_source,i);

		fprintf_tn(fID,1,"</PUnstructuredGrid>");
		fprintf_tn(fID,0,"</VTKFile>");

		fclose(fID);
	}

	strcpy(f_serial,"paraview/");
	strcat(f_serial,f_name);
	strcat(f_serial,MPIrank_c);
	strcat(f_serial,".vtu");

	if ((fID = fopen(f_serial,"w")) == NULL)
		printf("Error: File f_serial did not open.\n"), exit(1);

	fprintf_tn(fID,0,"<?xml version=\"1.0\"?>");
	fprintf_tn(fID,0,"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">");
	fprintf_tn(fID,1,"<UnstructuredGrid>\n");

	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		ELEMENT = get_ELEMENT_type(VOLUME->type);

		P = VOLUME->P;

		// Visualize with higher order to see curved geometry if applicable (P+1 is the maximum currently supported)
		if (!VOLUME->curved || DB.Adapt == ADAPT_0 || DB.Adapt == ADAPT_H)
			PP = P;
		else
			PP = P;
//			PP = min(PMax,P+1);

		connectivityE = ELEMENT->connectivityE[PP];

		Ne   = ELEMENT->Ne;
		NenP = PP+1;
		NvnP = ELEMENT->NvnP[PP];
		NvnG = VOLUME->NvnG;

		if (!VOLUME->curved)
			I_vG_vP = ELEMENT->I_vGs_vP[1][P][0];
		else
			I_vG_vP = ELEMENT->I_vGc_vP[P][PP][0];

		XYZ_vP = mm_Alloc_d(CBCM,CBT,CBNT,NvnP,d,NvnG,1.0,I_vG_vP,VOLUME->XYZ);     // free

		fprintf(fID,"\t\t<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n",NvnP,Ne);

			fprintf_tn(fID,3,"<Points>");
				fprintf(fID,"\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"%d\" format=\"ascii\">\n",3);
				for (i = 0; i < NvnP; i++) {
					fprintf(fID,"\t\t\t\t");
					for (dim = 0; dim < d; dim++)
						fprintf(fID,"% .4e ",XYZ_vP[dim*NvnP+i]);
					for (dim = d; dim < DMAX; dim++)
						fprintf(fID,"% .4e ",0.0);
					fprintf(fID,"\n");
				}
				fprintf_tn(fID,4,"</DataArray>");
			fprintf_tn(fID,3,"</Points>");

			fprintf_tn(fID,3,"<Cells>");

				fprintf_tn(fID,4,"<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">");
				for (i = 0; i < Ne; i++) {
					fprintf(fID,"\t\t\t\t");
					for (j = 0, jMax = NenP; j < jMax; j++)
						fprintf(fID,"%6d ",connectivityE[i*NenP+j]);
					fprintf(fID,"\n");
				}
				fprintf_tn(fID,4,"</DataArray>");

				fprintf_tn(fID,4,"<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">");
				fprintf(fID,"\t\t\t\t");
				for (i = 0, sum = 0; i < Ne; i++) {
					sum += NenP;
					fprintf(fID,"%6d ",sum);
					if ((i+1) % NenP == 0 && i != Ne-1)
						fprintf(fID,"\n\t\t\t\t");
					else if (i == Ne-1)
						fprintf(fID,"\n");
				}
				fprintf_tn(fID,4,"</DataArray>");

				fprintf_tn(fID,4,"<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">");
				fprintf(fID,"\t\t\t\t");
				for (i = 0; i < Ne; i++) {
					fprintf(fID,"%6d ",4); // This is the type for VTK_POLY_LINE
					if ((i+1) % NenP == 0 && i != Ne-1)
						fprintf(fID,"\n\t\t\t\t");
					else if (i == Ne-1)
						fprintf(fID,"\n");
				}
				fprintf_tn(fID,4,"</DataArray>");

			fprintf_tn(fID,3,"</Cells>");

		fprintf_tn(fID,2,"</Piece>\n");

		free(XYZ_vP);
	}

	fprintf_tn(fID,1,"</UnstructuredGrid>");
	fprintf(fID,"</VTKFile>");

	fclose(fID);
}


void output_to_paraview(const char *OutputType)
{
	if (strstr(OutputType,"Geom"))
		output_geom(OutputType);
	else if (strstr(OutputType,"Normals"))
		output_normals(OutputType);
	else if (strstr(OutputType,"Sol"))
		output_solution(OutputType);
	else if (strstr(OutputType,"MeshEdges"))
		output_mesh_edges(OutputType);
	else
		printf("Error: Unsupported OutputType in output_to_paraview.\n"), exit(1);

}
