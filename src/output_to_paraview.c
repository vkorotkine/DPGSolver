#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "database.h"
#include "parameters.h"
#include "functions.h"

/*
 *	Purpose:
 *		Output data to paraview for visualization.
 *
 *	Comments:
 *		Getting an error from paraview when using less than 3 components for the Points DataArray. (2016-04-21)
 *		Potentially place all of the elements in a single piece to reduce file sizes. (ToBeDeleted)
 *		Possibly convert output to binary if this is found to be slow when profiling. (ToBeDeleted)
 *		The following outputs are supported:
 *			G : Geometry only
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

static void output_geom()
{
	// Initialize database parameters
	char         *TestCase = DB.TestCase;
	unsigned int d         = DB.d;
	int          MPIrank   = DB.MPIrank,
	             MPIsize   = DB.MPIsize;

	// standard datatypes
	char MPIrank_c[STRLEN_MIN], f_name[STRLEN_MAX], f_parallel[STRLEN_MAX], f_serial[STRLEN_MAX];
	unsigned int i, iMax, j, jMax, dim, sum,
	             P, NE, NvnP, NvnG, NIn, NOut, NOut_Total, NCols,
				 NIn_SF[3], NOut_SF[3], Diag[3],
	             *connectivity, *types, *VTK_Ncorners;
	double *I_vG_vP, *XYZ_vP, *Input_SF, *OP_SF[3];
	FILE *fID;

	struct S_ELEMENT *ELEMENT, *ELEMENT_class[2];
	struct S_VOLUME *VOLUME;

	sprintf(MPIrank_c,"%d",MPIrank);
	strcpy(f_name,TestCase);
	strcat(f_name,"_geom");

	if (!DB.MPIrank) {
		strcpy(f_parallel,"paraview/");
		strcat(f_parallel,f_name);
		strcat(f_parallel,".pvtu");

		if ((fID = fopen(f_parallel,"w")) == NULL)
			printf("Error: File f_parallel did not open.\n"), exit(1);

		fprintf_tn(fID,0,"<?xml version=\"1.0\"?>");
		fprintf_tn(fID,0,"<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">");
		fprintf_tn(fID,1,"<PUnstructuredGrid GhostLevel=\"0\">\n");

		fprintf_tn(fID,2,"<PPointData Scalars=\"Scalars\" Vectors=\"Vectors\">");
		fprintf_tn(fID,3,"<PDataArray type=\"Float32\" Name=\"dummy\" format=\"ascii\"/>");
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

	/* To Do:
	 * Make a mixed ToBeCurved mesh in gmsh and
	 * run it through the setup_ToBeCurved. Then output to paraview and make sure that all is working correctly. Note
	 * the cubature_ES should also output a vector of element types matching the connectivity (required especially for
	 * TET/PYR).
	 */

	strcpy(f_serial,"paraview/");
	strcat(f_serial,f_name);
	strcat(f_serial,MPIrank_c);
	strcat(f_serial,".vtu");

	if ((fID = fopen(f_serial,"w")) == NULL)
		printf("Error: File f_serial did not open.\n"), exit(1);

	fprintf_tn(fID,0,"<?xml version=\"1.0\"?>");
	fprintf_tn(fID,0,"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">");
	fprintf_tn(fID,1,"<UnstructuredGrid>\n");

	for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next) {
		ELEMENT = get_ELEMENT_type(VOLUME->type);

		connectivity = ELEMENT->connectivity;
		types        = ELEMENT->connect_types;
		NE           = ELEMENT->connect_NE;

		P = VOLUME->P;

		if (VOLUME->Eclass == C_TP) {
			ELEMENT_class[0] = get_ELEMENT_Eclass(ELEMENT->type,0);

			NvnP         = ELEMENT_class[0]->NvnP;

			if (!VOLUME->curved) {
				NvnG = ELEMENT_class[0]->NvnGs[0];
				I_vG_vP = ELEMENT_class[0]->I_vGs_vP[0];
			} else {
				NvnG = ELEMENT_class[0]->NvnGc[P];
				I_vG_vP = ELEMENT_class[0]->I_vGc_vP[P];
			}

			Input_SF = VOLUME->XYZs;

			NIn = NvnG;
			for (dim = 0; dim < 3; dim++) {
				if (dim < d) NIn_SF[dim] = NIn;
				else         NIn_SF[dim] = 1;
			}
			
			NOut = NvnP;
			NOut_Total = 1;
			for (dim = 0; dim < 3; dim++) {
				if (dim < d) NOut_SF[dim] = NOut;
				else         NOut_SF[dim] = 1;
				NOut_Total *= NOut_SF[dim];
			}

			NCols = d*1; // d coordinates * 1 element
			OP_SF[0] = I_vG_vP;
			OP_SF[1] = OP_SF[0];
			OP_SF[2] = OP_SF[0];
			for (dim = 0; dim < 3; dim++) Diag[dim] = 0;

			XYZ_vP = malloc(NOut_Total*NCols * sizeof *XYZ_vP); // free
			sf_apply_d(Input_SF,XYZ_vP,NIn_SF,NOut_SF,NCols,OP_SF,Diag,d);

			NvnP = NOut_Total;
		} else if (VOLUME->Eclass == C_SI || VOLUME->Eclass == C_PYR) {
			NvnP = ELEMENT->NvnP;
			NvnG = VOLUME->NvnG;

			if (!VOLUME->curved)
				I_vG_vP = ELEMENT->I_vGs_vP[0];
			else
				I_vG_vP = ELEMENT->I_vGc_vP[P];

/*			if (VOLUME->indexg == 0)
				array_print_d(NvnP,NvnG,I_vG_vP,'R');

			printf("%d\n",VOLUME->indexg);
			array_print_d(NvnG,d,VOLUME->XYZs,'C');
*/

			XYZ_vP = mm_Alloc_d(CblasColMajor,CblasTrans,CblasNoTrans,NvnP,d,NvnG,1.0,I_vG_vP,VOLUME->XYZs); // free
		} else {
			printf("Add support: output_to_paraview (XYZ_vP)\n");
			exit(1);
		}

//array_print_d(NvnP,d,XYZ_vP,'C');

		fprintf(fID,"\t\t<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n",NvnP,NE);

			fprintf_tn(fID,3,"<Points>");
				fprintf(fID,"\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"%d\" format=\"ascii\">\n",3);
				for (i = 0; i < NvnP; i++) {
					fprintf(fID,"\t\t\t\t");
					for (dim = 0; dim < d; dim++)
						fprintf(fID,"% .3f ",XYZ_vP[dim*NvnP+i]);
					for (dim = d; dim < 3; dim++)
						fprintf(fID,"% .3f ",0.0);
					fprintf(fID,"\n");
				}
				fprintf_tn(fID,4,"</DataArray>");
			fprintf_tn(fID,3,"</Points>");

			fprintf_tn(fID,3,"<PointData Scalars=\"Scalars\" Vectors=\"Vectors\">");
				fprintf_tn(fID,4,"<DataArray type=\"Float32\" Name=\"dummy\" format=\"ascii\">");
				fprintf(fID,"\t\t\t\t");
				for (i = 0; i < NvnP; i++) {
					fprintf(fID,"% .3f ",0.0);
					if ((i+1) % 5 == 0 && i != NvnP-1)
						fprintf(fID,"\n\t\t\t\t");
					else if (i == NvnP-1)
						fprintf(fID,"\n");
				}
				fprintf_tn(fID,4,"</DataArray>");
			fprintf_tn(fID,3,"</PointData>");

			fprintf_tn(fID,3,"<Cells>");

			VTK_Ncorners = malloc(NE * sizeof * VTK_Ncorners); // free
			for (i = 0; i < NE; i++)
				VTK_Ncorners[i] = get_ELEMENT_Ncorners(types[i]);

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

void output_to_paraview(const char OutputType)
{
	if (OutputType == 'G')
		output_geom();

}
