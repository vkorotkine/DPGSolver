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
 *		The following outputs are supported:
 *			G : Geometry only
 *
 *	Notation:
 *
 *	References:
 */

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
	char *TestCase = DB.TestCase;
	int  MPIrank   = DB.MPIrank,
	     MPIsize   = DB.MPIsize;

	// standard datatypes
	char MPIrank_c[STRLEN_MIN], f_name[STRLEN_MAX], f_parallel[STRLEN_MAX], f_serial[STRLEN_MAX];
	unsigned int i, iMax,
	             d, P, NE, NvnP, NvnG,
	             *connectivity, *types;
	double *I_vG_vP, *XYZ_vP;
	FILE *fID;

	struct S_ELEMENT *ELEMENT;
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

/*
		fprintf_tn(fID,2,"<PPointData Scalars=\"Scalars\" Vectors=\"Vectors\">");
		fprintf_tn(fID,3,"<PDataArray type=\"Float32\" Name=\"Temperature\" format=\"ascii\"/>");
		fprintf_tn(fID,3,"<PDataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"ascii\"/>");
		fprintf_tn(fID,2,"</PPointData>\n");
*/

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
		P = VOLUME->P;
		NvnG = VOLUME->NvnG;
		ELEMENT = get_ELEMENT_type(VOLUME->type);

		d = ELEMENT->d;

		connectivity = ELEMENT->connectivity;
		types        = ELEMENT->connect_types;
		NvnP         = ELEMENT->NvnP;
		NE           = ELEMENT->connect_NE;

		if (!VOLUME->curved)
			I_vG_vP = ELEMENT->I_vGs_vP[0];
		else
			I_vG_vP = ELEMENT->I_vGc_vP[P];

		array_print_ui(NE,8,connectivity,'R');

		//array_print_d(NvnG,d,VOLUME->XYZs,'C');
		//array_print_d(NvnP,NvnG,I_vG_vP,'R');

		printf("%d %d %d\n",NvnP,d,NvnG);
		exit(1);

		XYZ_vP = mm_Alloc_d(CblasColMajor,CblasTrans,CblasNoTrans,NvnP,d,NvnG,1.0,I_vG_vP,VOLUME->XYZs); // free

		array_print_d(NvnG,d,VOLUME->XYZs,'C');
		array_print_d(NvnP,d,XYZ_vP,'C');
		exit(1);

		fprintf(fID,"\t\t<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n",1,2);

		fprintf_tn(fID,2,"</Piece>\n");

		free(XYZ_vP);
	}


	fclose(fID);

}

void output_to_paraview(const char OutputType)
{
	if (OutputType == 'G')
		output_geom();

}
