#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "database.h"
#include "parameters.h"
//#include "functions.h"

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
	unsigned int i, iMax;

	FILE *fID;

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
	 * Change cubature_TP to no longer output ES nodes, instead, make a cubature_ES script which takes care of all
	 * element types. It seems convenient to use the barycentric coordinates to find the nodes and it will thus be
	 * required to implement basisTET and basisPYR as well. After this is done, make a mixed ToBeCurved mesh in gmsh and
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

	fprintf(fID, "Entered serial\n");

	fclose(fID);

}

void output_to_paraview(const char OutputType)
{
	if (OutputType == 'G')
		output_geom();

}
