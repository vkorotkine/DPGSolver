// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "setup_Bezier.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "Macros.h"
#include "Parameters.h"
#include "S_DB.h"
#include "S_VOLUME.h"

#include "bases.h"
#include "cubature.h"
#include "select_functions.h"
#include "array_print.h"
#include "matrix_functions.h"

static void setup_XYZ_volume(struct S_VOLUME *VOLUME, double *ControlPoints, int nControlPoints,
		int *volumeConnectivity, int volIndex, double *ChiRefG, int NBF){

	/*
	Function used to compute the XYZ values for the given volume. That is,
	using the XYZ_hat vector, compute the XYZ vector (value at the geometry
	nodes). This is done because the rest of the code relies on the geometry to
	stored in this form when it comes to computing the metric terms, ...

	Input:
		- struct S_VOLUME *VOLUME = The pointer to the volume struct to build the XYZ 
			vector for.
		- ControlPoints = The array with all the control points for the mesh (stored in 
			column major form).
		- nControlPoints = The number of control points
		- volumeConnectivity = The array with the control point connectivity information
			for each volume. 
		- volIndex = The index of the volume in the linked list of volumes
		- ChiRefG = The operator which contains the basis functions evaluated at the geometry 
			node point positions (on the reference element). Multiply this to the XYZ_hat for
			each element to get the geometry node points.
		- NBF = Number of basis functions

	Output:
		None

	*/

	double *XYZ_hat, *XYZ;

	unsigned int NvnG, d, i, j, controlPointIndex;

	// NOTE: NvnG is the same for the number of control points for the given 
	//	volume.
	NvnG = VOLUME->NvnG;
	d = DB.d;

	// Build the XYZ_hat matrix
	XYZ_hat = malloc(NvnG*d * sizeof *XYZ_hat);  // free

	for (i = 0; i < NvnG; i++){
		// Loop over the control point connectivity for the given
		// volume index.

		// Get the index for which control point we are interested in the large list
		// of control points (subtract 1 since we have 0 based indexing for the array but
		// the connectivity file uses indexing starting from 1)
		controlPointIndex = volumeConnectivity[volIndex*NvnG + i] - 1;

		for (j = 0; j < d; j++){
			// Loop over the dimensions (x,y,z) and load the given control point's coordinates
			XYZ_hat[j*NvnG + i] = ControlPoints[j*nControlPoints + controlPointIndex];
		}
	}

	// There are an equal number of basis functions and geometry node points.
	XYZ = mm_Alloc_d(CBCM, CBT, CBNT, NvnG, d, NBF, 1.0, ChiRefG, XYZ_hat);
	VOLUME->XYZ = XYZ;


	if (!DB.MPIrank && !DB.Testing){
		// Compare the control points to the geometry node points location for each
		// volume

		printf("----------------------------------------------------- \n");
		printf("VOLUME : \n");

		printf("ChiRefG : \n");
		array_print_d(NvnG, NBF, ChiRefG, 'C');

		printf("XYZ_HAT : \n");
		array_print_d(NvnG, d, XYZ_hat, 'C');

		printf("XYZ : \n");
		array_print_d(NvnG, d, VOLUME->XYZ, 'C');

		printf("----------------------------------------------------- \n");

	}

	free(XYZ_hat);

}	


void setup_Bezier(){

	/*
	Function for setting up the Bezier elements using the Bezier extracted
	B-spline mesh. This function will do the following:
		1) Read the msh file to load in the control points and control 
			point connectivity for each element.
		2) Build the Chi_xyz operator which, when multiplied to the XYZ_hat vector, 
			will yield the XYZ vector for each element.
		3) Using the two sources of information, build the XYZ_hat vector
			for each element (no need to store it since we can get this 
			vector by using the inverse Chi operators on the XYZ vector).
			- Multiply the operator to this and store the XYZ values for the element.

	Input:
		None

	Output:
		None
	*/

	// TODO: Modify everywhere we have PGlobal to instead use PGc[PGlobal].
	//		We need to be able to handle the super/sub case.

	char *MeshFile = DB.MeshFile;
	unsigned int d = DB.d;

	unsigned int numControlPoints, numControlPointsPerVol, i, j;
	double *Bezier_XYZ_hat;
	int *BezierElem_Connectivity;
	char StringRead[STRLEN_MAX], *strings, *stringe;
	FILE *fID;

	long tmpl;
	double tmpd;

	// 1) Read the msh file to get the Bezier element information
	if ((fID = fopen(MeshFile,"r")) == NULL)
		printf("Mesh file: %s not present.\n",MeshFile), EXIT_MSG;

	// - Find numControlPoints. If it doesn't exist, then there is an issue
	numControlPoints = 0;
	while (fscanf(fID,"%[^\n]\n",StringRead) == 1) {
		if (strstr(StringRead,"$ControlPoints$")) {
			if (fscanf(fID,"%[^\n]\n",StringRead) == 1)
				sscanf(StringRead,"%d",&numControlPoints);
		}
	}
	rewind(fID);

	if (numControlPoints == 0){
		printf("BEZIER MESH FILE HAS NO CONTROL POINTS \n");
		EXIT_MSG;
	}

	// - Allocate the arrays to hold the Bezier element information (temporary)
	Bezier_XYZ_hat = malloc(numControlPoints*d * sizeof *Bezier_XYZ_hat);  // free
	// 		Each Bezier element (there are NV volumes) has (P+1) ^ d control points associated
	//		with it.
	numControlPointsPerVol = pow(DB.PGc[DB.PGlobal]+1, d);
	BezierElem_Connectivity = malloc(DB.NV*numControlPointsPerVol * sizeof *BezierElem_Connectivity);  // free

	// Read the file now and load the data into each array
	while (fscanf(fID,"%[^\n]\n",StringRead) == 1) {

		// Read the control point information
		if (strstr(StringRead,"$ControlPoints$")) {

			// First line is just the number of control points
			fscanf(fID,"%[^\n]\n",StringRead);

			for (i = 0; i < numControlPoints; i++){
				// Loop over all the control points (i index)

				if(fscanf(fID,"%[^\n]\n",StringRead) == 1){
					strings = StringRead;

					// Control Point Index
					tmpl = strtol(strings,&stringe,10); 
					strings = stringe;

					// Control Points:
					for (j = 0; j < d; j++){
						// For each ith control point, loop through the 
						// dimensions d.
						tmpd = strtod(strings, &stringe);
						strings = stringe;

						// Column major ordering of nodes
						Bezier_XYZ_hat[numControlPoints*j + i] = tmpd;
					}

				} else{
					printf("NOT ENOUGH CONTROL POINTS \n");
					EXIT_MSG;
				}
			}
		}

		// Read the connectivity information for the control points
		// and the elements
		if (strstr(StringRead,"$ControlPointConnectivity")) {

			for (i = 0; i < DB.NV; i++){
				// Loop over all the volumes

				if(fscanf(fID,"%[^\n]\n",StringRead) == 1){
					strings = StringRead;

					for (j = 0; j < numControlPointsPerVol; j++){
						// Loop over the number of control points for each volume
						tmpd = strtod(strings, &stringe);
						strings = stringe;

						// Place the (P+1)^d control points for each element in order
						// in the array (one after another that is)
						BezierElem_Connectivity[i*numControlPointsPerVol + j] = tmpd;
					}

				} else{
					printf("NOT ENOUGH CONNECTIVITY INFORMATION \n");
					EXIT_MSG;
				}
			}
		}
	}

	if (!DB.MPIrank && !DB.Testing) {
		// Now with all the bezier information loaded, print it out to verify
		// that it was loaded correctly
		printf("\n CONTROL POINTS : \n");
		for (i = 0; i < numControlPoints; i++){
			printf("%d %.15e %.15e \n", i+1, Bezier_XYZ_hat[i], Bezier_XYZ_hat[numControlPoints + i]);
		}

		printf("\n CONNECTIVITY CONTROL POINTS \n");
		for (i = 0; i < DB.NV; i++){
			for (j = 0; j < numControlPointsPerVol; j++){
				printf("%d ", BezierElem_Connectivity[i*numControlPointsPerVol + j]);
			}
			printf("\n");
		}
	}

	// 2) Build the Chi_Geo operator which is the Chi operator evaluated at the given 
	//	geometry node points. 

	// Basis functions (Chi) evaulated at the geometry node points
	// on the reference space.
	double *ChiRefG;  
	unsigned int NBFOut;
	
	// Function pointers
	cubature_tdef   cubature;
	basis_tdef      basis;
	grad_basis_tdef grad_basis;

	select_functions(&basis,&grad_basis,&cubature,QUAD);

	struct S_CUBATURE *CUBDATA = malloc(sizeof *CUBDATA); // free
	set_cubdata(CUBDATA,false,false,DB.NodeTypeG[0],DB.d,DB.PGc[DB.PGlobal],cubature); // free
	ChiRefG = basis(DB.PGc[DB.PGlobal], CUBDATA->rst, pow((DB.PGc[DB.PGlobal]+1), DB.d), &NBFOut, DB.d);

	// NOTE: ChiRefG is of dimension NBF x Nn. Take the transpose before multiplying it
	//	to any matrix with the coefficients (XYZ_hat).

	if (!DB.MPIrank && !DB.Testing){
		printf("ChiRefG : \n");
		array_print_d(pow((DB.PGc[DB.PGlobal]+1), DB.d), NBFOut, ChiRefG, 'C');
	}	

	// 3) Using the ChiRefG operator, compute the XYZ coordinates for each volume
	//	using their XYZ_hat values (from the control point list and the connectivity
	//	lists for each element).

	int volIndex = 0;
	struct S_VOLUME *VOLUME;
	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next){
		setup_XYZ_volume(VOLUME, Bezier_XYZ_hat, numControlPoints, BezierElem_Connectivity, volIndex++, ChiRefG, NBFOut);
	}

	// Free temporarily allocated arrays
	free(CUBDATA);

}





