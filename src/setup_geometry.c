#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "database.h"
#include "parameters.h"
#include "functions.h"

//#include "petscsys.h"

/*
 *	Purpose:
 *		Set up geometry.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
*/

struct S_ELEMENT *get_ELEMENT_type(const int type);

void setup_geometry()
{
	// Initialize DB Parameters
	int  ExactGeom = DB.ExactGeom,
	     d         = DB.d,
	     NV        = DB.NV,
		 *NE       = DB.NE,
		 *EToVe    = DB.EToVe,

	     Testing   = DB.Testing;

	double *VeXYZ  = DB.VeXYZ;

	int  PrintTesting = 0, MPIrank = DB.MPIrank;

	// Standard datatypes
	int i, ve, dim, v, P,
	    dE, Nve, Vs, PMax,
		*VeC;
	double *XYZc;

	struct S_ELEMENT *ELEMENT;
	struct S_VOLUME  *VOLUME;

	Vs = 0; for (i = 0; i < d; i++) Vs += NE[i];

	// Modify vertex locations if exact geometry is known
	if (ExactGeom) {
		if(!MPIrank) printf("    Modify vertex nodes if exact geometry is known\n");
		printf("Did not yet verify the implementation.\n");
		vertices_to_exact_geom();
	}

	// Set up global XYZ VOLUME coordinates at (s)tart (i.e. before curving)
	double **XYZs;

	VOLUME = DB.VOLUME;
	v = 0;
	while (VOLUME != NULL) {
		ELEMENT = get_ELEMENT_type(VOLUME->type);

		dE  = ELEMENT->d;
		Nve = ELEMENT->Nve;
		VeC = ELEMENT->VeC;

		XYZc = malloc (Nve*dE * sizeof *XYZc); // keep/free (tbd)

		for (ve = 0; ve < Nve; ve++) {
		for (dim = 0; dim < dE; dim++) {
			XYZc[ve*dE+dim] = VeXYZ[EToVe[(Vs+v)*8+VeC[ve]]*d+dim];
		}}

		XYZs = VOLUME->XYZs;
		for (P = 0; P <= PMax; P++) {
			if (!VOLUME->curved) {
				XYZs[P] = XYZc; // keep XYZc
			} else {
				// get operators

				if (VOLUME->Eclass == C_TP) {
					// write the sum factorization code

				} else if (VOLUME->Eclass == C_SI) {
					// write a routine to return pointer to element of correct type

					// figure out the matrix multiplication
					//XYZs = 
				}
			}
		}

		v++;
		VOLUME = VOLUME->next;
	}
}

struct S_ELEMENT *get_ELEMENT_type(const int type)
{
	struct S_ELEMENT *ELEMENT = DB.ELEMENT;
	
	while (ELEMENT != NULL) {
		if (type == ELEMENT->type)
			return ELEMENT;

		ELEMENT = ELEMENT->next;
	}

	printf("Error: Element type not found.\n"), exit(1);
}
