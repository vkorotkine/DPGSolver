// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "memory_free.h"

#include <stdlib.h>
#include <stdio.h>
 
#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_ELEMENT.h"
#include "S_VOLUME.h"
#include "S_FACE.h"

#include "array_free.h"
#include "memory_destructors.h"
#include "adaptation.h"

/*
 *	Purpose:
 *		Free remaining allocated memory.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

void memory_free(void)
{
	// Project to ML0 if h-adaptation is enabled.
	switch (DB.Adapt) {
	default: // ADAPT_H, ADAPT_HP
		mesh_to_level(0);
		break;
	case ADAPT_0:
	case ADAPT_P:
		// Don't do anything.
		break;
	}

	// DB Parameters

		// Initialization
		free(DB.TestCase);
		free(DB.MeshType);
		free(DB.Geometry);
		free(DB.MeshPath);
		free(DB.BumpOrder);
		free(DB.Form);
		free(DB.NodeType);
		free(DB.BasisType);
		free(DB.MeshFile);

		// Preprocessing

			// setup_parameters
			array_free2_c(NEC,DB.NodeTypeG);
			array_free3_c(DB.NP,NEC,DB.NodeTypeS);
			array_free3_c(DB.NP,NEC,DB.NodeTypeF);
			array_free3_c(DB.NP,NEC,DB.NodeTypeFrs);
			array_free3_c(DB.NP,NEC,DB.NodeTypeFrc);
			array_free3_c(DB.NP,NEC,DB.NodeTypeIfs);
			array_free3_c(DB.NP,NEC,DB.NodeTypeIfc);
			array_free3_c(DB.NP,NEC,DB.NodeTypeIvs);
			array_free3_c(DB.NP,NEC,DB.NodeTypeIvc);

			free(DB.PGc);
			free(DB.PF);
			array_free3_ui(DB.NP,2,DB.SF_BE);
			array_free2_ui(DB.NP,DB.PCs);
			array_free2_ui(DB.NP,DB.PCc);
			array_free2_ui(DB.NP,DB.PJs);
			array_free2_ui(DB.NP,DB.PJc);
			array_free2_ui(DB.NP,DB.PFrs);
			array_free2_ui(DB.NP,DB.PFrc);
			array_free2_ui(DB.NP,DB.PIfs);
			array_free2_ui(DB.NP,DB.PIfc);
			array_free2_ui(DB.NP,DB.PIvs);
			array_free2_ui(DB.NP,DB.PIvc);
			free(DB.VFPartUnity);

			// setup_mesh
			free(DB.PVe), free(DB.NE), free(DB.EType), free(DB.ETags), free(DB.EToVe), free(DB.EToPrt);
			free(DB.VToV), free(DB.VToF), free(DB.VToGF), free(DB.VToBC), free(DB.GFToVe), free(DB.VC), free(DB.GFC);
			free(DB.VeInfo); free(DB.VeXYZ);

			// setup_structures
			free(DB.NVgrp);
			free(DB.Vgrp);

		// Solving

			// Initialize_test_case
			free(DB.SolverType);

	// ELEMENTs
	struct S_ELEMENT *ELEMENT, *ELEMENTnext;

	for (ELEMENT = DB.ELEMENT; ELEMENT; ) {
		ELEMENTnext = ELEMENT->next;
		memory_destructor_E(ELEMENT);
		ELEMENT = ELEMENTnext;
	}

	// VOLUMEs
	struct S_VOLUME *VOLUME, *VOLUMEnext;
	for (VOLUME = DB.VOLUME; VOLUME; ) {
		VOLUMEnext = VOLUME->next;
		if (VOLUME->parent) {
			printf("Error: memory_free (VOLUMEs) requires that the mesh be projected to ML0 before executing.\n");
			EXIT_MSG;
		}
		memory_destructor_V(VOLUME);
		VOLUME = VOLUMEnext;
	}

	// FACEs
	struct S_FACE *FACE, *FACEnext;
	for (FACE = DB.FACE; FACE; ) {
		FACEnext = FACE->next;
		memory_destructor_F(FACE);
		FACE = FACEnext;
	}
}

void memory_free_children(void)
{
	/*
	 *	Purpose:
	 *		Free memory of non-leaf VOLUME and FACE children.
	 */

	// VOLUMEs
	struct S_VOLUME *VOLUME, *VOLUMEc, *VOLUMEcnext;
	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		VOLUMEc = VOLUME->child0;
		while (VOLUMEc) {
			VOLUMEcnext = VOLUMEc->next;
			memory_destructor_V(VOLUMEc);
			VOLUMEc = VOLUMEcnext;
		}
		VOLUME->child0 = NULL;
	}

	// FACEs
	struct S_FACE *FACE, *FACEc, *FACEcnext;
	for (FACE = DB.FACE; FACE; FACE = FACE->next) {
		FACEc = FACE->child0;
		while (FACEc) {
			FACEcnext = FACEc->next;
			memory_destructor_F(FACEc);
			FACEc = FACEcnext;
		}
		FACE->child0 = NULL;
	}
}
