/*
 *	Purpose:
 *		Define global parameters and objects.
 *
 *	Comments:
 *		The notation is presented in the first routine in which parameters appear.
 *		Memory for the S_DB structure is freed in memory_free.c. (ToBeModified)
 *
 *	Notation:
 *
 *	References:
 *
 */

#ifndef DPG__database_h__INCLUDED
#define DPG__database_h__INCLUDED

struct S_DB {
	// MPI and PETSC
	int MPIsize, MPIrank;

	// Initialization
	char *TestCase, *MeshType, *Form, *NodeType, *BasisType, *MeshFile;
	int  d, ML, Vectorized, EFE, Collocated, Adaptive, P, PMax, Restart, Testing;

	// Parameters
	char *Parametrization,
	     ***NodeTypeS,   ***NodeTypeF,   ***NodeTypeFrs, ***NodeTypeFrc,
		 ***NodeTypeIfs, ***NodeTypeIfc, ***NodeTypeIvs, ***NodeTypeIvc;
	int  NP, NDE, AC, ExactGeom, PR, PP, PGs,
	     *PGc, *PF,
	     **SF_BE, **PCs, **PCc, **PJs, **PJc, **PFrs, **PFrc, **PIfs, **PIfc, **PIvs, **PIvc;

	// Mesh
	int    NVe, NPVe, NfMax, NfveMax, NETotal, NV, NGF, NVC, NGFC,
	       *PVe, *NE, *EType, *ETags, *EToVe, *EToPrt, *VToV, *VToF, *VToGF, *VToBC, *GFToVe, *VC, *GFC;
	double *VeXYZ;

	// Structures
	struct S_ELEMENT *ELEMENT;
};
extern struct S_DB DB;

struct S_ELEMENT {
	// Mesh
	int present, type, d, Nve, Nf,
	    *Nfve, *VeC, *VeE, *VeF;

	struct S_ELEMENT *next;
};

#endif // DPG__database_h_INCLUDED
