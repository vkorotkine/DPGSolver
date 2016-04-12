#ifndef DPG__database_h__INCLUDED
#define DPG__database_h__INCLUDED

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

struct S_DB {
	// MPI and PETSC
	int MPIsize, MPIrank;

	// Initialization
	char         *TestCase, *MeshType, *Form, *NodeType, *BasisType, *MeshFile;
	unsigned int d, ML, Vectorized, EFE, Collocated, Adaptive, P, PMax, Testing;
	int          Restart;

	// Parameters
	char         *Parametrization,
	             ***NodeTypeS,   ***NodeTypeF,   ***NodeTypeFrs, ***NodeTypeFrc,
	             ***NodeTypeIfs, ***NodeTypeIfc, ***NodeTypeIvs, ***NodeTypeIvc;
	unsigned int NP, NEC, AC, ExactGeom, PR, PP, PGs,
	             *PGc, *PF,
	             **SF_BE, **PCs, **PCc, **PJs, **PJc, **PFrs, **PFrc, **PIfs, **PIfc, **PIvs, **PIvc;

	// Mesh
	unsigned int NVe, NPVe, NfMax, NfveMax, NETotal, NV, NGF, NVC, NGFC,
	             *PVe, *NE, *EType, *ETags, *EToVe, *EToPrt, *VToV, *VToF, *VToGF, *VToBC, *GFToVe, *VC, *GFC;
	double *VeXYZ;

	// Structures
	unsigned int NECgrp;

	// Structs
	struct S_ELEMENT *ELEMENT;
	struct S_VOLUME *VOLUME, **Vgrp;
};
extern struct S_DB DB;

struct S_ELEMENT {
	// Mesh
	unsigned int present, type, d, Nve, Nf,
	             *Nfve, *VeC, *VeE, *VeF;

	// Operators
	unsigned int *NvnGs, *NvnGc, *NvnCs, *NvnCc, *NvnJs, *NvnJc, *NvnS, *NvnF, *NvnFrs, *NvnFrc, *NvnIs, *NvnIc, *NvnP,
	             *NfnGc, *NfnIs, *NfnIc,
	             **Con_rst_vP;
	double       *nr,
	             **rst_vGs, **rst_vGc, **rst_vCs, **rst_vCc, **rst_vJs, **rst_vJc, **rst_vS, **rst_vF, **rst_vFrs, **rst_vFrc,
	             **rst_vIs, **rst_vIc, **rst_vP,
	             **wvIs, **wvIc,
	             ***rst_fGc, ***rst_fIs, ***rst_fIc, **wfIs, **wfIc,
	             **I_vGs_vGc;

	struct S_ELEMENT *next;
};

struct S_VOLUME {
	// Structures
	unsigned int indexl, indexg, P, type, Eclass, curved;

	// Geometry
	double *XYZc, *XYZs;

	// structs
	struct S_VOLUME *next, *grpnext;

};

#endif // DPG__database_h__INCLUDED
