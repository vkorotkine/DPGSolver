// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

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
	unsigned int d, ML, Vectorized, EFE, Collocated, Adaptive, P, PMax, Testing, *BumpOrder;
	int          Restart;

	// Parameters
	char         *Parametrization,
	             **NodeTypeG,
	             ***NodeTypeS,   ***NodeTypeF,   ***NodeTypeFrs, ***NodeTypeFrc,
	             ***NodeTypeIfs, ***NodeTypeIfc, ***NodeTypeIvs, ***NodeTypeIvc;
	unsigned int NP, NEC, AC, ExactGeom, PR, PP, PGs,
	             *PGc, *PF,
	             **SF_BE, **PCs, **PCc, **PJs, **PJc, **PFrs, **PFrc, **PIfs, **PIfc, **PIvs, **PIvc;

	// Mesh
	unsigned int NVe, NPVe, NfMax, NfveMax, NveMax, NfrefMax, NETotal, NV, NGF, NVC, NGFC,
	             *PVe, *NE, *EType, *ETags, *EToVe, *EToPrt, *VToV, *VToF, *VToGF, *VToBC, *GFToVe, *VC, *GFC;
	double *VeXYZ;

	// Structures
	unsigned int NECgrp;

	// Initialization
	unsigned int Nvar, Neq;

	// Structs
	struct S_ELEMENT *ELEMENT;
	struct S_VOLUME  *VOLUME, **Vgrp;
	struct S_FACET   *FACET;
};
extern struct S_DB DB;

struct S_ELEMENT {
	// Mesh
	unsigned int present, type, d, Nve, Nf,
	             *Nfve, *VeCGmsh, *VeFcon;

	// Operators
	unsigned int connect_NE, NvnP,
	             *NvnGs, *NvnGc, *NvnCs, *NvnCc, *NvnIs, *NvnIc,
	             **NfnIs, **NfnIc,
	             *Nfref, *NfMixed,
	             *connectivity, *connect_types;
	double       *VeF, *nr,
	//             **rst_vGs, **rst_vGc, **rst_vCs, **rst_vCc, **rst_vJs, **rst_vJc, **rst_vS, **rst_vF, **rst_vFrs, **rst_vFrc,
	//             **rst_vIs, **rst_vIc,
	//             **wvIs, **wvIc,
	//             ***rst_fGc, ***rst_fIs, ***rst_fIc, **wfIs, **wfIc,
	             **ICs, **ICc,
	             **I_vGs_vP, **I_vGs_vGc, **I_vGs_vCs, **I_vGs_vIs, **I_vGs_vIc, ****I_vGs_fIs, ****I_vGs_fIc,
	             **I_vGc_vP, **I_vGc_vCc,              **I_vGc_vIs, **I_vGc_vIc, ****I_vGc_fIs, ****I_vGc_fIc,
	             **I_vCs_vIs, **I_vCs_vIc,
	             **I_vCc_vIs, **I_vCc_vIc,
	             ****I_vCs_fIs, ****I_vCs_fIc,
	             ****I_vCc_fIs, ****I_vCc_fIc,
	             ***D_vGs_vCs, ***D_vGs_vIs,
	             ***D_vGc_vCc, ***D_vGc_vIc,
	             ***D_vCs_vCs,
	             ***D_vCc_vCc;

	struct S_ELEMENT *next;
	struct S_ELEMENT **ELEMENTclass;
};

struct S_VOLUME {
	// Structures
	unsigned int indexl, indexg, P, type, Eclass, curved, 
	             *Vneigh, *Fneigh;
	double *XYZ_vC;

	// Geometry
	unsigned int NvnG;
	double *XYZ_S, *XYZ, *detJV_vI, *C_vC, *C_vI, **C_vf;

	// structs
	struct S_VOLUME *next, *grpnext;

};

struct S_FACET {
	// Structures
	unsigned int P, VfIn, VfOut, indexg;

	// Geometry
	char   curved, typeInt;
	double *n, *nr;

	// structs
	struct S_VOLUME *VIn, *VOut;
	struct S_FACET  *next, *grpnext;
};

#endif // DPG__database_h__INCLUDED
