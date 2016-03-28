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
	int  NP, NEC, AC, ExactGeom, PR, PP, PGs,
	     *PGc, *PF,
	     **SF_BE, **PCs, **PCc, **PJs, **PJc, **PFrs, **PFrc, **PIfs, **PIfc, **PIvs, **PIvc;

	// Mesh
	int    NVe, NPVe, NfMax, NfveMax, NETotal, NV, NGF, NVC, NGFC,
	       *PVe, *NE, *EType, *ETags, *EToVe, *EToPrt, *VToV, *VToF, *VToGF, *VToBC, *GFToVe, *VC, *GFC;
	double *VeXYZ;

	// Structures
	struct S_ELEMENT *ELEMENT;
	struct S_VOLUME *VOLUME;
};
extern struct S_DB DB;

struct S_ELEMENT {
	// Mesh
	int present, type, d, Nve, Nf,
	    *Nfve, *VeC, *VeE, *VeF;

	// Operators
	int    *NvnGs, *NvnGc, *NvnCs, *NvnCc, *NvnJs, *NvnJc, *NvnS, *NvnF, *NvnFrs, *NvnFrc, *NvnIs, *NvnIc, *NvnP,
	       *NfnGc, *NfnIs, *NfnIc,
	       **Con_xir_vP;
	double *nr,
	       **xir_vGs, **xir_vGc, **xir_vCs, **xir_vCc, **xir_vJs, **xir_vJc, **xir_vS, **xir_vF, **xir_vFrs, **xir_vFrc,
	       **xir_vIs, **xir_vIc, **xir_vP,
	       **WvIs, **WvIc,
	       ***xir_fGc, ***xir_fIs, ***xir_fIc, **WfIs, **WfIc,
		   **I_vGs_vGc;

	struct S_ELEMENT *next;
};

struct S_VOLUME {
	// Structures
	int P, type, Eclass, curved;

	// Geometry
	double *XYZc, *XYZs;

	struct S_VOLUME *next;

};

#endif // DPG__database_h_INCLUDED
