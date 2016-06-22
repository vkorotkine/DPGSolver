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
	// Time
	double time_total;

	// MPI and PETSC
	int MPIsize, MPIrank;

	// Initialization
	char         *TestCase, *MeshType, *Form, *NodeType, *BasisType, *MeshFile;
	unsigned int d, ML, Vectorized, EFE, Collocated, Adapt, PGlobal, PMax, Testing, *BumpOrder;
	int          Restart;

	// Parameters
	char         *Parametrization,
	             **NodeTypeG,
	             ***NodeTypeS,   ***NodeTypeF,   ***NodeTypeFrs, ***NodeTypeFrc,
	             ***NodeTypeIfs, ***NodeTypeIfc, ***NodeTypeIvs, ***NodeTypeIvc;
	unsigned int NP, NEC, AC, ExactGeom, InviscidFluxType, ExplicitSolverType, PR, PP, PGs,
	             *PGc, *PF,
	             ***SF_BE, **PCs, **PCc, **PJs, **PJc, **PFrs, **PFrc, **PIfs, **PIfc, **PIvs, **PIvc;

	// Mesh
	unsigned int NVe, NPVe, NfMax, NfveMax, NveMax, NfrefMax, NETotal, NV, NGF, NVC, NGFC,
	             *PVe, *NE, *EType, *ETags, *EToVe, *EToPrt, *VToV, *VToF, *VToGF, *VToBC, *GFToVe, *VC, *GFC;
	double *VeXYZ;

	// Structures
	unsigned int NECgrp;

	// Initialization
	char         *SolverType;
	unsigned int Nvar, Neq, OutputInterval, DOF0;
	double       Xc, Yc, Rc, MInf, pInf, TInf, VInf, uInf, vInf, wInf, Rg, Cscale, PeriodL, PeriodFraction, FinalTime,
	             rIn, MIn, rhoIn, VIn;

	// Vectorization
	unsigned int update, NTVgrp, *NVgrp;

	// hp adaptation
	double DOFcap_frac, refine_frac, coarse_frac;

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
	unsigned int connect_NE, NvnP, *Nvve,
	             *NvnGs, *NvnGc, *NvnCs, *NvnCc, *NvnIs, *NvnIc, *NvnS,
	             **NfnIs, **NfnIc,
	             *Nfref, *NfMixed,
	             *connectivity, *connect_types,
	             ***nOrd_fIs, ***nOrd_fIc;
	double       **VeF, **VeV, *nr,
	             **w_vIs, **w_vIc,
	             ****ChiS_vP, ****ChiS_vIs, ****ChiS_vIc,
	             ****ChiInvS_vS,
	             ****ICs, ****ICc,
	             ****I_vGs_vP, ****I_vGs_vGc, ****I_vGs_vCs, ****I_vGs_vIs, ****I_vGs_vIc, ****I_vGs_vS,
	             ****I_vGc_vP,                ****I_vGc_vCc, ****I_vGc_vIs, ****I_vGc_vIc, ****I_vGc_vS,
	             ****I_vCs_vIs, ****I_vCs_vIc,
	             ****I_vCc_vIs, ****I_vCc_vIc,
	             *****D_vGs_vCs, *****D_vGs_vIs,
	             *****D_vGc_vCc, *****D_vGc_vIc,
	             *****D_vCs_vCs,
	             *****D_vCc_vCc,
	             ****ChiS_fIs, ****ChiS_fIc,
	             ****I_vGs_fIs, ****I_vGs_fIc,
	             ****I_vGc_fIs, ****I_vGc_fIc,
	             ****I_vCs_fIs, ****I_vCs_fIc,
	             ****I_vCc_fIs, ****I_vCc_fIc,
	             ****Is_Weak_VV, ****Ic_Weak_VV,
	             ****Is_Weak_FF, ****Ic_Weak_FF,
	             *****Ds_Weak_VV, *****Dc_Weak_VV;

	struct S_ELEMENT *next;
	struct S_ELEMENT **ELEMENTclass, **ELEMENT_FACET;
};

struct S_VOLUME {
	// Structures
	unsigned int indexl, indexg, P, type, Eclass, update, curved,
	             *Vneigh, *Fneigh;
	double *XYZ_vC;

	// Geometry
	unsigned int NvnG;
	double *XYZ_S, *XYZ, *detJV_vI, *C_vC, *C_vI, **C_vf;

	// Initialization
	unsigned int NvnS;
	double *What, *RES;

	// Solving
	double *RHS, *wdetJV_vI, *MInv;

	// hp adaptivity
	unsigned int Vadapt, adapt_type, PNew;
//	double       minRES, maxRES;

	// structs
	struct S_VOLUME *next, *grpnext;

};

struct S_FACET {
	// Structures
	unsigned int P, type, VfIn, VfOut, indexg, BC, IndOrdInOut, IndOrdOutIn;

	// Geometry
	char   curved, typeInt;
	double *XYZ_fI, *n_fI, *detJF_fI;

	// Solving
	double *RHSIn, *RHSOut;

	// structs
	struct S_VOLUME *VIn, *VOut;
	struct S_FACET  *next, *grpnext;
};

struct S_OpCSR {
	unsigned int NRows, NVals, *rowIndex, *columns;
	double       *values;
};

#endif // DPG__database_h__INCLUDED
