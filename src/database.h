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
	unsigned int d, ML, Vectorized, EFE, Collocated, Adapt, PGlobal, PMax, LevelsMax, Testing, *BumpOrder;
	int          Restart;

	// Parameters
	char         *Parametrization,
	             **NodeTypeG,
	             ***NodeTypeS,   ***NodeTypeF,   ***NodeTypeFrs, ***NodeTypeFrc,
	             ***NodeTypeIfs, ***NodeTypeIfc, ***NodeTypeIvs, ***NodeTypeIvc;
	unsigned int NP, NEC, AC, ExactGeom, InviscidFluxType, ExplicitSolverType, PR, PP, PGs,
	             *PGc, *PF, *VFPartUnity,
	             ***SF_BE, **PCs, **PCc, **PJs, **PJc, **PFrs, **PFrc, **PIfs, **PIfc, **PIvs, **PIvc;

	// Mesh
	unsigned int NVe, NPVe, NfMax, NfveMax, NveMax, NfrefMax, NETotal, NV, NVglobal, NGF, NVC, NGFC,
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
	unsigned int present, type, d, Nve, Nf, Nvref, NvrefSF, Eclass,
	             *Nfve, *VeCGmsh, *VeFcon, *NrefV;

	// Operators
	unsigned int *connect_NE, *NvnP, *Nvve,
	             *NvnGs, *NvnGc, *NvnCs, *NvnCc, *NvnIs, *NvnIc, *NvnS,
	             **NfnS, **NfnIs, **NfnIc,
	             *Nfref, *NfMixed,
	             **connectivity, **connect_types,
	             ***nOrd_fS, ***nOrd_fIs, ***nOrd_fIc;
	double       **VeF, **VeV, *nr,
	             **w_vIs, **w_vIc,
	             ****ChiS_vP, ****ChiS_vS, ****ChiS_vIs, ****ChiS_vIc,
	             ****ChiInvS_vS,
	             ****ICs, ****ICc,
	             ****I_vGs_vP, ****I_vGs_vGs, ****I_vGs_vGc, ****I_vGs_vCs, ****I_vGs_vS, ****I_vGs_vIs, ****I_vGs_vIc,
	             ****I_vGc_vP,                               ****I_vGc_vCc, ****I_vGc_vS, ****I_vGc_vIs, ****I_vGc_vIc,
	             ****I_vCs_vS, ****I_vCs_vIs, ****I_vCs_vIc,
	             ****I_vCc_vS, ****I_vCc_vIs, ****I_vCc_vIc,
	             ****Ihat_vS_vS,
	             *****D_vGs_vCs, *****D_vGs_vIs,
	             *****D_vGc_vCc, *****D_vGc_vIc,
	             *****D_vCs_vCs,
	             *****D_vCc_vCc,
	             ****ChiS_fS, ****ChiS_fIs, ****ChiS_fIc,
	             ****I_vGs_fS, ****I_vGs_fIs, ****I_vGs_fIc,
	             ****I_vGc_fS, ****I_vGc_fIs, ****I_vGc_fIc,
	             ****I_vCs_fS, ****I_vCs_fIs, ****I_vCs_fIc,
	             ****I_vCc_fS, ****I_vCc_fIs, ****I_vCc_fIc,
	             ****Is_Weak_VV, ****Ic_Weak_VV,
	             ****Is_Weak_FF, ****Ic_Weak_FF,
	             *****Ds_Weak_VV, *****Dc_Weak_VV,
	             ****L2hat_vS_vS,
	             ****GfS_fIs, ****GfS_fIc;

	struct S_OpCSR ****ChiS_fIs_sp, ****ChiS_fIc_sp,
	               *****Ds_Weak_VV_sp, *****Dc_Weak_VV_sp,
	               ****Is_Weak_FF_sp, ****Ic_Weak_FF_sp;

	struct S_ELEMENT *next;
	struct S_ELEMENT **ELEMENTclass, **ELEMENT_FACET;
};

struct S_VOLUME {
	// Structures
	unsigned int indexl, indexg, P, type, Eclass, update, curved, level,
	             *NsubF, *neigh, *neigh_f;
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
	unsigned int refine_current, Vadapt, adapt_type, PNew, hrefine_type;
//	double       minRES, maxRES;

	// structs
	struct S_VOLUME *next, *grpnext, *child0, *parent;
	struct S_FACET  **FACET;
};

struct S_FACET {
	// Structures
	unsigned int P, type, VfIn, VfOut, indexg, BC, IndOrdInOut, IndOrdOutIn, level, update, adapt_type;

	// Geometry
	unsigned int curved, typeInt;
	double       *XYZ_fI, *XYZ_fS, *n_fI, *n_fS, *detJF_fI, *detJF_fS;

	// Solving
	double *RHSIn, *RHSOut;

	// structs
	struct S_VOLUME *VIn, *VOut;
	struct S_FACET  *next, *grpnext, *child0, *parent;
};

struct S_OpCSR {
	unsigned int NRows, NVals, *rowIndex, *columns;
	double       *values;
};

#endif // DPG__database_h__INCLUDED
