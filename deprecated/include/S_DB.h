// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__S_DB_h__INCLUDED
#define DPG__S_DB_h__INCLUDED

#include <stdbool.h>

struct S_DB {
	// Time
	double time_total;

	// MPI and PETSC
	int MPIsize, MPIrank;

	// Initialization
	char         *TestCase, *MeshType, *MeshPath, *Form, *NodeType, *BasisType, *MeshFile, *Geometry, *GeomSpecifier,
	             *PDE, *PDESpecifier;
	unsigned int d, ML, Vectorized, EFE, Collocated, Adapt, PGlobal, PMax, LevelsMax, Testing, *BumpOrder, Method;
	unsigned int BezierBasis, BezierMesh;
	int          Restart;

	// Parameters
	bool         AllowSparseVOL, AllowSparseFACE;
	char         **NodeTypeG,
	             ***NodeTypeS,   ***NodeTypeF,   ***NodeTypeFrs, ***NodeTypeFrc,
	             ***NodeTypeIfs, ***NodeTypeIfc, ***NodeTypeIvs, ***NodeTypeIvc;
	unsigned int NP, AC, ExactGeom, InviscidFluxType, ViscousFluxType, ExplicitSolverType, PR, PP, PGs,
	             *PGc, *PF, *VFPartUnity, Blending, Blending_HO, Parametrization,
	             ***SF_BE, **PCs, **PCc, **PJs, **PJc, **PFrs, **PFrc, **PIfs, **PIfc, **PIvs, **PIvc, *PTRS, *PTRF;

	// Mesh
	unsigned int NVe, NPVe, NfMax, NfveMax, NveMax, NfrefMax, NETotal, NV, NVglobal, NGF, NVC, NGFC,
	             *PVe, *NE, *EType, *ETags, *EToVe, *EToPrt, *VToV, *VToF, *VToGF, *VToBC, *GFToVe, *VC, *GFC, *VeInfo;
	double *VeXYZ;

	// Structures
	unsigned int NECgrp;

	// Initialization
	bool         Inviscid, Viscous, Const_mu, Symmetric, Fv_func_of_W;
	char         *SolverType, imex_type;
	unsigned int Nvar, Neq, OutputInterval, DOF0, SourcePresent, PDE_index, MMS;
	double       Xc, Yc, Rc, MInf, rhoInf, pInf, TInf, cInf, VInf, uInf, vInf, wInf, Rg, Cscale, pBack, p_Total, T_Total,
	             PeriodL, PeriodFraction, FinalTime, Poisson_scale, ADV_b[3], ADV_XYZB[3],
	             rIn, MIn, rhoIn, pIn, VIn, rOut, TIn, uIn, omega, mu, Pr, Cp, kappa,
	             Q0, KMin, KMax, GBa, GBb, GBc, aIn, bIn, cIn, aOut, bOut, cOut, lHR,
	             NSc, NSt, NS0, NS1, NS2, NS3, NS4, JSa, JSl, JSxL, rho_store[7], p_store[7], u_store[7], v_store[7],
	             geo_store[4], f_store[7];

	// Solving
	unsigned int dof;
	bool         init_with_L2;

	// Vectorization
	unsigned int update, NTVgrp, *NVgrp;

	// hp adaptation
	unsigned int TETrefineType;
	double DOFcap_frac, refine_frac, coarse_frac;

	// Testing
	unsigned int TestL2projection;

	// Structs
	struct S_ELEMENT *ELEMENT;
	struct S_VOLUME  *VOLUME, **Vgrp;
	struct S_FACE   *FACE;
};
extern struct S_DB DB;

#endif // DPG__S_DB_h__INCLUDED
