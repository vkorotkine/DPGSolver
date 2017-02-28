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
	char         *TestCase, *MeshType, *MeshPath, *Form, *NodeType, *BasisType, *MeshFile, *Geometry;
	unsigned int d, ML, Vectorized, EFE, Collocated, Adapt, PGlobal, PMax, LevelsMax, Testing, *BumpOrder;
	int          Restart;

	// Parameters
	char         **NodeTypeG,
	             ***NodeTypeS,   ***NodeTypeF,   ***NodeTypeFrs, ***NodeTypeFrc,
	             ***NodeTypeIfs, ***NodeTypeIfc, ***NodeTypeIvs, ***NodeTypeIvc;
	unsigned int NP, AC, ExactGeom, InviscidFluxType, ViscousFluxType, ExplicitSolverType, PR, PP, PGs,
	             *PGc, *PF, *VFPartUnity, Blending, Blending_HO, Parametrization,
	             ***SF_BE, **PCs, **PCc, **PJs, **PJc, **PFrs, **PFrc, **PIfs, **PIfc, **PIvs, **PIvc;

	// Mesh
	unsigned int NVe, NPVe, NfMax, NfveMax, NveMax, NfrefMax, NETotal, NV, NVglobal, NGF, NVC, NGFC,
	             *PVe, *NE, *EType, *ETags, *EToVe, *EToPrt, *VToV, *VToF, *VToGF, *VToBC, *GFToVe, *VC, *GFC, *VeInfo;
	double *VeXYZ;

	// Structures
	unsigned int NECgrp;

	// Initialization
	bool         Viscous;
	char         *SolverType;
	unsigned int Nvar, Neq, OutputInterval, DOF0, SourcePresent;
	double       Xc, Yc, Rc, MInf, rhoInf, pInf, TInf, cInf, VInf, uInf, vInf, wInf, Rg, Cscale, pBack, p_Total, T_Total,
	             PeriodL, PeriodFraction, FinalTime,
	             rIn, MIn, rhoIn, pIn, VIn, rOut, Q0, KMin, KMax, GBa, GBb, GBc, aIn, bIn, cIn, aOut, bOut, cOut, lHR,
	             NSc, NSt, NS0, NS1, NS2, NS3, NS4, JSa, JSl, JSxL;

	// Solving
	unsigned int dof;

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
