// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

/*
 *	Comments:
 *		The definitions below should be the same as those in $(ROOT)/include/Parameters.h
 */

// Default characteristic length (Not used for transfinite meshes)
lc = 1.0;


// Element types (Gmsh convention)
POINT = 15;
LINE  = 1;
TRI   = 2;
QUAD  = 3;
TET   = 4;
HEX   = 5;
WEDGE = 6;
PYR   = 7;

MIXED2D    = 20;
MIXED3D    = 21;
MIXED3D_TP = 22;
MIXED3D_HW = 23;


// PDE Names
ADVECTION    = 0;
POISSON      = 1;
EULER        = 2;
NAVIERSTOKES = 3;


// PDE Specifiers
PERIODIC = 0;


// GeomSpecifiers
EXTENSION_DISABLED = 0;
EXTENSION_ENABLED  = 1;

GEOM_AR_1 = 1;
GEOM_AR_2 = 2;
GEOM_AR_3 = 3;

GEOM_ADV_NONE = 0; // Dummy
GEOM_ADV_YL   = 1; // (ADV)ection (Y)-coord (L)eft


// MeshCurving Specifiers
STRAIGHT   = 0;
CURVED     = 1;
TOBECURVED = 2;



// Boundary conditions
BC_STEP_SC      = 10000;
BC_PERIODIC_MIN = 50;
PERIODIC_XL     = 51;
PERIODIC_XR     = 52;
PERIODIC_YL     = 53;
PERIODIC_YR     = 54;
PERIODIC_ZL     = 55;
PERIODIC_ZR     = 56;

GMSH_XLINE_MIN  = 1001;
GMSH_YLINE_MIN  = 2001;
GMSH_ZLINE_MIN  = 3001;
GMSH_XYFACE_MIN = 4001;
GMSH_XZFACE_MIN = 5001;
GMSH_YZFACE_MIN = 6001;
GMSH_XYZVOL_MIN = 7001;


BC_RIEMANN        = 1; // Euler
BC_SLIPWALL       = 2;
BC_BACKPRESSURE   = 3;
BC_TOTAL_TP       = 4;
BC_SUPERSONIC_IN  = 5;
BC_SUPERSONIC_OUT = 6;

BC_NOSLIP_T         = 7; // Navier-Stokes
BC_NOSLIP_ADIABATIC = 8;

BC_DIRICHLET    = 11; // Poisson
BC_NEUMANN      = 12;

BC_INFLOW       = 13; // Advection
BC_OUTFLOW      = 14;
