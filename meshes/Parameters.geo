// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

/*
 *	Comments:
 *		The definitions below should be the same as those in $(ROOT)/include/Parameters.h
 */

// Default characteristic length (Not used for transfinite meshes)
lc = 1.0;


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


BC_RIEMANN        = 1;
BC_SLIPWALL       = 2;
BC_BACKPRESSURE   = 3;
BC_TOTAL_TP       = 4;
BC_SUPERSONIC_IN  = 5;
BC_SUPERSONIC_OUT = 6;

BC_DIRICHLET    = 11;
BC_NEUMANN      = 12;


// Element types (Gmsh convention)
POINT = 15;
LINE  = 1;
TRI   = 2;
QUAD  = 3;
TET   = 4;
HEX   = 5;
WEDGE = 6;
PYR   = 7;
