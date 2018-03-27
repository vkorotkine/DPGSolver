/*
 * The boundary condition definitions below **must** be the same as those in $(ROOT)/src/simulation/test_case/boundary/definitions_bc.h
 */

// Default characteristic length (Not used for transfinite meshes)
lc = 1.0;


EPS = 1e-15;


// Dummy value
GMSH_DUMMY = 314159265;


// Element types (Gmsh convention)
POINT = 15;
LINE  = 1;
TRI   = 2;
QUAD  = 3;
TET   = 4;
HEX   = 5;
WEDGE = 6;
PYR   = 7;

MIXED    = 20;
MIXED_TP = 21;
MIXED_HW = 22;


// PDE Names
ADVECTION     = 1;
DIFFUSION     = 2;
EULER         = 3;
NAVIER_STOKES = 4;


// PDE Specifiers
STEADY_DEFAULT                = 1;
STEADY_SUPERSONIC_VORTEX      = 2;
PERIODIC_PERIODIC_VORTEX      = 3;
STEADY_SUBSONIC_GAUSSIAN_BUMP = 4;
STEADY_PLANE_COUETTE          = 5;
STEADY_TAYLOR_COUETTE         = 6;
STEADY_JOUKOWSKI              = 7;
STEADY_VORTEX                 = 8;
STEADY_FREE_STREAM            = 9;
DEFAULT_STEADY                = 100; // Replace with STEADY_DEFAULT.


// GeomSpecifiers
EXTENSION_DISABLED = 0;
EXTENSION_ENABLED  = 1;

GEOM_NONE = 0; // Dummy value

GEOM_ADV_YL   = 1; // (ADV)ection (Y)-coord (L)eft
GEOM_ADV_XL   = 2; // (ADV)ection (X)-coord (L)eft
GEOM_ADV_XLR  = 5; // (ADV)ection (X)-coord (L)eft/(R)ight
GEOM_ADV_XYL  = 3; // (ADV)ection (X)(Y)-coords (L)eft
GEOM_ADV_XYZL = 4; // (ADV)ection (X)(Y)(Z)-coords (L)eft

GEOM_ADV_PERIODIC = 11; // (ADV)ection (PERIODIC) through all faces

GEOM_BC_ADIABATIC_O = 31; /// (B)oundary (C)ondition (ADIABATIC) (O)uter
GEOM_BC_DIABATIC_O  = 32; /// (B)oundary (C)ondition (DIABATIC) (O)uter

GEOM_2BEXP_0 = 0;
GEOM_2BEXP_1 = 1;
GEOM_2BEXP_2 = 2;

GEOM_CONFORMAL_HALF = 41;
GEOM_CONFORMAL_FULL = 42;


// MeshCurving Specifiers
STRAIGHT   = 0;
BLENDED    = 1;
PARAMETRIC = 2;



// Boundary conditions
BC_STEP_SC      = 10000;

BC_PERIODIC_MIN = 50;
PERIODIC_XL     = 51;
PERIODIC_XR     = 52;
PERIODIC_YL     = 53;
PERIODIC_YR     = 54;
PERIODIC_ZL     = 55;
PERIODIC_ZR     = 56;

PERIODIC_XL_REFLECTED_Y = 61;
PERIODIC_XR_REFLECTED_Y = 62;

GMSH_XLINE_MIN  = 1001;
GMSH_YLINE_MIN  = 2001;
GMSH_ZLINE_MIN  = 3001;
GMSH_XYFACE_MIN = 4001;
GMSH_XZFACE_MIN = 5001;
GMSH_YZFACE_MIN = 6001;
GMSH_XYZVOL_MIN = 7001;


// See comments in $DPG_ROOT/src/simulation/test_case/boundary/definitions_bc.h for when "ALT"ernate BC values are required.

BC_INFLOW       = 1; // Advection
BC_INFLOW_ALT1  = 2;
BC_INFLOW_ALT2  = 3;
BC_OUTFLOW      = 11;
BC_OUTFLOW_ALT1 = 12;
BC_OUTFLOW_ALT2 = 13;

BC_DIRICHLET      = 21; // Diffusion
BC_DIRICHLET_ALT1 = 22;
BC_NEUMANN        = 31;
BC_NEUMANN_ALT1   = 32;

BC_RIEMANN        = 101; // Euler
BC_SLIPWALL       = 102;
BC_BACKPRESSURE   = 103;
BC_TOTAL_TP       = 104;
BC_SUPERSONIC_IN  = 105;
BC_SUPERSONIC_OUT = 106;

BC_NOSLIP_ADIABATIC    = 111; // Navier-Stokes
BC_NOSLIP_DIABATIC     = 112;
BC_NOSLIP_ALL_ROTATING = 113;
