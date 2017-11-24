/*
 * The boundary condition definitions below **must** be the same as those in $(ROOT)/src/constants/constants_bc.h
 */

// Default characteristic length (Not used for transfinite meshes)
lc = 1.0;


EPS = 1e-15;


// Dummy value
GMSH_DUMMY = -999;


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
ADVECTION    = 1;
POISSON      = 2;
EULER        = 3;
NAVIERSTOKES = 4;


// PDE Specifiers
STEADY_DEFAULT                = 1;
STEADY_SUPERSONIC_VORTEX      = 2;
PERIODIC_PERIODIC_VORTEX      = 3;
STEADY_SUBSONIC_GAUSSIAN_BUMP = 4;
STEADY_PLANE_COUETTE          = 5;
STEADY_TAYLOR_COUETTE         = 6;


// GeomSpecifiers
EXTENSION_DISABLED = 0;
EXTENSION_ENABLED  = 1;

GEOM_NONE = 0; // Dummy value

GEOM_AR_1 = 1;
GEOM_AR_2 = 2;
GEOM_AR_3 = 3;

GEOM_ADV_YL   = 1; // (ADV)ection (Y)-coord (L)eft
GEOM_ADV_XL   = 2; // (ADV)ection (X)-coord (L)eft
GEOM_ADV_XYL  = 3; // (ADV)ection (X)(Y)-coords (L)eft
GEOM_ADV_XYZL = 4; // (ADV)ection (X)(Y)(Z)-coords (L)eft

GEOM_2BEXP_0 = 0;
GEOM_2BEXP_1 = 1;
GEOM_2BEXP_2 = 2;


// MeshCurving Specifiers
STRAIGHT   = 0;
CURVED     = 1;
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
