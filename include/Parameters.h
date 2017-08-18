// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__Parameters_h__INCLUDED
#define DPG__Parameters_h__INCLUDED

/*
 *	Purpose:
 *		Define parameters.
 *
 *	Comments:
 *		See gmsh_reader.c for further discussion regarding periodic related mesh file parameters.
 *
 *	Notation:
 *		NEC             : (N)umber of (E)lement (C)lasses: TP, SI, PYR
 *		BC_PERIODIC_MIN : (B)oundary(C)ondition (PERIODIC) (MIN)imum.
 *		NVEMAX          : (MAX)imum (N)umber of (VE)rtices for an element.
 *		NEVEMAX         : (MAX)imum (N)umber of (E)DGE (VE)rtices for an element.
 *		NEREFMAX        : (MAX)imum (N)umber of (E)DGE (REF)inements.
 *		NFREFMAX        : (MAX)imum (N)umber of (F)ACE (REF)inements.
 *		NVREFMAX        : (MAX)imum (N)umber of (V)OLUME (REF)inements.
 *		NVREFSFMAX      : (MAX)imum (N)umber of (V)OLUME (REF)inements if using (S)um (F)actorized operators.
 *		NEMAX           : (MAX)imum (N)umber of (E)DGEs for an element.
 *		NFMIXEDMAX      : (MAX)imum (N)umber of (MIXED) (F)ACEs for an element.
 *		NESUBCMAX       : (MAX)imum (N)umber of (E)lement (SUB)(C)lasses
 *		NFORDMAX        : (MAX)imum (N)umber of (F)ACE (ORD)ering possibilities
 *		NREFVVARMAX     : (MAX)imum (N)umber of h-adaptive (REF)ined (V)olume (VAR)iations
 *		NSUBFMAX        : (MAX)imum (N)umber of h-adaptive (SUB)-(F)aces (on each FACE).
 *		NVISUBFMAX      : (MAX)imum (N)umber of h-adaptive (V)OLUME (I)nternal (SUB)-(F)aces (within the VOLUME).
 *		NSIBMAX         : (MAX)imum (N)umber of (SIB)lings on the same level after h-refinement.
 *		NEHREFMAX       : (MAX)imum (N)umber of (E)LEMENT types present in (H)-(REF)ined ELEMENT
 *		NVEINFO         : (N)umber of pieces of (INFO)rmation associated with each (VE)rtex.
 *
 *	References:
 *
 */

#include "constants_core.h"
#include "constants_elements.h"
#include "constants_bc.h"

#ifndef TEST
	#define TEST 0
#endif // TEST

// Alternate names
#define CBRM CblasRowMajor
#define CBCM CblasColMajor
#define CBT  CblasTrans
#define CBNT CblasNoTrans

// Magic numbers
#define NEC             3

// Geometry related parameters
#define GORDON_HALL       1
#define SZABO_BABUSKA     2
#define SCOTT             3
#define NIELSON           4

#define ARC_LENGTH        1
#define RADIAL_PROJECTION 2
#define NORMAL            3
#define ORDER_H           4

// ELEMENT related numbers
#define NVEMAX          8  // HEX
#define NEVEMAX         2  // LINE
#define NEREFMAX        3  // LINE
#define NFREFMAX        9  // QUAD
#define NVREFMAX        27 // HEX
#define NVREFSFMAX      5  // TRI
#define NEMAX           12 // HEX
#define NFMIXEDMAX      2  // WEDGE/PYR (TRI + QUAD)
#define NESUBCMAX       2  // WEDGE (TRI + LINE)
#define NFORDMAX        8  // QUAD
#define NREFVVARMAX     7  // HEX
#define NSUBFMAX        4  // QUAD/TRI (Isotropic refinement)
#define NVISUBFMAX      16 // TET (Isotropic refinement to 12 TETs)
#define NSIBMAX         12 // TET (12 TET)
#define NVEINFO         4

// Cubature related numbers
#define PIvcMaxTET 10
#define PIvcMaxPYR 10

// PDE Names
#define PDE_ADVECTION    1
#define PDE_POISSON      2
#define PDE_EULER        3
#define PDE_NAVIERSTOKES 4


// Solver related parameters
#define RK3_SSP 0
#define RK4_LS  1
#define EULER   2

// Allowed adaptation options
#define ADAPT_0  0
#define ADAPT_P  1
#define ADAPT_H  2
#define ADAPT_HP 3

// Adaptation flags
#define HREFINE  1
#define PREFINE  2
#define HCOARSE  3
#define PCOARSE  4
#define HPREFINE 5
#define HPCOARSE 6
#define HDELETE  10

// h-refinement related numbers
#define NEHREFMAX    2 // PYR (PYR + TET)

#define TET8         0
#define TET12        1
#define TET6         2

#define NREFMAXPOINT 1
#define NREFMAXLINE  3
#define NREFMAXTRI   5
#define NREFMAXQUAD  9
#define NREFMAXTET   13
#define NREFMAXHEX   27
#define NREFMAXWEDGE 15
#define NREFMAXPYR   11

// Common variables
#define PI    3.1415926535897932
#define GAMMA 1.4
#define GM1   0.4
#define GM3   -1.6

// Element classes
#define C_TP    0
#define C_SI    1
#define C_PYR   2
#define C_WEDGE 3

// Flux types
#define FLUX_LF     0 // Change to 1-3 (ToBeModified)
#define FLUX_ROE    1
#define FLUX_UPWIND 2

#define FLUX_BR2  11
#define FLUX_CDG2 12

// Solver Methods
#define METHOD_DG   1 /// \todo ToBeDeleted (Only need this when simulation is included -> defined in Simulation.h)
#define METHOD_HDG  2

// Tolerances
#define NODETOL      1.0e-10
#define NODETOL_MESH 1.0e-5

#define EPS        1.0e-15
#define EPS_PHYS   1.0e-13
#define SQRT_EPS   3.162277660168379e-08
#define CX_STEP    1.0e-30
#define REFINE_TOL 1.0e-10 // Decrease for additional refinement
#define COARSE_TOL 1.0e-4 // Increase for additional coarsening

// Min/Max String Lengths
// Min/Max length slightly less than 2^6, 2^9
#define STRLEN_MIN 60
#define STRLEN_MAX 508







///\{ \deprecated
#define BC_PERIODIC_MIN 50

#define GMSH_XLINE_MIN  1001
#define GMSH_YLINE_MIN  2001
#define GMSH_ZLINE_MIN  3001
#define GMSH_XYFACE_MIN 4001
#define GMSH_XZFACE_MIN 5001
#define GMSH_YZFACE_MIN 6001
#define GMSH_XYZVOL_MIN 7001
///\}

#endif // DPG__Parameters_h__INCLUDED
