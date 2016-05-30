// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__parameters_h__INCLUDED
#define DPG__parameters_h__INCLUDED

/*
 *	Purpose:
 *		Define parameters.
 *
 *	Comments:
 *		See gmsh_reader.c for further discussion regarding periodic related mesh file parameters.
 *
 *	Notation:
 *		DMAX            : (MAX)imum (D)imension
 *		BC_STEP_SC      : (B)oundary(C)ondition step between (S)traight and (C)urved BCs.
 *		BC_PERIODIC_MIN : (B)oundary(C)ondition (PERIODIC) (MIN)imum.
 *		NVEMAX          : (MAX)imum (N)umber of (VE)rtices for an element.
 *		NFVEMAX         : (MAX)imum (N)umber of (F)ACET (VE)rtices for an element.
 *		NFREFMAX        : (MAX)imum (N)umber of (F)ACET (REF)inements.
 *		NFMAX           : (MAX)imum (N)umber of (F)ACET for an element.
 *		NFMIXEDMAX      : (MAX)imum (N)umber of (MIXED) (F)ACETs for an element.
 *		NESUBCMAX       : (MAX)imum (N)umber of (E)lement (SUB)(C)lasses
 *		NFORDMAX        : (MAX)imum (N)umber of (F)ACET (ORD)ering possibilities
 *
 *	References:
 *
 */

// Magic numbers
#define DMAX            3

#define BC_STEP_SC      10000
#define BC_PERIODIC_MIN 50
#define PERIODIC_XL     51
#define PERIODIC_XR     52
#define PERIODIC_YL     53
#define PERIODIC_YR     54
#define PERIODIC_ZL     55
#define PERIODIC_ZR     56

#define GMSH_XLINE_MIN  1001
#define GMSH_YLINE_MIN  2001
#define GMSH_ZLINE_MIN  3001
#define GMSH_XYFACE_MIN 4001
#define GMSH_XZFACE_MIN 5001
#define GMSH_YZFACE_MIN 6001
#define GMSH_XYZVOL_MIN 7001

#define NVEMAX          8 // HEX
#define NFVEMAX         4 // QUAD
#define NFREFMAX        9 // QUAD
#define NFMAX           6 // HEX
#define NFMIXEDMAX      2 // WEDGE/PYR (TRI + QUAD)
#define NESUBCMAX       2 // WEDGE (TRI + LINE)
#define NFORDMAX        8 // QUAD

#define NFREFMAXPOINT   1
#define NFREFMAXLINE    3
#define NFREFMAXTRI     5
#define NFREFMAXQUAD    9

// Common variables
#define PI    3.1415926535897932
#define GAMMA 1.4
#define GM1   0.4

// Element types (Gmsh convention)
#define POINT 15
#define LINE  1
#define TRI   2
#define QUAD  3
#define TET   4
#define HEX   5
#define WEDGE 6
#define PYR   7

// Element classes
#define C_TP    0
#define C_SI    1
#define C_PYR   2
#define C_WEDGE 3

// Node Tolerance (for physical coordinate comparison)
#define NODETOL      1.0e-10
#define NODETOL_MESH 1.0e-5

// Value close to double machine zero
#define EPS     1.0e-15

// Min/Max String Lengths
// Min/Max length slightly less than 2^6, 2^9
#define STRLEN_MIN 60
#define STRLEN_MAX 508

// Macros
#define max(a,b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b); _a > _b ? _a : _b; })
#define min(a,b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b); _a < _b ? _a : _b; })

#define sign(a) ({ __typeof__ (a) _a = (a); (_a > 0) ? 1 : ((_a < 0) ? -1 : 0); })



#endif // DPG__parameters_h__INCLUDED
