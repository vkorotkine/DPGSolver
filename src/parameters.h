// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__parameters_h__INCLUDED
#define DPG__parameters_h__INCLUDED

/*
 *	Purpose:
 *		Define parameters.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 *
 */

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
