// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__definitions_bc_h__INCLUDED
#define DPG__definitions_bc_h__INCLUDED
/**	\file
 *	\brief Provides the definitions related to the boundary conditions.
 *
 *	The naming convention for the periodic boundary conditions omits the 'BC_' prefix as they are only used for
 *	preprocessing and not as actual boundary conditions.
 *
 *	\note The value of these constants must be identical to those in the $ROOT/input/meshes/parameters.geo file.
 */

///\{ \name The step in the boundary conditions between straight and curved BCs.
#define BC_STEP_SC 10000
///\}

///\{ \name The values for the various boundary conditions after taking the modulus with \ref BC_STEP_SC.

// General
#define PERIODIC_XL 51
#define PERIODIC_XR 52
#define PERIODIC_YL 53
#define PERIODIC_YR 54
#define PERIODIC_ZL 55
#define PERIODIC_ZR 56

// Advection
#define BC_INFLOW  13
#define BC_OUTFLOW 14

// Poisson
#define BC_DIRICHLET 11
#define BC_NEUMANN   12

// Euler
#define BC_RIEMANN        1
#define BC_SLIPWALL       2
#define BC_BACKPRESSURE   3
#define BC_TOTAL_TP       4
#define BC_SUPERSONIC_IN  5
#define BC_SUPERSONIC_OUT 6

// Navier-Stokes
#define BC_NOSLIP_T         7
#define BC_NOSLIP_ADIABATIC 8
///\}


#endif // DPG__definitions_bc_h__INCLUDED
