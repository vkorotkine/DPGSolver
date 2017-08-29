// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__test_integration_Advection_h__INCLUDED
#define DPG__test_integration_Advection_h__INCLUDED
/// \file

/**	\brief Test various aspects of the Advection solver implementation.
 *
 *	This currently includes integration tests for:
 *		- Linearization;
 *		- Optimal convergence orders;
 *
 *	\attention Convergence order tests failed (TRI/QUAD) meshes when using a GLL nodal basis (but succeeded for WSH/GL
 *	           nodal and all modal). The failure was due to KSPConvergedReason = -11. Aditya tried with a
 *	           pseudotimestepping code and obtained optimal orders for P2 EQ nodes. \todo INVESTIGATE
 */
void test_integration_Advection
	(int nargc,  ///< Standard.
	 char **argv ///< Standard.
	);

#endif // DPG__test_integration_Advection_h__INCLUDED
