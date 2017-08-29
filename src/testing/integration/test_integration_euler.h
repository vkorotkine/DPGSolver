// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__test_integration_euler_h__INCLUDED
#define DPG__test_integration_euler_h__INCLUDED
/** \file
 *	Provides functionality for Euler integration testing.
 */

/** \brief Tests various components of the Euler solver.
 *
 *	Integration tests are currently provided for:
 *		- Equivalence between real and complex versions of functions;
 *		- Equivalence between running using different algorithms;
 *		- Linearizations;
 *		- Expected convergence orders;
 *
 *	\attention Optimal convergence orders for 3D curved meshes which are not associated with extruded 2D meshes have so
 *	           far not been obtained. This is potentially a result of mesh regularity issues or because memory
 *	           constraints do not allow for the attainment of the asymptotic regime.
 */
void test_integration_euler
	(struct Test_Info*const test_info /// \ref Test_Info.
	);

#endif // DPG__test_integration_euler_h__INCLUDED
