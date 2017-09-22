// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__test_integration_Euler_h__INCLUDED
#define DPG__test_integration_Euler_h__INCLUDED
/// \file

/**	\brief Test various aspects of the Euler solver implementation.
 *
 *	This currently includes integration tests for:
 *		- Equivalence between real and complex versions of functions;
 *		- Equivalence between running using different algorithms (with different flop counts);
 *		- Linearization;
 *		- Optimal convergence orders;
 *
 *	\attention Optimal convergence orders for 3D curved meshes which are not associated with extruded 2D meshes has so
 *	           far not been obtained. This is potentially a result of mesh regularity issues or because memory
 *	           constraints do not allow for the attainment of the asymptotic regime.
 */
void test_integration_Euler
	(int nargc,  ///< Standard.
	 char **argv ///< Standard.
	);

#endif // DPG__test_integration_Euler_h__INCLUDED
