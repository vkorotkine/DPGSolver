// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__test_integration_NavierStokes_h__INCLUDED
#define DPG__test_integration_NavierStokes_h__INCLUDED
/// \file

/**	\brief Test various aspects of the Poisson solver implementation.
 *
 *	This currently includes integration tests for:
 *		- Equivalence between real and complex versions of functions;
 *		- Equivalence between running using different algorithms (with different flop counts);
 *		- Linearization;
 *		- Optimal convergence orders;
 *
 *	\note Both the collocated and non-collocated versions of the scheme are checked below as the operators are different
 *	      in each case.
 *
 */
void test_integration_NavierStokes
	(int nargc,  ///< Standard.
	 char **argv ///< Standard.
	);

#endif // DPG__test_integration_NavierStokes_h__INCLUDED
