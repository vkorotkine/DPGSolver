// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "test_integration_Poisson.h"

#include <stdlib.h>
#include <stdio.h>

#include "Parameters.h"

#include "test_code_integration_conv_order.h"
#include "test_code_integration_linearization.h"
#include "test_support.h"

#include "array_free.h"

/*
 *	Purpose:
 *		Test various aspects of the Poisson solver implementation:
 *			- Linearization;
 *			- Optimal convergence orders.
 *
 *	Comments:
 *
 *		It was very difficult to find a case where it was clear that blending of a curved boundary was leading to a loss
 *		of optimal convergence, despite the potentially unbounded mapping derivatives with h-refinement required for the
 *		optimal error estimate in Ciarlet(1972) (Theorem 5). As noted in Scott(1973) (p. 54), the mapping function and
 *		all of its derivatives are bounded IN TERMS OF THE BOUNDARY CURVATUVE AND ITS DERIVATIVES, which motivated the
 *		implementation of cases with geometry possessing high element curvatuve on coarse meshes. For these cases,
 *		optimal convergence is lost until the mesh has been refined "enough" (such that the element curvature is small)
 *		at which point it is recovered. Perhaps even more significant, the L2 error of the projection of the exact
 *		solution sometimes increased with mesh refinement on coarse meshes (run with Compute_L2proj = 1), a trend which
 *		was not observed for the computed solution. There was no modification found which could serve to fix this issue
 *		(such as the use of alternate blending functions or generalizations of the high-order blending used to fix the
 *		NIELSON blending as proposed by Lenoir(1986)). The mathematical and numerical support for this observation can
 *		be found in Zwanenburg(2017).
 *
 *		As a result of the conclusions of the study performed in Zwanenburg(2017) (and based on my current
 *		understanding), it seems that any of the optimal blending functions should give analogous results (SCOTT,
 *		SZABO_BABUSKA, LENOIR) and that the NORMAL surface parametrization is best (certainly in the 2D case). For the
 *		extension to 3D blending, the approach to be taken would be to extend the SB blending to PYR elements (using a
 *		combination of GH and SB for the LINE and TRI parts of the WEDGE element). The surface parametrization for the
 *		non-edge geometry nodes could use the NORMAL parametrization. An investigation into the optimal EDGE
 *		parametrization is still required.
 *
 *		It may also be noted that, despite converging at the same rate, the L2 error of the computed solution is
 *		significantly higher than that of the L2 projected exact solution for certain polynomial orders. Comparison with
 *		results obtained based on the DPG solver may be interesting here. (ToBeModified)
 *
 *
 *	*** IMPORTANT ***   Convergence Order Testing   *** IMPORTANT ***
 *
 *		It was found that optimal convergence was not possible to obtain using a series of uniformly refined TET meshes
 *		based on the "refine by splitting" algorithm in gmsh. However, optimal orders were recovered when a series of
 *		unstructured meshes consisting of TETs of decreasing volume was used.
 *
 *		Gmsh(2.14)'s "refine by splitting" algorithm splits each TET into 8 TETs, in a manner identical to the TET8
 *		h-refinement algorithm implemented in the code. Assuming an initially regular TET, with all edges having length
 *		= 2.0, is refined in this manner, the result will be 4 regular TETs with edge length = 1.0 and 4 other TETs with
 *		5 edges having length = 1.0 and 1 edge having length = sqrt(2.0). Taking a measure of the regularity of the mesh
 *		to be the ratio of the spheres enclosing and enclosed by the TETs, the regularity bound is violated through this
 *		refinement process if the splitting direction of the internal octohedron is taken at random, but not if taken
 *		consistently along the same axis (See Lenoir(1986) for the regularity requirement). Hence an alternative method
 *		of generating the refined mesh sequence must be employed (as compared to gmsh) to achieve the optimal
 *		convergence:
 *
 *			1) A sequence of unstructured (non-nested) meshes (gmsh) with diminishing volume gave optimal orders.
 *			2) TET8 refinement along a consistent axis (this code) gave optimal orders.
 *			3) TET6 and TET12 refinement algorithms both gave sub-optimal orders (to date).
 *				- Continued investigation may be made after a mesh quality improvement algorithm is implemented (e.g.
 *				  Peiro(2013)-Defining_Quality_Measures_for_Validation_and_Generation_of_High-Order_Tetrahedral_Meshes).
 *				  ToBeModified.
 *
 *	*** IMPORTANT ***   Convergence Order Testing   *** IMPORTANT ***
 *
 *	Notation:
 *
 *	References:
 *		Ciarlet(1972)-Interpolation_Theory_Over_Curved_Elements,_with_Applications_to_Finite_Element_Methods
 *		Scott(1973)-Finite_Element_Techniques_for_Curved_Boundaries
 *		Lenoir(1986)-Optimal_Isoparametric_Finite_Elements_and_Error_Estimates_for_Domains_Involving_Curved_Boundaries
 *		Zwanenburg(2017)-A_Necessary_High-Order_Meshing_Constraint_when_using_Polynomial_Blended_Geometry_Elements_for_
 *		                 Curved_Boundary_Representation
 */

void test_integration_Poisson(int nargc, char **argv)
{
	bool const RunTests_linearization = 1,
	           RunTests_conv_order    = 1;

	char **argvNew, *PrintName;

	argvNew    = malloc(2          * sizeof *argvNew);  // free
	argvNew[0] = malloc(STRLEN_MAX * sizeof **argvNew); // free
	argvNew[1] = malloc(STRLEN_MAX * sizeof **argvNew); // free
	PrintName  = malloc(STRLEN_MAX * sizeof *PrintName); // free

	// silence
	strcpy(argvNew[0],argv[0]);

	// **************************************************************************************************** //
	// Linearization Testing
	// **************************************************************************************************** //
	if (RunTests_linearization) {
		struct S_linearization *data_l = calloc(1 , sizeof *data_l); // free

		data_l->nargc     = nargc;
		data_l->argvNew   = argvNew;
		data_l->PrintName = PrintName;

		// 2D (Mixed TRI/QUAD mesh)
		test_linearization(data_l,"Poisson_MIXED2D");

		// 3D (TET mesh)
//		test_linearization(data_l,"Poisson_TET");
		test_print_warning("Poisson 3D testing needs to be updated");
		// Revisit when parametrized 3D geometry is available. (ToBeDeleted)

		free(data_l);
	} else {
		test_print_warning("Poisson linearization testing currently disabled");
	}

	// **************************************************************************************************** //
	// Convergence Order Testing
	// **************************************************************************************************** //
	if (RunTests_conv_order) {
		struct S_convorder *data_c = calloc(1 , sizeof *data_c); // free

		data_c->nargc     = nargc;
		data_c->argvNew   = argvNew;
		data_c->PrintName = PrintName;

//		test_conv_order(data_c,"Poisson_n-Ellipsoid_HollowSection_TRI");
//		test_conv_order(data_c,"Poisson_n-Ellipsoid_HollowSection_QUAD");
		test_conv_order(data_c,"Poisson_n-Ellipsoid_HollowSection_MIXED2D");

		free(data_c);
	} else {
		test_print_warning("Poisson convergence order testing currently disabled");
	}

	array_free2_c(2,argvNew);
	free(PrintName);
}
