// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__test_integration_Poisson_h__INCLUDED
#define DPG__test_integration_Poisson_h__INCLUDED
/// \file

/**	\brief Test various aspects of the Poisson solver implementation.
 *
 *	This currently includes integration tests for:
 *		- Linearization;
 *		- Optimal convergence orders;
 *
 *	\section s1_poisson General Comments
 *
 *	\todo Likely move this section (and subsections) to another place in the code.
 *
 *	\subsection s1_1_poisson Optimal Convergence for Curved Elements
 *
 *	It was very difficult to find a case where it was clear that blending of a curved boundary was leading to a loss of
 *	optimal convergence, despite the potentially unbounded mapping derivatives with h-refinement required for the
 *	optimal error estimate of Ciarlet et al. (Theorem 5) \cite Ciarlet1972. As noted by Scott (p. 54)
 *	\cite Scott_thesis, the mapping function and all of its derivatives are bounded IN TERMS OF THE BOUNDARY CURVATURE
 *	AND ITS DERIVATIVES, which motivated the implementation of cases with geometry possessing high element curvatuve on
 *	coarse meshes. For these cases, optimal convergence is lost until the mesh has been refined "enough" (such that the
 *	element curvature is "small") at which point it is recovered. Perhaps even more significant, the L2 error of the
 *	projection of the exact solution sometimes increased with mesh refinement on coarse meshes (enable computation of
 *	the L2 projection), a trend which was not observed for the computed solution. The mathematical and numerical support
 *	for this observation can be found in Zwanenburg(2017).
 *
 *	As a result of the conclusions of the study performed in Zwanenburg(2017) (and based on my current understanding),
 *	it seems that any of the optimal blending functions should give analogous results (SCOTT, SZABO_BABUSKA, LENOIR) and
 *	that the NORMAL surface parametrization is best (certainly in the 2D case). For the extension to 3D blending, the
 *	approach to be taken would be to extend the SB blending to PYR elements (using a combination of GH and SB for the
 *	LINE and TRI parts of the WEDGE element). The surface parametrization for the non-edge geometry nodes could use the
 *	NORMAL parametrization. An investigation into the optimal EDGE parametrization is still required.
 *
 *	It may also be noted that, despite converging at the same rate, the L2 error of the computed solution is
 *	significantly higher than that of the L2 projected exact solution for certain polynomial orders. Comparison with
 *	results obtained based on the DPG solver may be interesting here. \todo Update this comment after performing the
 *	study.
 *
 *	\todo Add reference to Zwanenburg(2017) when submitted.
 *
 *	\subsection s1_2_poisson Convergence Order Testing
 *
 *	It was found that optimal convergence was not possible to obtain using a series of uniformly refined TET meshes
 *	based on the "refine by splitting" algorithm in gmsh. However, optimal orders were recovered when a series of
 *	unstructured meshes consisting of TETs of decreasing volume was used.
 *
 *	Gmsh(2.14)'s "refine by splitting" algorithm splits each TET into 8 TETs, in a manner identical to the TET8
 *	h-refinement algorithm implemented in the code. Assuming an initially regular TET, with all edges having length =
 *	2.0, is refined in this manner, the result will be 4 regular TETs with edge length = 1.0 and 4 other TETs with 5
 *	edges having length = 1.0 and 1 edge having length = sqrt(2.0). Taking a measure of the regularity of the mesh to be
 *	the ratio of the radii of the spheres enclosing and enclosed by the TETs, the regularity bound is violated through
 *	this refinement process if the splitting direction of the internal octohedron is taken at random, but not if taken
 *	consistently along the same axis (See Lenoir \cite Lenoir1986 for the regularity requirement). Hence an alternative
 *	method of generating the refined mesh sequence must be employed (as compared to gmsh) to achieve the optimal
 *	convergence:
 *
 *	1. A sequence of unstructured (non-nested) meshes (gmsh) with diminishing volume gave optimal orders.
 *	2. TET8 refinement along a consistent axis (method employed in this code) gave optimal orders.
 *	3. TET6 and TET12 refinement algorithms both gave sub-optimal orders.
 *		- Continued investigation may be made after a mesh quality improvement algorithm is implemented (e.g. see
 *		  Gargallo-Peiro et al. \cite Gargallo-Peiro2014). \todo Continue this investigation.
 */
void test_integration_Poisson
	(int nargc,  ///< Standard.
	 char **argv ///< Standard.
	);

#endif // DPG__test_integration_Poisson_h__INCLUDED
