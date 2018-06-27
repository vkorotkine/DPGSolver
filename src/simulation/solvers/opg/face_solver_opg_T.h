/* {{{
This file is part of DPGSolver.

DPGSolver is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or any later version.

DPGSolver is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
Public License for more details.

You should have received a copy of the GNU General Public License along with DPGSolver.  If not, see
<http://www.gnu.org/licenses/>.
}}} */
/** \file
 *  \brief Provides the interface for the templated \ref OPG_Solver_Face_T container and associated functions.
 *
 *  These faces are needed by the 'O'ptimal 'P'etrov 'G'alerkin solver functions.
 */

#include "def_templates_face_solver_opg.h"
#include "def_templates_matrix.h"
#include "def_templates_penalty_opg.h"
#include "def_templates_compute_rlhs.h"
#include "def_templates_numerical_flux.h"

struct Flux_Ref_T;
struct Numerical_Flux_T;
struct Solver_Face_T;
struct Solver_Storage_Implicit;

/** \brief Pointer to functions adding the penalty term to required boundary faces.
 *
 *  \param flux     \ref Flux_T.
 *  \param num_flux \ref Numerical_Flux_T.
 *  \param s_face   \ref Solver_Face_T.
 *  \param ssi      \ref Solver_Storage_Implicit.
 *
 *  As the (bilinear) form obtained only with the dg-like terms results in the specification of the test function along
 *  each charateristic only to within an arbitrary constant, an additional penalty term is added to ensure that a
 *  sufficient number of boundary conditions are set for the test function such that the system is solvable. The penalty
 *  term (added to the rhs) takes the form:
 *
 *  \f$ \eps^{-1} <v,g-w>_{\Gamma^\text{characteristic out}} \forall v \f$
 *
 *  where:
 *  - \f$ \eps \f$ is chosen as \f$ \eps = C h^{p_\text{test}+1} > 0 \f$ following the \cite Barrett1986;
 *  - \f$ v \f$ represents the test function which is being used to test the equation;
 *  - \f$ w \f$ represents the test function which is being solved for;
 *  - \f$ g \f$ represents the (arbitrary) boundary value;
 *  - \f$ \Gamma^\text{characteristic out}} \f$ represents boundaries having outgoing characteristics.
 *
 *  Note that this penalty __must not__ be applied on every face with outgoing characteristics; doing so results in
 *  overconstraining the problem. The simplest example of this was discussed by Brunken et al., section 3.2
 *  \cite Brunken2018 where it was noted that applying Dirichlet BCs for the test function along all outflow boundaries
 *  for linear advection restricts the value of the solution at the outflow corner to zero; consequently, the
 *  non-trivial constant solution cannot even be recovered... They advocate extending the domain to remedy the problem
 *  but this is addressed here by eliminating the imposed boundary conditions which result in non-physical solution
 *  restrictions. The correct method for imposing the boundary conditions is still under investigation.
 *
 *  In the case of the linear advection equation, the term can be simplified by choosing \f$ g = 0\f$ and \f$ w_0 = 0\f$
 *  (i.e. choosing the initial test function values as being zero on the outflow boundary). It thus has no contribution
 *  to the rhs but is present in the lhs (linearization) in this case.
 */
typedef void (*constructor_rlhs_f_b_test_penalty_T)
	(const struct Flux_Ref_T*const flux_r,
	 const struct Numerical_Flux_T*const num_flux,
	 struct Solver_Face_T*const s_face,
	 struct Solver_Storage_Implicit*const ssi
	 );

/// \brief Container for data relating to the OPG solver faces.
struct OPG_Solver_Face_T {
	struct Solver_Face_T face; ///< The base \ref Solver_Face_T.

	const struct const_Matrix_T* m_inv; ///< The inverse mass matrix.

	/// Version of \ref constructor_rlhs_f_b_test_penalty_T for rhs (and lhs if applicable) term(s).
	constructor_rlhs_f_b_test_penalty_T constructor_rlhs_penalty;

	/// The type of 'b'oundary 'c'ondition to set for \ref Solver_Volume_T::test_s_coef on the current face.
	int bc_test_s;
};

/// \brief Constructor for a derived \ref OPG_Solver_Face_T.
void constructor_derived_OPG_Solver_Face_T
	(struct Face* face_ptr,       ///< Pointer to the face.
	 const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Destructor for a derived \ref OPG_Solver_Face_T.
void destructor_derived_OPG_Solver_Face_T
	(struct Face* face_ptr ///< Pointer to the face.
	);

/** \brief Get the pointer to the appropriate \ref OPG_Solver_Element::cv0_vt_fc operator.
 *  \return See brief. */
const struct Operator* get_operator__cv0_vt_fc_T
	(const int side_index,                           ///< The index of the side of the face under consideration.
	 const struct OPG_Solver_Face_T*const opg_s_face ///< Standard.
		);

/** \brief Get the pointer to the appropriate \ref OPG_Solver_Element::cv1_vt_fc operator.
 *  \return See brief. */
struct Multiarray_Operator get_operator__cv1_vt_fc_T
	(const int side_index,                           ///< The index of the side of the face under consideration.
	 const struct OPG_Solver_Face_T*const opg_s_face ///< Standard.
	);

#include "undef_templates_face_solver_opg.h"
#include "undef_templates_matrix.h"
#include "undef_templates_penalty_opg.h"
#include "undef_templates_compute_rlhs.h"
#include "undef_templates_numerical_flux.h"
