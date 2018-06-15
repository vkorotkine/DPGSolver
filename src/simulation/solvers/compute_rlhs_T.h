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
 *  \brief Provides templated functions used for computing the all contributions to the right and left-hand side
 *         (rlhs) terms of supported schemes.
 */

#include "def_templates_compute_rlhs.h"
#include "def_templates_multiarray.h"
#include "def_templates_flux.h"

struct Flux_T;

/// \brief Container for the reference flux related parameters.
struct Flux_Ref_T {
	const struct const_Multiarray_T* fr;       ///< The reference flux.
	const struct const_Multiarray_T* dfr_ds;   ///< The reference flux Jacobians with respect to the solution.
	const struct const_Multiarray_T* dfr_dg;   ///< The reference flux Jacobians with respect to the solution gradients.
	const struct const_Multiarray_T* d2fr_ds2; ///< The reference flux Hessians  with respect to the solution.
};

/** \brief Constructor for a \ref Flux_Ref_T container.
 *  \return See brief.
 *
 *  This function constructs the reference fluxes by multiplying the fluxes with the appropriate metric terms as
 *  specified by Zwanenburg et al. (eq. (B.3), \cite Zwanenburg2016).
 *
 *  The memory layout of the reference flux and reference flux Jacobian terms is that of the corresponding physical flux
 *  terms with the second extent (dim) moved to the last extent. For example:
 *  - df_ds: (nodes,dim,eq,var) -> dfr_ds: (nodes,eq,var,dim).
 *  This memory layout was chosen such that terms are grouped by dimension, allowing for differentiation operators to be
 *  applied efficiently; note that **this is not the same ordering as that used for the physical flux**. Please consult
 *  \ref compute_geometry_volume_T for the ordering of the metric terms if desired.
 */
struct Flux_Ref_T* constructor_Flux_Ref_T
	(const struct const_Multiarray_T*const m, ///< The metric terms.
	 const struct Flux_T*const flux           ///< The physical \ref Flux_T.
		);

/// \brief Destructor for a \ref Flux_Ref_T container.
void destructor_Flux_Ref_T
	(struct Flux_Ref_T*const flux_ref ///< Standard.
		);

#include "undef_templates_compute_rlhs.h"
#include "undef_templates_multiarray.h"
#include "undef_templates_flux.h"
