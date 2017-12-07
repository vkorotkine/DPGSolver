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

#ifndef DPG__test_complex_numerical_flux_h__INCLUDED
#define DPG__test_complex_numerical_flux_h__INCLUDED
/** \file
 *  \brief Provides `complex` versions of containers and functions defined in \ref numerical_flux.h.
 */

#include "test_complex_boundary.h"

#include "def_templates_type_dc.h"
#include "def_templates_multiarray.h"
#include "def_templates_boundary.h"
#include "def_templates_numerical_flux.h"
#include "numerical_flux_T.h"
#include "undef_templates_type.h"
#include "undef_templates_multiarray.h"
#include "undef_templates_boundary.h"
#include "undef_templates_numerical_flux.h"

#if 0
struct Numerical_Flux_Input_c;
struct mutable_Numerical_Flux_c;
struct Simulation;

#include <stdbool.h>
#include "numerical_flux.h"
#include "test_complex_boundary.h"

/** \brief `complex` version of \ref compute_Numerical_Flux_fptr.
 *  \return See brief.
 *
 *  \param num_flux_i See brief.
 *  \param num_flux   See brief.
 */
typedef void (*compute_Numerical_Flux_c_fptr)
	(const struct Numerical_Flux_Input_c* num_flux_i,
	 struct mutable_Numerical_Flux_c* num_flux
	);

/// \brief Derived `complex` version of \ref Numerical_Flux_Input.
struct Numerical_Flux_Input_c {
	struct Numerical_Flux_Input num_flux_i; ///< Base \ref Numerical_Flux_Input.

	const bool has_complex_J; ///< Flag for whether the `complex` Jacobian terms should be computed.

	struct Boundary_Value_Input_c bv_l; ///< See brief.
	struct Boundary_Value_c       bv_r; ///< See brief.

	compute_Numerical_Flux_c_fptr compute_Numerical_Flux;     ///< See brief.
	compute_Numerical_Flux_c_fptr compute_Numerical_Flux_1st; ///< See brief.
	compute_Numerical_Flux_c_fptr compute_Numerical_Flux_2nd; ///< See brief.
};

/// \brief `complex` version of \ref Numerical_Flux.
struct Numerical_Flux_c {
	const struct const_Multiarray_c* nnf; ///< See brief.

	/// \brief See brief.
	struct Neigh_Info_NF_c {
		const struct const_Multiarray_c* dnnf_ds; ///< See brief.
		const struct const_Multiarray_c* dnnf_dg; ///< See brief.
	} neigh_info[2]; ///< See brief.
};

/// \brief `mutable` version of \ref Numerical_Flux_c.
struct mutable_Numerical_Flux_c {
	struct Multiarray_c* nnf; ///< See brief.

	/// See brief.
	struct m_Neigh_Info_NF_c {
		struct Multiarray_c* dnnf_ds; ///< See brief.
		struct Multiarray_c* dnnf_dg; ///< See brief.
	} neigh_info[2]; ///< See brief.
};

// Interface functions ********************************************************************************************** //

/** \brief Constructor for a \ref Numerical_Flux_Input_c container.
 *  \return See brief. */
struct Numerical_Flux_Input_c* constructor_Numerical_Flux_Input_c
	(const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Destructor for a \ref Numerical_Flux_Input_c container.
void destructor_Numerical_Flux_Input_c
	(struct Numerical_Flux_Input_c* num_flux_i ///< Standard.
	);

/** \brief Constructor for a \ref Numerical_Flux_c container.
 *  \return See brief. */
struct Numerical_Flux_c* constructor_Numerical_Flux_c
	(const struct Numerical_Flux_Input_c* num_flux_i ///< \ref Numerical_Flux_Input_c.
	);

/// \brief Destructor for a \ref Numerical_Flux_c container.
void destructor_Numerical_Flux_c
	(struct Numerical_Flux_c* num_flux ///< Standard.
	);

/// \brief `complex` version of \ref compute_Numerical_Flux_1.
void compute_Numerical_Flux_c_1
	(const struct Numerical_Flux_Input_c* num_flux_i, ///< See brief.
	 struct mutable_Numerical_Flux_c* num_flux        ///< See brief.
	);

/// \brief `complex` version of \ref compute_Numerical_Flux_12.
void compute_Numerical_Flux_c_12
	(const struct Numerical_Flux_Input_c* num_flux_i, ///< See brief.
	 struct mutable_Numerical_Flux_c* num_flux        ///< See brief.
	);
#endif

#endif // DPG__test_complex_numerical_flux_h__INCLUDED
