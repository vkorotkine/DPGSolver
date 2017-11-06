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

#ifndef DPG__geometry_element_h__INCLUDED
#define DPG__geometry_element_h__INCLUDED
/** \file
 *  \brief Provides the interface for the derived \ref Geometry_Element container and associated functions.
 *
 *  \note `const` and non-`const` versions of \ref Geometry_Element must have identical members and layout.
 */

#include "element.h"

struct Simulation;

/// \brief Container for data relating to the geometry elements.
struct Geometry_Element {
	struct const_Element element; ///< Base \ref const_Element.

/// \todo Replace the s/c with [2].
	const struct Multiarray_Operator* vc0_vgc_vgc; ///< See notation in \ref element_operators.h.

	const struct Multiarray_Operator* cv1_vgs_vcs; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* cv1_vgc_vcc; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* cv1_vgs_vms; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* cv1_vgc_vmc; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* vv0_vms_vcs; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* vv0_vmc_vcc; ///< See notation in \ref element_operators.h.

	const struct Multiarray_Operator* cv0_vgs_fcs; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* cv0_vgs_fcc; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* cv0_vgc_fcs; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* cv0_vgc_fcc; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* vv0_vms_fcs; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* vv0_vms_fcc; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* vv0_vmc_fcs; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* vv0_vmc_fcc; ///< See notation in \ref element_operators.h.

	// Tensor-product sub-operators.
	const struct Multiarray_Operator* cv0_vgs_vcs; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* cv0_vgc_vcc; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* cv0_vgs_vms; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* cv0_vgc_vmc; ///< See notation in \ref element_operators.h.

	const struct Multiarray_Operator* cv0_vgs_vcc; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* cv0_vgc_vcs; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* vv0_vms_vcc; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* vv0_vmc_vcs; ///< See notation in \ref element_operators.h.
};

/// \brief `const` version of the \ref Geometry_Element container.
struct const_Geometry_Element {
	const struct const_Element element; ///< Defined in \ref Geometry_Element.

	const struct Multiarray_Operator*const vc0_vgc_vgc; ///< Defined in \ref Geometry_Element.

	const struct Multiarray_Operator*const cv1_vgs_vcs; ///< Defined in \ref Geometry_Element.
	const struct Multiarray_Operator*const cv1_vgc_vcc; ///< Defined in \ref Geometry_Element.
	const struct Multiarray_Operator*const cv1_vgs_vms; ///< Defined in \ref Geometry_Element.
	const struct Multiarray_Operator*const cv1_vgc_vmc; ///< Defined in \ref Geometry_Element.
	const struct Multiarray_Operator*const vv0_vms_vcs; ///< Defined in \ref Geometry_Element.
	const struct Multiarray_Operator*const vv0_vmc_vcc; ///< Defined in \ref Geometry_Element.

	const struct Multiarray_Operator*const cv0_vgs_fcs; ///< Defined in \ref Geometry_Element.
	const struct Multiarray_Operator*const cv0_vgs_fcc; ///< Defined in \ref Geometry_Element.
	const struct Multiarray_Operator*const cv0_vgc_fcs; ///< Defined in \ref Geometry_Element.
	const struct Multiarray_Operator*const cv0_vgc_fcc; ///< Defined in \ref Geometry_Element.
	const struct Multiarray_Operator*const vv0_vms_fcs; ///< Defined in \ref Geometry_Element.
	const struct Multiarray_Operator*const vv0_vms_fcc; ///< Defined in \ref Geometry_Element.
	const struct Multiarray_Operator*const vv0_vmc_fcs; ///< Defined in \ref Geometry_Element.
	const struct Multiarray_Operator*const vv0_vmc_fcc; ///< Defined in \ref Geometry_Element.

	// Tensor-product sub-operators.
	const struct Multiarray_Operator*const cv0_vgs_vcs; ///< Defined in \ref Geometry_Element.
	const struct Multiarray_Operator*const cv0_vgc_vcc; ///< Defined in \ref Geometry_Element.
	const struct Multiarray_Operator*const cv0_vgs_vms; ///< Defined in \ref Geometry_Element.
	const struct Multiarray_Operator*const cv0_vgc_vmc; ///< Defined in \ref Geometry_Element.

	const struct Multiarray_Operator*const cv0_vgs_vcc; ///< Defined in \ref Geometry_Element.
	const struct Multiarray_Operator*const cv0_vgc_vcs; ///< Defined in \ref Geometry_Element.
	const struct Multiarray_Operator*const vv0_vms_vcc; ///< Defined in \ref Geometry_Element.
	const struct Multiarray_Operator*const vv0_vmc_vcs; ///< Defined in \ref Geometry_Element.
};

// Constructor/Destructor functions ********************************************************************************* //

/// \brief Constructor for a derived \ref Geometry_Element.
void constructor_derived_Geometry_Element
	(struct Element* element,     ///< \ref Geometry_Element.
	 const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Destructor for a \ref Geometry_Element.
void destructor_derived_Geometry_Element
	(struct Element* element ///< Standard.
	);

#endif // DPG__geometry_element_h__INCLUDED
