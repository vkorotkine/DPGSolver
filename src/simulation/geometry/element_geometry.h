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

#ifndef DPG__element_geometry_h__INCLUDED
#define DPG__element_geometry_h__INCLUDED
/** \file
 *  \brief Provides the interface for the derived \ref Geometry_Element container and associated functions.
 */

#include "element.h"

struct Simulation;

/// \brief Container for data relating to the geometry elements.
struct Geometry_Element {
	const struct const_Element element; ///< Base \ref const_Element.

	const struct Multiarray_Operator* vv0_vv_vg[2]; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* vc0_vgc_vgc;  ///< See notation in \ref element_operators.h.

	const struct Multiarray_Operator* cv1_vg_vc[2]; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* cv1_vg_vm[2]; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* vv0_vm_vc[2]; ///< See notation in \ref element_operators.h.

	const struct Multiarray_Operator* cv0_vgs_fc[2]; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* cv0_vgc_fc[2]; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* vv0_vms_fc[2]; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* vv0_vmc_fc[2]; ///< See notation in \ref element_operators.h.

	// Blending/curved element operators
	const struct Multiarray_Operator* vv0_fv_vgc;  ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* vv0_vv_fgc;  ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* vv0_vv_fv;   ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* vv0_vv_fcc;  ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* vv0_fv_fgc;  ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* vv0_fgc_vgc; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* vv0_vgc_fgc; ///< See notation in \ref element_operators.h.

	// Non-conforming testing operators
	const struct Multiarray_Operator* cv0_vgc_fgc; ///< See notation in \ref element_operators.h.

	// Tensor-product sub-operators.
	const struct Multiarray_Operator* cv0_vg_vc[2]; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* cv0_vg_vm[2]; ///< See notation in \ref element_operators.h.

	const struct Multiarray_Operator* cv0_vgs_vcc; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* cv0_vgc_vcs; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* vv0_vms_vcc; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* vv0_vmc_vcs; ///< See notation in \ref element_operators.h.

	const struct Multiarray_Operator* vv0_vgc_vgc; ///< See notation in \ref element_operators.h.

	const struct Multiarray_Operator* cv0_vgc_vgc; ///< See notation in \ref element_operators.h.
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

#endif // DPG__element_geometry_h__INCLUDED
