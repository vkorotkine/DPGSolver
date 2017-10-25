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
/// \file

#include "element_error.h"

#include "macros.h"
#include "definitions_elements.h"

#include "multiarray.h"

#include "element_operators.h"
#include "element_operators_tp.h"
#include "multiarray_operator.h"
#include "simulation.h"

// Static function declarations ************************************************************************************* //

/// \brief Constructor for a derived \ref Error_Element using the standard operators.
static void constructor_derived_Error_Element_std
	(struct Element* element_ptr, ///< Defined for \ref constructor_derived_Error_Element.
	 const struct Simulation* sim ///< Defined for \ref constructor_derived_Error_Element.
	);

/// \brief Constructor for a derived \ref Error_Element using the tensor-product of sub-element operators.
static void constructor_derived_Error_Element_tp
	(struct Element* element_ptr, ///< Defined for \ref constructor_derived_Error_Element.
	 const struct Simulation* sim ///< Defined for \ref constructor_derived_Error_Element.
	);

/// \brief Constructor for the common members of a derived \ref Error_Element.
static void constructor_derived_Error_Element_common
	(struct Element* element_ptr, ///< Defined for \ref constructor_derived_Error_Element.
	 const struct Simulation* sim ///< Defined for \ref constructor_derived_Error_Element.
	);

// Interface functions ********************************************************************************************** //

void constructor_derived_Error_Element (struct Element* element_ptr, const struct Simulation* sim)
{
	switch (element_ptr->type) {
	case LINE: case TRI: case TET: case PYR:
		constructor_derived_Error_Element_std(element_ptr,sim);
		break;
	case QUAD: case HEX: case WEDGE:
		constructor_derived_Error_Element_tp(element_ptr,sim);
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}
	constructor_derived_Error_Element_common(element_ptr,sim);
}

void destructor_derived_Error_Element (struct Element* element_ptr)
{
	struct Error_Element* element = (struct Error_Element*) element_ptr;

	destructor_const_Multiarray_Vector_d(element->w_vc[0]);
	destructor_const_Multiarray_Vector_d(element->w_vc[1]);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void constructor_derived_Error_Element_std (struct Element* element_ptr, const struct Simulation* sim)
{
	UNUSED(element_ptr);
	UNUSED(sim);
}

static void constructor_derived_Error_Element_tp (struct Element* element_ptr, const struct Simulation* sim)
{
	UNUSED(element_ptr);
	UNUSED(sim);
}

static void constructor_derived_Error_Element_common (struct Element* element_ptr, const struct Simulation* sim)
{
	struct const_Element* b_e = (struct const_Element*)element_ptr;
	struct Error_Element* e   = (struct Error_Element*)element_ptr;

	e->w_vc[0] = constructor_operators_w("vcs","vcs","H_1_P_PM0",sim->p_s_v,b_e,sim); // destructed
	e->w_vc[1] = constructor_operators_w("vcc","vcc","H_1_P_PM0",sim->p_s_v,b_e,sim); // destructed
}
