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

#include "element_solver_dg.h"

#include "macros.h"
#include "definitions_elements.h"

#include "multiarray.h"

#include "element_operators.h"
#include "element_operators_tp.h"
#include "multiarray_operator.h"
#include "simulation.h"

// Static function declarations ************************************************************************************* //

/// \brief Constructor for a derived \ref DG_Solver_Element using the standard operators.
static void constructor_derived_DG_Solver_Element_std
	(struct Element* element_ptr, ///< Defined for \ref constructor_derived_DG_Solver_Element.
	 const struct Simulation* sim ///< Defined for \ref constructor_derived_DG_Solver_Element.
	);

/// \brief Constructor for a derived \ref DG_Solver_Element using the tensor-product of sub-element operators.
static void constructor_derived_DG_Solver_Element_tp
	(struct Element* element_ptr, ///< Defined for \ref constructor_derived_DG_Solver_Element.
	 const struct Simulation* sim ///< Defined for \ref constructor_derived_DG_Solver_Element.
	);

// Interface functions ********************************************************************************************** //

void constructor_derived_DG_Solver_Element (struct Element* element_ptr, const struct Simulation* sim)
{
	switch (element_ptr->type) {
	case LINE: case TRI: case TET: case PYR:
		constructor_derived_DG_Solver_Element_std(element_ptr,sim);
		break;
	case QUAD: case HEX: case WEDGE:
		constructor_derived_DG_Solver_Element_tp(element_ptr,sim);
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}
}

void destructor_derived_DG_Solver_Element (struct Element* element_ptr)
{
	struct DG_Solver_Element* dg_s_e = (struct DG_Solver_Element*) element_ptr;
	UNUSED(dg_s_e);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void constructor_derived_DG_Solver_Element_std (struct Element* element_ptr, const struct Simulation* sim)
{
	struct DG_Solver_Element* dg_s_e = (struct DG_Solver_Element*) element_ptr;
	struct const_Element* e = (struct const_Element*)element_ptr;
	UNUSED(dg_s_e);
	UNUSED(e);
	UNUSED(sim);
}

static void constructor_derived_DG_Solver_Element_tp (struct Element* element_ptr, const struct Simulation* sim)
{
	struct DG_Solver_Element* dg_s_e = (struct DG_Solver_Element*) element_ptr;
	const struct const_Element* e      = (const struct const_Element*) element_ptr;
	UNUSED(dg_s_e);
	UNUSED(e);
	UNUSED(sim);
}
