// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/**	\file
 *	\brief Set up operators for the base element.
 *
 *	\todo update comments.
 *
 *	\section s1_suo Notation
 *
 *	Operators names are chosen according to the following template: `(0)(D){1}(r)_{2}{3}(4)_{5}{6}(7)` where elements in
 *	curly braces {} are required and those in round brackets () are optional.
 *		- (0): Special operator prefix for fused operations
 *			- op_DW_ (\todo ToBeModified: Move this to the derived element class comments where relevant)
 *		- (D): Perform a (D)ifferentiation as part of the operation.
 *		- {1}: Letter to denote the type of operation to be performed:
 *			- {T}ransform:   coefficients to coefficients
 *			- {E}valuate:    coefficients to values
 *			- {I}nterpolate: values       to values
 *			- {P}roject:     values       to coefficients
 *		- (r): Denotes the (r)eference basis (from which all operators are built).
 *		- {2/5}: Element entity
 *			- {v}olume
 *			- {f}ace
 *			- {e}dge
 *		- {3/6}: Type of the quantity for which the operator should be used.
 *			- {s}olution
 *			- {c}ubature
 *			- {g}eometry
 *			- {m}etric
 *			- {p}lotting
 *		- {4/7}: Indication for whether the operator is for {s}traight or {c}urved elements.
 */

#include "set_up_operators.h"
#include "Simulation.h"

#include <stdlib.h>
#include <stdio.h>

#include "Macros.h"

void set_up_operators (const struct Simulation*const simulation)
{
UNUSED(simulation);
/*	for (struct Element* element = simulation->element_head; element; element = element->next) {
		if (!element->tensor_product)
			;
		else
			;
	}
*/
}
