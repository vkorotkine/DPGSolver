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

#ifndef DPG__cubature_operators_h__INCLUDED
#define DPG__cubature_operators_h__INCLUDED
/** \file
 *  \brief Provides the interface to functions computing the coordinates (and associated cubature weights if relevant)
 *         to be used for the computation of element operators.
 *
 *  The cubature functions provided here have the possibility of returning the cubature nodes/weights displaced to
 *  sub-regions of the reference elements for use in h-adaptive operators.
 */

/** \brief Constructor for a \ref Cubature container of the given super type.
 *  \return Standard.
 *
 *  Cubature containers can be generated for the output of the operator (which is the standard usage) as well as for the
 *  input basis (which is required for operator transformation such as from the reference to the arbitrary basis). In
 *  the case of operators requiring transformation of the output as well (such as for "cc" and "vc" operators), the
 *  cubature nodes returned must necessarily form a basis for the output coefficient space.
 */
const struct const_Cubature* constructor_const_Cubature_h
	(const int ind_io,                   ///< Index indicating whether an input/output cubature is being generated.
	 const struct Op_IO[2] op_io,        ///< \ref Op_IO.
	 const int p_op,                     /**< The polynomial order index of the associated operator (**Not the order
	                                      *   of the cubature rule**). */
	 const struct const_Element* element ///< \ref const_Element.
	 const struct Simulation* sim        ///< \ref Simulation.
	);

#endif // DPG__cubature_operators_h__INCLUDED
