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

#include "approximate_nearest_neighbor.h"

#include "macros.h"

#include "matrix.h"
#include "vector.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

const struct const_Vector_i* constructor_ann_indices
	(const struct const_Matrix_d*const nodes_b, const struct const_Matrix_d*const nodes_s)
{
	/** The implementation is based on the minimalist [ANN algorithm of Chan].
	 *
	 *  <!-- References: -->
	 *  [ANN algorithm of Chan]: ann/Chan2006_A_Minimalist's_Implementation_of_an_Approximate_Nearest_Neighbor_Algorithm_in_Fixed_Dimensions.pdf
	 */

EXIT_ADD_SUPPORT; UNUSED(nodes_b); UNUSED(nodes_s);
// Might need floating point less than operator from thesis below. If they used the conversion to unsigned long long
// however, it would be simpler to just imagine the Matrix_d*s and Matrix_ull*s and directly use Chan's code.
// https://pdfs.semanticscholar.org/fbab/cbdac7002cd147c582915f8491771e254c3f.pdf

}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
