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
 *  \brief Provides the macro definitions used for c-style templating related to the volume rlhs computing functions.
 */

#if TYPE_RC == TYPE_REAL

///\{ \name Data types
#define Flux_Ref_T Flux_Ref
///\}

///\{ \name Function names
#define constructor_Flux_Ref_T constructor_Flux_Ref
#define destructor_Flux_Ref_T destructor_Flux_Ref
///\}

///\{ \name Static names
#define constructor_flux_ref_piece_T constructor_flux_ref_piece_T
///\}

#elif TYPE_RC == TYPE_COMPLEX

///\{ \name Data types
#define Flux_Ref_T Flux_Ref_c
///\}

///\{ \name Function names
#define constructor_Flux_Ref_T constructor_Flux_Ref_c
#define destructor_Flux_Ref_T destructor_Flux_Ref_c
///\}

///\{ \name Static names
#define constructor_flux_ref_piece_T constructor_flux_ref_piece_c
///\}

#endif
