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

#ifndef DPG__definitions_core_h__INCLUDED
#define DPG__definitions_core_h__INCLUDED
/** \file
 *  \brief Provides the definitions of core constants.
 */

/// \todo Potentially remove DMAX.
///\{ \name The maximum number of dimensions
#define DMAX 3
///\}

///\{ \name Supported adaptation options.
#define ADAPT_0  0
#define ADAPT_P  1
#define ADAPT_H  2
#define ADAPT_HP 3
///\}

#define DIM 3 ///< 'DIM'ension of this build. Must be passed as a CMake -D parameter.

#define NVAR_EULER (DIM+2) ///< 'N'umber of 'VAR'iables for the Euler/Navier-Stokes equations.
#define NEQ_EULER  (DIM+2) ///< 'N'umber of 'EQ'uations for the Euler/Navier-Stokes equations.


///\{ \name Macro used to form arrays of size DIM.
#if DIM == 1
	#define ARRAY_DIM(a,b,c) { (a), }
#elif DIM == 2
	#define ARRAY_DIM(a,b,c) { (a), (b), }
#elif DIM == 3
	#define ARRAY_DIM(a,b,c) { (a), (b), (c), }
#endif
///\}

#endif // DPG__definitions_core_h__INCLUDED
