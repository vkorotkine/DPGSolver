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

#ifndef DPG__definitions_slepc_h__INCLUDED
#define DPG__definitions_slepc_h__INCLUDED
/** \file
 *  \brief Provides the definitions relating to SLEPc library.
 *
 *  \note Due to the \ref EPS macro having been previously defined, an alternate name for the SLEPc EPS is required.
 */

typedef struct _p_EPS* SlepcEPS; ///< Typedef to avoid conflict with \ref EPS.

#endif // DPG__definitions_slepc_h__INCLUDED
