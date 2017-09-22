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

#ifndef DPG__definitions_intrusive_h__INCLUDED
#define DPG__definitions_intrusive_h__INCLUDED
/**	\file
 *	\brief Provides the definitions relating to the intrusive lists base/derived classes.
 */

///\{ \name The \ref Element list names.
#define IL_ELEMENT          100
#define IL_GEOMETRY_ELEMENT 101
#define IL_SOLVER_ELEMENT   102
///\}

///\{ \name The \ref Volume list names.
#define IL_VOLUME        200
#define IL_SOLVER_VOLUME 201
///\}

///\{ \name The \ref Face list names.
#define IL_FACE        300
#define IL_SOLVER_FACE 301
///\}

#endif // DPG__definitions_intrusive_h__INCLUDED
