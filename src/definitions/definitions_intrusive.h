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
/** \file
 *  \brief Provides the definitions relating to the intrusive lists base/derived classes.
 */

///\{ \name Invalid list name
#define IL_INVALID 0
///\}

///\{ \name The intrusive list category names. Used when multiple lists must be updated together.
#define IL_BASE                1000
#define IL_SOLVER              1010
#define IL_SOLVER_DG           1011
#define IL_SOLVER_DPG          1012
#define IL_SOLVER_DG_COMPLEX   1021
#define IL_SOLVER_DPG_COMPLEX  1022
///\}

///\{ \name The \ref Element list names.
#define IL_ELEMENT             100
#define IL_ELEMENT_GEOMETRY    101
#define IL_ELEMENT_PLOTTING    102
#define IL_ELEMENT_SOLUTION    103
#define IL_ELEMENT_ERROR       104
#define IL_ELEMENT_SOLVER      105
#define IL_ELEMENT_SOLVER_DG   106
#define IL_ELEMENT_SOLVER_DPG  107
///\}

///\{ \name The \ref Volume list names.
#define IL_VOLUME                    200
#define IL_VOLUME_SOLVER             210
#define IL_VOLUME_SOLVER_DG          211
#define IL_VOLUME_SOLVER_DPG         212
#define IL_VOLUME_SOLVER_DG_COMPLEX  221
#define IL_VOLUME_SOLVER_DPG_COMPLEX 222
///\}

///\{ \name The \ref Face list names.
#define IL_FACE                     300
#define IL_FACE_SOLVER              310
#define IL_FACE_SOLVER_DG           311
#define IL_FACE_SOLVER_DPG          312
#define IL_FACE_SOLVER_DG_COMPLEX   321
#define IL_FACE_SOLVER_DPG_COMPLEX  322
///\}

#endif // DPG__definitions_intrusive_h__INCLUDED
