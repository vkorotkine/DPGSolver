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

#ifndef DPG__definitions_h_ref_h__INCLUDED
#define DPG__definitions_h_ref_h__INCLUDED
/** \file
 *  \brief Provides the definitions of h-refinement related constants.
 *
 *  The notation used here to define the various h-refinement related elements and sub-elements is: H_[0][1]_[2](3)[4]
 *  where all parameters are required:
 *  - [0]: Master element type (i.e. type of the element on which the h-refinement is performed).
 *  - [1]: Number of sub-elements. A value of 1 indicates no h-refinement (i.e. the standard unrefined element).
 *  - [2]: The computational element indicator in relation to the element type. Options:
 *  	- [V]olume
 *  	- [F]ace
 *  	- [E]dge
 *  - (3): The index associated with the parent computational element (present before the h-refinement).
 *  - [4]: The index associated with the child computational element(s) (present after the h-refinement).
 *
 *  The values of the constants are used as indices for Multiarray_Matrix_d\* operators and should consequently not be
 *  modified.
 *
 *  \todo Ensure that all magic numbers are replaced with the appropriate constants defined here.
 *  \todo Add reference to figures accompanying the definitions.
 */

// LINE element ***************************************************************************************************** //

///\{ \name LINE element volumes.
#define H_LINE1_V0 0 ///< Vertice(s): [0,2].

#define H_LINE2_V0 1 ///< Vertice(s): [0,1].
#define H_LINE2_V1 2 ///< Vertice(s): [1,2].
///\}

///\{ \name LINE element faces.
#define H_LINE1_F0 0 ///< Vertice(s): [0].
#define H_LINE1_F1 1 ///< Vertice(s): [1].
///\}


// TRI element ****************************************************************************************************** //

///\{ \name TRI element volumes.
#define H_TRI1_V0 0 ///< Vertice(s): [0,2,5].

#define H_TRI4_V0 1 ///< Vertice(s): [0,1,3].
#define H_TRI4_V1 2 ///< Vertice(s): [1,2,4].
#define H_TRI4_V2 3 ///< Vertice(s): [3,4,5].
#define H_TRI4_V3 4 ///< Vertice(s): [4,3,2].
///\}

///\{ \name TRI element faces.
#define H_TRI1_F0 0 ///< Vertice(s): [2,5].
#define H_TRI1_F1 1 ///< Vertice(s): [0,5].
#define H_TRI1_F2 2 ///< Vertice(s): [0,2].

#define H_TRI4_F00 3 ///< Vertice(s): [2,4].
#define H_TRI4_F01 4 ///< Vertice(s): [4,5].
#define H_TRI4_F10 5 ///< Vertice(s): [0,3].
#define H_TRI4_F11 6 ///< Vertice(s): [3,5].
#define H_TRI4_F20 7 ///< Vertice(s): [0,1].
#define H_TRI4_F21 8 ///< Vertice(s): [1,2].
///\}


// QUAD element ***************************************************************************************************** //

///\{ \name QUAD element volumes.
#define H_QUAD1_V0 0 ///< Vertice(s): [0,2,6,8].

#define H_QUAD4_V0 1 ///< Vertice(s): [0,1,3,4].
#define H_QUAD4_V1 2 ///< Vertice(s): [1,2,4,5].
#define H_QUAD4_V2 3 ///< Vertice(s): [3,4,6,7].
#define H_QUAD4_V3 4 ///< Vertice(s): [4,5,7,8].
///\}

///\{ \name QUAD element faces.
#define H_QUAD1_F0 0 ///< Vertice(s): [0,6].
#define H_QUAD1_F1 1 ///< Vertice(s): [2,8].
#define H_QUAD1_F2 2 ///< Vertice(s): [0,2].
#define H_QUAD1_F3 3 ///< Vertice(s): [6,8].

#define H_QUAD4_F00 4  ///< Vertice(s): [0,3].
#define H_QUAD4_F01 5  ///< Vertice(s): [3,6].
#define H_QUAD4_F10 6  ///< Vertice(s): [2,5].
#define H_QUAD4_F11 7  ///< Vertice(s): [5,8].
#define H_QUAD4_F20 8  ///< Vertice(s): [0,1].
#define H_QUAD4_F21 9  ///< Vertice(s): [1,2].
#define H_QUAD4_F30 10 ///< Vertice(s): [6,7].
#define H_QUAD4_F31 11 ///< Vertice(s): [7,8].
///\}


// HEX element ****************************************************************************************************** //

///\{ \name HEX element volumes.
#define H_HEX1_V0 0 ///< Vertice(s): [ 0, 2, 6, 8,18,20,24,26].

#define H_HEX8_V0 1 ///< Vertice(s): [ 0, 1, 3, 4, 9,10,12,13].
#define H_HEX8_V1 2 ///< Vertice(s): [ 1, 2, 4, 5,10,11,13,14].
#define H_HEX8_V2 3 ///< Vertice(s): [ 3, 4, 6, 7,12,13,15,16].
#define H_HEX8_V3 4 ///< Vertice(s): [ 4, 5, 7, 8,13,14,16,17].
#define H_HEX8_V4 5 ///< Vertice(s): [ 9,10,12,13,18,19,21,22].
#define H_HEX8_V5 6 ///< Vertice(s): [10,11,13,14,19,20,22,23].
#define H_HEX8_V6 7 ///< Vertice(s): [12,13,15,16,21,22,24,25].
#define H_HEX8_V7 8 ///< Vertice(s): [13,14,16,17,22,23,25,26].
///\}

///\{ \name HEX element faces.
#define H_HEX1_F0 0 ///< vertice(s): [ 0, 6,18,24].
#define H_HEX1_F1 1 ///< vertice(s): [ 2, 8,20,26].
#define H_HEX1_F2 2 ///< vertice(s): [ 0, 2,18,20].
#define H_HEX1_F3 3 ///< vertice(s): [ 6, 8,24,26].
#define H_HEX1_F4 4 ///< vertice(s): [ 0, 2, 6, 8].
#define H_HEX1_F5 5 ///< vertice(s): [18,20,24,26].

#define H_HEX8_F00 6  ///< vertice(s): [ 0, 3, 9,12].
#define H_HEX8_F01 7  ///< vertice(s): [ 3, 6,12,15].
#define H_HEX8_F02 8  ///< vertice(s): [ 9,12,18,21].
#define H_HEX8_F03 9  ///< vertice(s): [12,15,21,24].
#define H_HEX8_F10 10 ///< vertice(s): [ 2, 5,11,14].
#define H_HEX8_F11 11 ///< vertice(s): [ 5, 8,14,17].
#define H_HEX8_F12 12 ///< vertice(s): [11,14,20,23].
#define H_HEX8_F13 13 ///< vertice(s): [14,17,23,26].
#define H_HEX8_F20 14 ///< vertice(s): [ 0, 1, 9,10].
#define H_HEX8_F21 15 ///< vertice(s): [ 1, 2,10,11].
#define H_HEX8_F22 16 ///< vertice(s): [ 9,10,18,19].
#define H_HEX8_F23 17 ///< vertice(s): [10,11,19,20].
#define H_HEX8_F30 18 ///< vertice(s): [ 6, 7,15,16].
#define H_HEX8_F31 19 ///< vertice(s): [ 7, 8,16,17].
#define H_HEX8_F32 20 ///< vertice(s): [15,16,24,25].
#define H_HEX8_F33 21 ///< vertice(s): [16,17,25,26].
#define H_HEX8_F40 22 ///< vertice(s): [ 0, 1, 3, 4].
#define H_HEX8_F41 23 ///< vertice(s): [ 1, 2, 4, 5].
#define H_HEX8_F42 24 ///< vertice(s): [ 3, 4, 6, 7].
#define H_HEX8_F43 25 ///< vertice(s): [ 4, 5, 7, 8].
#define H_HEX8_F50 26 ///< vertice(s): [18,19,21,22].
#define H_HEX8_F51 27 ///< vertice(s): [19,20,22,23].
#define H_HEX8_F52 28 ///< vertice(s): [21,22,24,25].
#define H_HEX8_F53 29 ///< vertice(s): [22,23,25,26].
///\}


// WEDGE element **************************************************************************************************** //

///\{ \name WEDGE element volumes.
#define H_WEDGE1_V0 0 ///< Vertice(s): [ 0, 2, 5,12,14,17].

#define H_WEDGE8_V0 1 ///< Vertice(s): [ 0, 1, 3, 6, 7, 9].
#define H_WEDGE8_V1 2 ///< Vertice(s): [ 1, 2, 4, 7, 8,10].
#define H_WEDGE8_V2 3 ///< Vertice(s): [ 3, 4, 5, 9,10,11].
#define H_WEDGE8_V3 4 ///< Vertice(s): [ 4, 3, 1,10, 9, 7].
#define H_WEDGE8_V4 5 ///< Vertice(s): [ 6, 7, 9,12,13,15].
#define H_WEDGE8_V5 6 ///< Vertice(s): [ 7, 8,10,13,14,16].
#define H_WEDGE8_V6 7 ///< Vertice(s): [ 9,10,11,15,16,17].
#define H_WEDGE8_V7 8 ///< Vertice(s): [10, 9, 7,16,15,13].
///\}

///\{ \name WEDGE element faces.
#define H_WEDGE1_F0 0 ///< vertice(s): [ 2, 5,14,17].
#define H_WEDGE1_F1 1 ///< vertice(s): [ 0, 5,12,17].
#define H_WEDGE1_F2 2 ///< vertice(s): [ 0, 2,12,14].
#define H_WEDGE1_F3 3 ///< vertice(s): [ 0, 2, 5].
#define H_WEDGE1_F4 4 ///< vertice(s): [12,14,17].

#define H_WEDGE8_F00 5  ///< vertice(s): [ 2, 4, 8,10].
#define H_WEDGE8_F01 6  ///< vertice(s): [ 4, 5,10,11].
#define H_WEDGE8_F02 7  ///< vertice(s): [ 8,10,14,16].
#define H_WEDGE8_F03 8  ///< vertice(s): [10,11,16,17].
#define H_WEDGE8_F10 9  ///< vertice(s): [ 0, 3, 6, 9].
#define H_WEDGE8_F11 10 ///< vertice(s): [ 3, 5, 9,11].
#define H_WEDGE8_F12 11 ///< vertice(s): [ 6, 9,12,15].
#define H_WEDGE8_F13 12 ///< vertice(s): [ 9,11,15,17].
#define H_WEDGE8_F20 13 ///< vertice(s): [ 0, 1, 6, 7].
#define H_WEDGE8_F21 14 ///< vertice(s): [ 1, 2, 7, 8].
#define H_WEDGE8_F22 15 ///< vertice(s): [ 6, 7,12,13].
#define H_WEDGE8_F23 16 ///< vertice(s): [ 7, 8,13,14].
#define H_WEDGE8_F30 17 ///< vertice(s): [ 0, 1, 3].
#define H_WEDGE8_F31 18 ///< vertice(s): [ 1, 2, 4].
#define H_WEDGE8_F32 19 ///< vertice(s): [ 3, 4, 5].
#define H_WEDGE8_F33 20 ///< vertice(s): [ 4, 3, 1].
#define H_WEDGE8_F40 21 ///< vertice(s): [12,13,15].
#define H_WEDGE8_F41 22 ///< vertice(s): [13,14,16].
#define H_WEDGE8_F42 23 ///< vertice(s): [15,16,17].
#define H_WEDGE8_F43 24 ///< vertice(s): [16,15,13].
///\}

#endif // DPG__definitions_h_ref_h__INCLUDED
