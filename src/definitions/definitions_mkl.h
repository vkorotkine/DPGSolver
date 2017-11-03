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

#ifndef DPG__definitions_mkl_h__INCLUDED
#define DPG__definitions_mkl_h__INCLUDED
/** \file
 *  \brief Provides the definitions relating to MKL library function calls.
 */

///\{ \name MKL layout parameters.
#define CBRM CblasRowMajor
#define CBCM CblasColMajor
///\}

///\{ \name MKL transpose parameters.
#define CBT  CblasTrans
#define CBNT CblasNoTrans
///\}

///\{ \name Long ugly name indicating that fortran 1-based indexing should be used.
#define ADD_ONE_BASED_INDEXING 1
///\}

/**\{ \name Redefine mkl complex types for compatibility with c99 standard types.
 *    This was done following the suggestion in
 *    [this][https://software.intel.com/en-us/forums/intel-math-kernel-library/topic/285810] intel post.
 */
#define MKL_Complex16 double complex
///\}

#endif // DPG__definitions_mkl_h__INCLUDED
