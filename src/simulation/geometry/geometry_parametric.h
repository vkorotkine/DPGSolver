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

#ifndef DPG__geometry_parametric_h__INCLUDED
#define DPG__geometry_parametric_h__INCLUDED
/** \file
 *  \brief Provides the interface to real functions used for parametric geometry processing.
 */

#include "def_templates_type_d.h"
#include "def_templates_geometry_parametric.h"
#include "def_templates_volume_solver.h"
#include "def_templates_multiarray.h"
#include "geometry_parametric_T.h"
#include "undef_templates_type.h"
#include "undef_templates_geometry_parametric.h"
#include "undef_templates_volume_solver.h"
#include "undef_templates_multiarray.h"

/// \brief Container for function data.
struct Function_Data_GP {
	double scale; ///< Scaling parameter
};

/** \brief Return the value of the Gaussian Bump function of given differentiation degree at the input coordinate.
 *  \return See brief.
 *
 *  \note An additional parameter is provided such that the bump can be scaled by a linear function to remove its
 *        symmetry.
 */
double f_gaussian_bump
	(const double x,                              ///< x-coordinate.
	 const int diff_degree,                     ///< Differentiation degree.
	 const struct Function_Data_GP*const f_data ///< \ref Function_Data_GP.
	);

#endif // DPG__geometry_parametric_h__INCLUDED
