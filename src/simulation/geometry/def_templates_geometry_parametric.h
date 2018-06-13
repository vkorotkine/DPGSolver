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
 *  \brief Provides the macro definitions used for c-style templating related to the parametric geometry functions.
 */

#if TYPE_RC == TYPE_REAL

///\{ \name Function names
#define constructor_xyz_fixed_cube_parametric_T                 constructor_xyz_fixed_cube_parametric
#define constructor_xyz_cylinder_parametric_T                   constructor_xyz_cylinder_parametric
#define constructor_xyz_trigonometric_cube_parametric_T         constructor_xyz_trigonometric_cube_parametric
#define constructor_xyz_trigonometric_cube_parametric_xl_T      constructor_xyz_trigonometric_cube_parametric_xl
#define constructor_xyz_trigonometric_cube_parametric_xl_oct1_T constructor_xyz_trigonometric_cube_parametric_xl_oct1
#define constructor_xyz_joukowski_parametric_T                  constructor_xyz_joukowski_parametric
#define constructor_xyz_gaussian_bump_parametric_T              constructor_xyz_gaussian_bump_parametric
///\}

///\{ \name Static names
#define Geo_Data Geo_Data
#define get_geo_data get_geo_data
#define romberg romberg
#define f_al_gaussian_bump f_al_gaussian_bump
#define read_data_joukowski read_data_joukowski
#define read_data_gaussian_bump read_data_gaussian_bump
#define read_data_fixed_cube read_data_fixed_cube
///\}

#elif TYPE_RC == TYPE_COMPLEX

///\{ \name Function names
#define constructor_xyz_fixed_cube_parametric_T                 constructor_xyz_fixed_cube_parametric_c
#define constructor_xyz_cylinder_parametric_T                   constructor_xyz_cylinder_parametric_c
#define constructor_xyz_trigonometric_cube_parametric_T         constructor_xyz_trigonometric_cube_parametric_c
#define constructor_xyz_trigonometric_cube_parametric_xl_T      constructor_xyz_trigonometric_cube_parametric_xl_c
#define constructor_xyz_trigonometric_cube_parametric_xl_oct1_T constructor_xyz_trigonometric_cube_parametric_xl_oct1_c
#define constructor_xyz_joukowski_parametric_T                  constructor_xyz_joukowski_parametric_c
#define constructor_xyz_gaussian_bump_parametric_T              constructor_xyz_gaussian_bump_parametric_c
///\}

///\{ \name Static names
#define Geo_Data Geo_Data_c
#define get_geo_data get_geo_data_c
#define romberg romberg_c
#define f_al_gaussian_bump f_al_gaussian_bump_c
#define read_data_joukowski read_data_joukowski_c
#define read_data_gaussian_bump read_data_gaussian_bump_c
#define read_data_fixed_cube read_data_fixed_cube_c
///\}

#endif
