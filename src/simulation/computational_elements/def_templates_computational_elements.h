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
 *  \brief Provides the macro definitions used for c-style templating related to the computational elements.
 */

#if TYPE_RC == TYPE_REAL

///\{ \name Function names
#define constructor_derived_computational_elements_T constructor_derived_computational_elements
#define destructor_derived_computational_elements_T  destructor_derived_computational_elements
///\}

#define Derived_Comp_Elements_Info                 Derived_Comp_Elements_Info
#define get_c_Derived_Comp_Elements_Info           get_c_Derived_Comp_Elements_Info
#define get_d_Derived_Comp_Elements_Info           get_d_Derived_Comp_Elements_Info
#define update_computational_element_list_pointers update_computational_element_list_pointers
#define update_volume_list_pointers                update_volume_list_pointers
#define update_face_list_pointers                  update_face_list_pointers
#define constructor_base_Intrusive_Link            constructor_base_Intrusive_Link
#define get_list_category                          get_list_category
#define update_volume_pointers                     update_volume_pointers
#define update_face_pointers                       update_face_pointers

#elif TYPE_RC == TYPE_COMPLEX

///\{ \name Function names
#define constructor_derived_computational_elements_T constructor_derived_computational_elements_c
#define destructor_derived_computational_elements_T  destructor_derived_computational_elements_c
///\}

#define Derived_Comp_Elements_Info                 Derived_Comp_Elements_Info_c
#define get_c_Derived_Comp_Elements_Info           get_c_Derived_Comp_Elements_Info_c
#define get_d_Derived_Comp_Elements_Info           get_d_Derived_Comp_Elements_Info_c
#define update_computational_element_list_pointers update_computational_element_list_pointers_c
#define update_volume_list_pointers                update_volume_list_pointers_c
#define update_face_list_pointers                  update_face_list_pointers_c
#define constructor_base_Intrusive_Link            constructor_base_Intrusive_Link_c
#define get_list_category                          get_list_category_c
#define update_volume_pointers                     update_volume_pointers_c
#define update_face_pointers                       update_face_pointers_c

#endif
