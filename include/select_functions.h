// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__select_functions_h__INCLUDED
#define DPG__select_functions_h__INCLUDED

#include "cubature.h"
#include "bases.h"
#include "matrix_structs.h"

extern void select_functions            (basis_tdef *basis, grad_basis_tdef *grad_basis, cubature_tdef *cubature,
                                         const unsigned int type);
extern void select_functions_basis      (basis_tdef *basis,           const unsigned int type);
extern void select_functions_grad_basis (grad_basis_tdef *grad_basis, const unsigned int type);
extern void select_functions_cubature   (cubature_tdef *cubature,     const unsigned int type);
extern void select_functions_basis_Bezier (basis_tdef *basis, const unsigned int type);

#endif // DPG__select_functions_h__INCLUDED
