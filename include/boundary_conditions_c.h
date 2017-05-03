// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__boundary_conditions_c_h__INCLUDED
#define DPG__boundary_conditions_c_h__INCLUDED

#include <complex.h>

#include "boundary_conditions.h"

extern void set_BC_from_BType            (struct S_BC *const BCdata, char const *const BType);
extern void correct_XYZ_for_exact_normal (struct S_BC *const BCdata, char const *const BType);

extern void compute_boundary_values_c    (struct S_BC *const BCdata);
extern void boundary_Riemann_c           (struct S_BC *const BCdata);
extern void boundary_SlipWall_c          (struct S_BC *const BCdata);
extern void boundary_BackPressure_c      (struct S_BC *const BCdata);
extern void boundary_Total_TP_c          (struct S_BC *const BCdata);
extern void boundary_SupersonicInflow_c  (struct S_BC *const BCdata);
extern void boundary_SupersonicOutflow_c (struct S_BC *const BCdata);
extern void boundary_NoSlip_Dirichlet_c  (struct S_BC *const BCdata);
extern void boundary_NoSlip_Adiabatic_c  (struct S_BC *const BCdata);

#endif // DPG__boundary_conditions_c_h__INCLUDED
