// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__jacobian_boundary_conditions_h__INCLUDED
#define DPG__jacobian_boundary_conditions_h__INCLUDED

#include "boundary_conditions.h"

extern void compute_jacobian_boundary_values    (struct S_BC *const BCdata);
extern void jacobian_boundary_Riemann           (struct S_BC *const BCdata);
extern void jacobian_boundary_SlipWall          (struct S_BC *const BCdata);
extern void jacobian_boundary_BackPressure      (struct S_BC *const BCdata);
extern void jacobian_boundary_Total_TP          (struct S_BC *const BCdata);
extern void jacobian_boundary_SupersonicInflow  (struct S_BC *const BCdata);
extern void jacobian_boundary_SupersonicOutflow (struct S_BC *const BCdata);
extern void jacobian_boundary_NoSlip_Dirichlet  (struct S_BC *const BCdata);
extern void jacobian_boundary_NoSlip_Adiabatic  (struct S_BC *const BCdata);

#endif // DPG__jacobian_boundary_conditions_h__INCLUDED
