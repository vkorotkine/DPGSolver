// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__jacobian_boundary_conditions_h__INCLUDED
#define DPG__jacobian_boundary_conditions_h__INCLUDED

#include "boundary_conditions.h"

extern void compute_jacobian_boundary_values    (struct S_BC *const BCdata);
extern void jacobian_boundary_Riemann           (const unsigned int Nn, const unsigned int Nel, const double *XYZ,
                                                 const double *WL, double *WOut, double *dWdW, const double *nL,
                                                 const unsigned int d, const unsigned int Neq);
extern void jacobian_boundary_SlipWall          (const unsigned int Nn, const unsigned int Nel, const double *WL,
                                                 double *dWdW, const double *nL, const unsigned int d,
                                                 const unsigned int Neq);
extern void jacobian_boundary_BackPressure      (const unsigned int Nn, const unsigned int Nel, const double *WL,
                                                 double *dWdW, const double *nL, const unsigned int d,
                                                 const unsigned int Neq);
extern void jacobian_boundary_Total_TP          (const unsigned int Nn, const unsigned int Nel, const double *XYZ,
                                                 const double *WL, double *dWdW, const double *nL, const unsigned int d,
                                                 const unsigned int Neq);
extern void jacobian_boundary_SupersonicInflow  (const unsigned int Nn, const unsigned int Nel, const double *XYZ,
                                                 const double *WL, double *dWdW, const double *nL, const unsigned int d,
                                                 const unsigned int Neq);
extern void jacobian_boundary_SupersonicOutflow (const unsigned int Nn, const unsigned int Nel, const double *XYZ,
                                                 const double *WL, double *dWdW, const double *nL, const unsigned int d,
                                                 const unsigned int Neq);
extern void jacobian_boundary_NoSlip_Dirichlet  (struct S_BC *const BCdata);
extern void jacobian_boundary_NoSlip_Adiabatic  (struct S_BC *const BCdata);

#endif // DPG__jacobian_boundary_conditions_h__INCLUDED
