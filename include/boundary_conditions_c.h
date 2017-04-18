// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__boundary_conditions_c_h__INCLUDED
#define DPG__boundary_conditions_c_h__INCLUDED

#include <complex.h>

#include "boundary_conditions.h"

extern void boundary_Riemann_c           (const unsigned int Nn, const unsigned int Nel, const double *const XYZ,
                                          const double complex *const WL, double complex *const WOut,
                                          double complex *const WB, const double *const nL, const unsigned int d);
extern void boundary_SlipWall_c          (const unsigned int Nn, const unsigned int Nel, const double complex *const WL,
                                          double complex *const WB, const double *const nL, const unsigned int d);
extern void boundary_BackPressure_c      (const unsigned int Nn, const unsigned int Nel, const double complex *const WL,
                                          double complex *const WB, const double *const nL, const unsigned int d,
                                          const unsigned int Neq);
extern void boundary_Total_TP_c          (const unsigned int Nn, const unsigned int Nel, const double *const XYZ,
                                          const double complex *const WL, double complex *const WB,
                                          const double *const nL, const unsigned int d, const unsigned int Nvar);
extern void boundary_SupersonicInflow_c  (const unsigned int Nn, const unsigned int Nel, const double *const XYZ,
                                          const double complex *const WL, double complex *const WB,
                                          const double *const nL, const unsigned int d, const unsigned int Nvar);
extern void boundary_SupersonicOutflow_c (const unsigned int Nn, const unsigned int Nel, const double *const XYZ,
                                          const double complex *const WL, double complex *const WB,
                                          const double *const nL, const unsigned int d, const unsigned int Nvar);
extern void boundary_NoSlip_Dirichlet_c  (struct S_BC *const BCdata);
extern void boundary_NoSlip_Adiabatic_c  (struct S_BC *const BCdata);

#endif // DPG__boundary_conditions_c_h__INCLUDED
