// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__boundary_conditions_c_h__INCLUDED
#define DPG__boundary_conditions_c_h__INCLUDED

#include <complex.h>


extern void boundary_Riemann_c  (const unsigned int Nn, const unsigned int Nel, double *XYZ, double complex *WL,
                                 double complex *WOut, double complex *WB, double *nL, const unsigned int d);
extern void boundary_SlipWall_c (const unsigned int Nn, const unsigned int Nel, double complex *WL, double complex *WB,
                                 double *nL, const unsigned int d);

#endif // DPG__boundary_conditions_c_h__INCLUDED
