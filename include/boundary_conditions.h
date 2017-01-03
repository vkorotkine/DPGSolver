// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__boundary_conditions_h__INCLUDED
#define DPG__boundary_conditions_h__INCLUDED

extern void boundary_Riemann  (const unsigned int Nn, const unsigned int Nel, double *XYZ, double *WL, double *WOut,
                               double *WB, double *nL, const unsigned int d);
extern void boundary_SlipWall (const unsigned int Nn, const unsigned int Nel, double *WL, double *WB, double *nL,
                               const unsigned int d);

#endif // DPG__boundary_conditions_h__INCLUDED
