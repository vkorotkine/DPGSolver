// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__jacobian_boundary_conditions_h__INCLUDED
#define DPG__jacobian_boundary_conditions_h__INCLUDED

//extern void jacobian_boundary_Riemann  (const unsigned int Nn, const unsigned int Nel, double *XYZ, double *WL, double *WOut,
 //                                       double *dWdW, double *nL, const unsigned int d);
extern void jacobian_boundary_SlipWall (const unsigned int Nn, const unsigned int Nel, double *WL, double *dWdW, double *nL,
                                        const unsigned int d);

#endif // DPG__jacobian_boundary_conditions_h__INCLUDED
