// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__solver_Poisson_h__INCLUDED
#define DPG__solver_Poisson_h__INCLUDED

#include <stdbool.h>

extern void solver_Poisson        (bool PrintEnabled);
extern void implicit_info_Poisson (void);
extern void project_to_sphere     (const unsigned int Nn, double *XYZIn, double *XYZOut, const unsigned int BCcurved);
extern void boundary_Poisson      (const unsigned int Nn, const unsigned int Nel, double *XYZ, double *normals,
                                   double *uL, double *uR, double *graduL, double *graduR, const unsigned int BC,
                                   const unsigned int BCcurved);
extern void jacobian_flux_coef    (const unsigned int Nn, const unsigned int Nel, const double *nIn, const double *h,
                                   const unsigned int P, double *gradu_avg, double *u_jump, const unsigned int d,
                                   const unsigned int flux_type);

#endif // DPG__solver_Poisson_h__INCLUDED
