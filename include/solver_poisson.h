// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__solver_poisson_h__INCLUDED
#define DPG__solver_poisson_h__INCLUDED

extern void solver_Poisson        (void);
extern void implicit_info_Poisson (void);
extern void boundary_Dirichlet    (const unsigned int Nn, const unsigned int Nel, double *XYZ, double *uL, double *uR);
extern void trace_coef            (const unsigned int Nn, const unsigned int Nel, const double *nL, double *u_avg, double *jump_u,        
                                   const unsigned int d, const char *trace_type);
extern void jacobian_flux_coef    (const unsigned int Nn, const unsigned int Nel, const double *nIn, const double *h,
                                   const unsigned int P, double *gradu_avg, double *u_jump, const unsigned int d, char *flux_type,
                                   char side);

#endif // DPG__solver_poisson_h__INCLUDED
