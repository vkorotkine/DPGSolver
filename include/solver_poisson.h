// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__solver_poisson_h__INCLUDED
#define DPG__solver_poisson_h__INCLUDED

extern void solver_Poisson        (void);
extern void implicit_info_Poisson (void);
extern void boundary_Dirichlet    (const unsigned int Nn, const unsigned int Nel, double *XYZ, double *uL, double *uR);
extern void trace_IP              (const unsigned int Nn, const unsigned int Nel, double *uL, double *uR, double *uNum);
extern void flux_IP               (const unsigned int Nn, const unsigned int Nel, double *uL, double *uR, double *grad_uL, double *grad_uR,
                                   double *qL, double *qR, double *h, const unsigned int P, double *nqNum, double *nL, const unsigned int d);

#endif // DPG__solver_poisson_h__INCLUDED
