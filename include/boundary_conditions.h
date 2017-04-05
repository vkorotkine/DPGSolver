// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__boundary_conditions_h__INCLUDED
#define DPG__boundary_conditions_h__INCLUDED

extern void get_boundary_values        (const double X, const double Y, double *rho, double *u, double *v, double *w,
                                         double *p);
extern void boundary_Riemann           (const unsigned int Nn, const unsigned int Nel, const double *XYZ,
                                        const double *WL, double *WOut, double *WB, const double *nL,
                                        const unsigned int d);
extern void boundary_SlipWall          (const unsigned int Nn, const unsigned int Nel, const double *WL, double *WB,
                                        const double *nL, const unsigned int d);
extern void boundary_BackPressure      (const unsigned int Nn, const unsigned int Nel, const double *WL, double *WB,
                                        const double *nL, const unsigned int d, const unsigned int Neq);
extern void boundary_Total_TP          (const unsigned int Nn, const unsigned int Nel, const double *XYZ,
                                        const double *WL, double *WB, const double *nL, const unsigned int d,
                                        const unsigned int Nvar);
extern void boundary_SupersonicInflow  (const unsigned int Nn, const unsigned int Nel, const double *XYZ,
                                        const double *WL, double *WB, const double *nL, const unsigned int d,
                                        const unsigned int Nvar);
extern void boundary_SupersonicOutflow (const unsigned int Nn, const unsigned int Nel, const double *XYZ,
                                        const double *WL, double *WB, const double *nL, const unsigned int d,
                                        const unsigned int Nvar);
extern void boundary_NoSlip_Dirichlet  (const unsigned int Nn, const unsigned int Nel, const double *XYZ,
                                        const double *WL, double *WB, const double *nL, const unsigned int d,
                                        const unsigned int Nvar);
extern void boundary_NoSlip_Adiabatic  (const unsigned int Nn, const unsigned int Nel, const double *XYZ,
                                        const double *WL, double *WB, const double *nL, const unsigned int d,
                                        const unsigned int Nvar);

#endif // DPG__boundary_conditions_h__INCLUDED
