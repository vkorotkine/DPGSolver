// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__vertices_to_exact_geom_h__INCLUDED
#define DPG__vertices_to_exact_geom_h__INCLUDED

#include "S_VOLUME.h"


extern void   vertices_to_exact_geom (void);
extern void   Ringleb_boundary       (double *xStore, double *yStore, double qIn, double kIn, const char RinglebType);
extern double f_gaussian_bump        (const double x, const double y, const unsigned int d);
extern double f_naca_symmetric       (const double x, const double y, const unsigned int d);

#endif // DPG__vertices_to_exact_geom_h__INCLUDED
