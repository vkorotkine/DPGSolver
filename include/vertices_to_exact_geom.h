// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__vertices_to_exact_geom_h__INCLUDED
#define DPG__vertices_to_exact_geom_h__INCLUDED

#include "S_VOLUME.h"


extern void vertices_to_exact_geom (void);
extern void Ringleb_boundary       (double *xStore, double *yStore, double qIn, double kIn, const char RinglebType);

#endif // DPG__vertices_to_exact_geom_h__INCLUDED
