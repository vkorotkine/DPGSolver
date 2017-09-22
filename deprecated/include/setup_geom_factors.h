// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__setup_geom_factors_h__INCLUDED
#define DPG__setup_geom_factors_h__INCLUDED

#include "S_VOLUME.h"
#include "S_FACE.h"


extern void setup_geom_factors (struct S_VOLUME *VOLUME);
extern void compute_detJV      (const unsigned int Nn, double *J, double *detJV);

#endif // DPG__setup_geom_factors_h__INCLUDED
