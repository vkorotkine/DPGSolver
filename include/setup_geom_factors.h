// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__setup_geom_factors_h__INCLUDED
#define DPG__setup_geom_factors_h__INCLUDED

#include "S_VOLUME.h"
#include "S_FACE.h"


extern void setup_geom_factors           (struct S_VOLUME *VOLUME);
extern void setup_geom_factors_highorder (struct S_FACE *FACE);

#endif // DPG__setup_geom_factors_h__INCLUDED
