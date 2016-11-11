// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__setup_geometry_h__INCLUDED
#define DPG__setup_geometry_h__INCLUDED

#include "S_VOLUME.h"
#include "S_FACET.h"


extern void setup_geometry  (void);
extern void setup_straight  (struct S_VOLUME *VOLUME);
extern void setup_FACET_XYZ (struct S_FACET *FACET);

#endif // DPG__setup_geometry_h__INCLUDED
