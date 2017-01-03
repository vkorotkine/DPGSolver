// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__setup_geometry_h__INCLUDED
#define DPG__setup_geometry_h__INCLUDED

#include "S_VOLUME.h"
#include "S_FACE.h"


extern void setup_geometry     (void);
extern void setup_straight     (struct S_VOLUME *VOLUME);
extern void setup_FACE_XYZ     (struct S_FACE *FACE);

#endif // DPG__setup_geometry_h__INCLUDED
