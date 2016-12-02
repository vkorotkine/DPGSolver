// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__setup_Curved_h__INCLUDED
#define DPG__setup_Curved_h__INCLUDED

#include "S_VOLUME.h"


extern void compute_plane(const double *XYZ1, const double *XYZ2, const double *XYZ3, double *n, double *d_p);
extern void setup_Curved (struct S_VOLUME *VOLUME);

#endif // DPG__setup_Curved_h__INCLUDED
