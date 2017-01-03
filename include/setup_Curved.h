// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__setup_Curved_h__INCLUDED
#define DPG__setup_Curved_h__INCLUDED

#include "S_VOLUME.h"


extern void setup_Curved                (struct S_VOLUME *VOLUME);
extern void setup_Curved_vertices       (struct S_VOLUME *VOLUME);
extern void compute_plane               (const double *XYZ1, const double *XYZ2, const double *XYZ3, double *n, double *d_p);
extern void compute_normal_displacement (const unsigned int Nn, const unsigned int curved_normal, const double *XYZ_S,
                                         const double *normals, double *XYZ_CmS, const unsigned int BC);
extern void get_abc_ellipse             (const unsigned int Nn, double *XYZ, double *abc);

#endif // DPG__setup_Curved_h__INCLUDED
