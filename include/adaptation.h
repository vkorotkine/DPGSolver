// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__adaptation_h__INCLUDED
#define DPG__adaptation_h__INCLUDED

#include <stddef.h>

#include "S_VOLUME.h"


extern void         adapt_hp          (void);
extern void         mesh_update       (void);
extern void         mesh_to_level     (const unsigned int level);
extern void         mesh_to_order     (const unsigned int order);
extern void         mesh_h_adapt      (const unsigned int Nadapt, const char h_adapt_type);
extern void         ensure_1irregular (unsigned int *hp_update);

extern void         get_PS_range     (size_t *const PSMin, size_t *const PSMax);
extern void         get_Pb_range     (unsigned int const P, size_t *const PbMin, size_t *const PbMax);
extern void         get_vh_range     (const struct S_VOLUME *VOLUME, unsigned int *vhMin, unsigned int *vhMax);
extern void         get_fh_range     (const struct S_VOLUME *VOLUME, const unsigned int f, unsigned int *fhMin, unsigned int *fhMax);
extern unsigned int get_VOLUMEc_type (const unsigned int VType, const unsigned int vh);
extern unsigned int get_IndEhref     (const unsigned int VType, const unsigned int vh);

#endif // DPG__adaptation_h__INCLUDED
