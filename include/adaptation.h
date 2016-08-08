// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__adaptation_h__INCLUDED
#define DPG__adaptation_h__INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "Database.h"
#include "Macros.h"

#include "element_functions.h"
#include "array_sort.h"
#include "update_VOLUMEs.h"
#include "update_FACETs.h"
#include "memory_free.h"


extern void         adapt_hp         (void);
extern void         mesh_update      (void);
extern void         mesh_to_level    (const unsigned int level);

extern void         get_PS_range     (unsigned int *PSMin, unsigned int *PSMax);
extern void         get_Pb_range     (const unsigned int P, unsigned int *PbMin, unsigned int *PbMax);
extern void         get_vh_range     (const struct S_VOLUME *VOLUME, unsigned int *vhMin, unsigned int *vhMax);
extern void         get_fh_range     (const struct S_VOLUME *VOLUME, const unsigned int f, unsigned int *fhMin, unsigned int *fhMax);
extern unsigned int get_VOLUMEc_type (const unsigned int VType, const unsigned int vh);
extern unsigned int get_IndEhref     (const unsigned int VType, const unsigned int vh);

#endif // DPG__adaptation_h__INCLUDED
