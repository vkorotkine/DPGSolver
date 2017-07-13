// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__update_VOLUMEs_h__INCLUDED
#define DPG__update_VOLUMEs_h__INCLUDED

#include "S_VOLUME.h"


extern void update_VOLUME_hp       (void);
extern void update_VOLUME_list     (void);
extern void update_Vgrp            (void);
extern void compute_inverse_mass   (struct S_VOLUME *VOLUME);
extern void update_VOLUME_Ops      (void);
extern void update_VOLUME_finalize (void);
extern void update_memory_VOLUMEs  (void);

#endif // DPG__update_VOLUMEs_h__INCLUDED
