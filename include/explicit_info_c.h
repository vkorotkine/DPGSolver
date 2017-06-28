// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__template_h__INCLUDED
#define DPG__template_h__INCLUDED

#include "Parameters.h"
#include "S_VOLUME.h"

 struct S_LOCAL_VOLUMES {
	unsigned int size;
	struct S_VOLUME *Vlist[NFMAX*NSUBFMAX];
};

extern struct S_LOCAL_VOLUMES compute_local_VOLUME_list (struct S_VOLUME const *const VOLUME);
extern bool is_VOLUME_in_local_list (struct S_VOLUME const *const VOLUME, struct S_LOCAL_VOLUMES const *const local_VOLUMEs);

#endif // DPG__template_h__INCLUDED
