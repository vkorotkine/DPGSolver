// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__template_h__INCLUDED
#define DPG__template_h__INCLUDED

#include "Parameters.h"
#include "S_VOLUME.h"
#include "S_FACE.h"

 struct S_LOCAL_MESH_ELEMENTS {
	unsigned int NVl, NFl;
	struct S_VOLUME *Vlist[NFMAX*NSUBFMAX+1];
	struct S_FACE   *Flist[NFMAX*NSUBFMAX];
};

extern struct S_LOCAL_MESH_ELEMENTS compute_local_ELEMENT_list (struct S_VOLUME const *const VOLUME, char const entity);
extern bool   is_VOLUME_in_local_list (struct S_VOLUME const *const VOLUME,
                                       struct S_LOCAL_MESH_ELEMENTS const *const local_ELEMENTs);
extern bool   is_FACE_in_local_list (struct S_FACE const *const FACE,
                                     struct S_LOCAL_MESH_ELEMENTS const *const local_ELEMENTs);

#endif // DPG__template_h__INCLUDED
