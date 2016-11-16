// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__memory_destructors_h__INCLUDED
#define DPG__memory_destructors_h__INCLUDED

#include "S_ELEMENT.h"
#include "S_VOLUME.h"
#include "S_FACE.h"


extern void memory_destructor_E             (struct S_ELEMENT *ELEMENT);
extern void memory_destructor_V             (struct S_VOLUME *VOLUME);
extern void memory_destructor_F             (struct S_FACE *FACE);
extern void memory_destructor_L2_projection (const unsigned int EType);

#endif // DPG__memory_destructors_h__INCLUDED
