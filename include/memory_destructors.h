// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__memory_destructors_h__INCLUDED
#define DPG__memory_destructors_h__INCLUDED

#include <stdlib.h>

#include "Database.h"
#include "Parameters.h"

#include "array_free.h"
#include "element_functions.h"
#include "adaptation.h"


extern void memory_destructor_E             (struct S_ELEMENT *ELEMENT);
extern void memory_destructor_V             (struct S_VOLUME *VOLUME);
extern void memory_destructor_F             (struct S_FACET *FACET);
extern void memory_destructor_L2_projection (const unsigned int EType);

#endif // DPG__memory_destructors_h__INCLUDED
