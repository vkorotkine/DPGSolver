// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__memory_free_h__INCLUDED
#define DPG__memory_free_h__INCLUDED

#include <stdlib.h>
#include <stdio.h>
 
#include "Database.h"
#include "Parameters.h"

#include "array_free.h"
#include "memory_destructors.h"
#include "adaptation.h"


extern void memory_free          (void);
extern void memory_free_children (void);

#endif // DPG__memory_free_h__INCLUDED
