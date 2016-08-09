// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__update_FACETs_h__INCLUDED
#define DPG__update_FACETs_h__INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "Database.h"
#include "Parameters.h"
#include "Macros.h"

#include "element_functions.h"
#include "adaptation.h"
#include "setup_geometry.h"
#include "setup_normals.h"
#include "memory_constructors.h"
#include "memory_destructors.h"


extern void update_FACET_hp (void);

#endif // DPG__update_FACETs_h__INCLUDED
