// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__update_VOLUMEs_h__INCLUDED
#define DPG__update_VOLUMEs_h__INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <limits.h>

//#include "Database.h"
//#include "Parameters.h"


extern void update_VOLUME_hp       (void);
extern void update_VOLUME_list     (void);
extern void update_Vgrp            (void);
extern void update_VOLUME_Ops      (void);
extern void update_VOLUME_finalize (void);

#endif // DPG__update_VOLUMEs_h__INCLUDED
