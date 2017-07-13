// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__implicit_GradW_h__INCLUDED
#define DPG__implicit_GradW_h__INCLUDED

extern void implicit_GradW          (void);
extern void implicit_GradW_VOLUME   (void);
extern void implicit_GradW_FACE     (void);
extern void implicit_GradW_finalize (void);

#endif // DPG__implicit_GradW_h__INCLUDED
