// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__solver_Poisson_c_h__INCLUDED
#define DPG__solver_Poisson_c_h__INCLUDED

extern void compute_What_VOLUME_c (void);
extern void compute_What_FACE_c   (void);
extern void correct_collocated_for_symmetry_c (void);

#endif // DPG__solver_Poisson_c_h__INCLUDED
