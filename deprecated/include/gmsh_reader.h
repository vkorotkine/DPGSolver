// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__gmsh_reader_h__INCLUDED
#define DPG__gmsh_reader_h__INCLUDED

extern void gmsh_reader               (void);
extern void find_periodic_connections (unsigned int *Pve, unsigned int *pvePointer, const unsigned int VeMax);

#endif // DPG__gmsh_reader_h__INCLUDED
