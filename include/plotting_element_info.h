// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__plotting_element_info_h__INCLUDED
#define DPG__plotting_element_info_h__INCLUDED

extern void plotting_element_info (double **rst, unsigned int **connect, unsigned int **types, unsigned int **connectE,
                                   unsigned int *Nn, unsigned int *NE, const unsigned int P, const unsigned int typeIn);

#endif // DPG__plotting_element_info_h__INCLUDED
