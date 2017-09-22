// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__sum_factorization_h__INCLUDED
#define DPG__sum_factorization_h__INCLUDED

extern void   get_sf_parameters  (const unsigned int NIn0, const unsigned int NOut0, double const *const OP0,
                                  const unsigned int NIn1, const unsigned int NOut1, double const *const OP1,
                                  unsigned int *const NIn_SF, unsigned int *const NOut_SF, double const **const OP_SF,
                                  const unsigned int d, const unsigned int dim1, const unsigned int Eclass);
extern void   get_sf_parametersV (const unsigned int NIn0, const unsigned int NOut0, double const *const *const OP0,
                                  const unsigned int NIn1, const unsigned int NOut1, double const *const *const OP1,
                                  unsigned int *const NIn_SF, unsigned int *const NOut_SF, double const **const OP_SF,
                                  const unsigned int d, const unsigned int vh, const unsigned int Eclass);
extern void   get_sf_parametersF (const unsigned int NIn0, const unsigned int NOut0, double const *const *const OP0,
                                  const unsigned int NIn1, const unsigned int NOut1, double const *const *const OP1,
                                  unsigned int *const NIn_SF, unsigned int *const NOut_SF, double const **const OP_SF,
                                  const unsigned int d, const unsigned int Vf, const unsigned int Eclass);
extern void   get_sf_parametersFd (const unsigned int NIn0, const unsigned int NOut0, double const *const *const OP0,
                                   const unsigned int NIn1, const unsigned int NOut1,
                                   double const *const *const *const OP1, unsigned int *const NIn_SF,
                                   unsigned int *const NOut_SF, double const **const OP_SF, const unsigned int d,
                                   const unsigned int Vf, const unsigned int Eclass, const unsigned int dimF,
                                   const unsigned int dimD);
extern void   get_sf_parametersE (const unsigned int NIn0, const unsigned int NOut0, double const *const *const OP0,
                                  const unsigned int NIn1, const unsigned int NOut1, double const *const *const OP1,
                                  unsigned int *const NIn_SF, unsigned int *const NOut_SF, double const **const OP_SF,
                                  const unsigned int d, const unsigned int Ve, const unsigned int Eclass);
extern void   sf_operate_d       (const unsigned int NOut, const unsigned int NCols, const unsigned int NIn,
                                  const unsigned int BRowMaxIn, double *OP, double *Input, double *Output);
extern void   sf_swap_d          (double *Input, const unsigned int NRows, const unsigned int NCols,
                                  const unsigned int iBound, const unsigned int jBound, const unsigned int kBound,
                                  const unsigned int iStep, const unsigned int jStep, const unsigned int kStep);
extern void   sf_apply_d         (const double *Input, double *Output, const unsigned int NIn[3],
                                  const unsigned int NOut[3], const unsigned int NCols, double const *const *const OP,
                                  const unsigned int Diag[3], const unsigned int d);
extern double *sf_assemble_d     (const unsigned int NIn[3], const unsigned int NOut[3], const unsigned int d,
                                  double const **const BOP);

#endif // DPG__sum_factorization_h__INCLUDED
