" This file allows for additional syntax colouring when using vim. It can be loaded by adding the following line to the
" .vimrc file:
" autocmd BufNewFile,BufRead *DPGSolver/* source ~/FULL_PATH_TO_ROOT/DPGSolver/vim/c_structs.vim

" To see the available colour groups: `:highlight`.

syntax keyword Type MKL_INT
syntax keyword Type lapack_int
syntax keyword Type PetscInt PetscScalar PetscReal Vec Mat KSP PC IS

syntax keyword Special Vector_i const_Vector_i
syntax keyword Special Vector_d const_Vector_d
syntax keyword Special Vector_c const_Vector_c

syntax keyword Special Matrix_i const_Matrix_i
syntax keyword Special Matrix_d const_Matrix_d
syntax keyword Special Matrix_c const_Matrix_c

syntax keyword Special Matrix_CSR_d const_Matrix_CSR_d

syntax keyword Special Multiarray_i const_Multiarray_i
syntax keyword Special Multiarray_d const_Multiarray_d
syntax keyword Special Multiarray_c const_Multiarray_c
syntax keyword Special Multiarray_Vector_i const_Multiarray_Vector_i
syntax keyword Special Multiarray_Vector_d const_Multiarray_Vector_d
syntax keyword Special Multiarray_Matrix_d const_Multiarray_Matrix_d
syntax keyword Special Multiarray_Matrix_c const_Multiarray_Matrix_c
syntax keyword Special Multiarray_Operator Operator

syntax keyword Special Intrusive_List const_Intrusive_List
syntax keyword Special Intrusive_Link const_Intrusive_Link


syntax keyword Identifier Element const_Element
syntax keyword Identifier Basis
syntax keyword Identifier Geometry_Element
syntax keyword Identifier Plotting_Element
syntax keyword Identifier Solution_Element
syntax keyword Identifier Adaptation_Element
syntax keyword Identifier Solver_Element
syntax keyword Identifier DG_Solver_Element
syntax keyword Identifier DPG_Solver_Element

syntax keyword Identifier Volume
syntax keyword Identifier Solver_Volume_T         Solver_Volume     Solver_Volume_c
syntax keyword Identifier DG_Solver_Volume_T   DG_Solver_Volume  DG_Solver_Volume_c
syntax keyword Identifier DPG_Solver_Volume_T DPG_Solver_Volume DPG_Solver_Volume_c
syntax keyword Identifier Adaptive_Solver_Volume

syntax keyword Identifier Face
syntax keyword Identifier Solver_Face_T         Solver_Face     Solver_Face_c
syntax keyword Identifier DG_Solver_Face_T   DG_Solver_Face  DG_Solver_Face_c
syntax keyword Identifier DPG_Solver_Face_T DPG_Solver_Face DPG_Solver_Face_c
syntax keyword Identifier Adaptive_Solver_Face



syntax keyword Type Type Real Complex Index
syntax keyword Special Vector_T const_Vector_T
syntax keyword Special Vector_R const_Vector_R
syntax keyword Special Matrix_T const_Matrix_T
syntax keyword Special Matrix_R const_Matrix_R
syntax keyword Special Matrix_C const_Matrix_C
syntax keyword Special Multiarray_T const_Multiarray_T
syntax keyword Special Multiarray_R const_Multiarray_R
syntax keyword Special Multiarray_Vector_T const_Multiarray_Vector_T
syntax keyword Special Multiarray_Matrix_T const_Multiarray_Matrix_T

syntax keyword Special Matrix_CSR_T const_Matrix_CSR_T
