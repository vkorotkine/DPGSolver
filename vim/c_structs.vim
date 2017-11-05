" This file allows for additional syntax colouring when using vim. It can be loaded by adding the following line to the
" .vimrc file:
" autocmd BufNewFile,BufRead *DPGSolver/* source ~/FULL_PATH_TO_ROOT/DPGSolver/vim/c_structs.vim

" To see the available colour groups: `:highlight`.

syntax keyword Type lapack_int
syntax keyword Type PetscInt PetscScalar PetscReal Vec Mat KSP PC

syntax keyword Special Vector_i const_Vector_i
syntax keyword Special Vector_d const_Vector_d

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
syntax keyword Special Multiarray_Operator Operator
syntax keyword Special Multiarray_c

syntax keyword Special Intrusive_List const_Intrusive_List
syntax keyword Special Intrusive_Link const_Intrusive_Link


syntax keyword Identifier Element const_Element
syntax keyword Identifier Geometry_Element const_Geometry_Element
syntax keyword Identifier Plotting_Element const_Plotting_Element
syntax keyword Identifier Solution_Element const_Solution_Element
syntax keyword Identifier Error_Element const_Error_Element
syntax keyword Identifier DG_Solver_Element const_DG_Solver_Element

syntax keyword Identifier Volume
syntax keyword Identifier Solver_Volume
syntax keyword Identifier DG_Solver_Volume
syntax keyword Identifier Complex_DG_Solver_Volume

syntax keyword Identifier Face
syntax keyword Identifier Solver_Face
syntax keyword Identifier DG_Solver_Face
syntax keyword Identifier Complex_DG_Solver_Face
