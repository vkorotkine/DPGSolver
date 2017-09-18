# Usage of the vim Text Editor

Several aids are provided when working with the vim text editor.

## Customized Syntax Colouring

Customized syntax colouring can be enabled by loading the provided [c_structs file](c_structs) by adding the following
line to the `.vimrc` file:
```vim
autocmd BufNewFile,BufRead *DPGSolver/* source ~/FULL_PATH_TO_ROOT/DPGSolver/vim/c_structs.vim
```

If it is desired to use alternate colour groups, those which are available can be seen with `:highlight`.

## Text Folding

Text folding can be achieved using vim by fencing code blocks between the `{{{` and `}}}` fold markers. You can then
enter `zm` to fold all and `zr` to unfold all. Folding by default upon opening files can be enabled by adding the
following line to the `.vimrc` file:
```vim
set foldmethod=marker
```

## Using the System Clipboard

Using the system clipboard instead of the vim clipboard allows for yanking and pasting across different vim sessions.
This is of particular value due to the distribution of the source across multiple directories. To use the system
clipboard, first ensure that the version of vim which you are using was compiled with the clipboard option:

1. Ensure that the version of vim which you are using was compiled with the clipboard option; you should see
`+clipboard` in the list when typing `:version` in a vim session.
2. Usage of the system clipboard is then enabled by adding the following line to the `.vimrc` file:
```vim
set clipboard=unnamed
```
