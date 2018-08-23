# Using Ctags for Code Navigation

[Exuberant Ctags](http://ctags.sourceforge.net/) can be used to generate tags which can be used to
navigate through the code.

## Generating the TAGS File

To (re-)generate the tags, make the 'tags' target in the [build_targets directory](../../build_targets):
```sh
PROJECT_ROOT/build_targets$ make tags
```

Please see the [Ctags documentation](http://ctags.sourceforge.net/ctags.html) if you are interested
in the details of the parameters passed to ctags. The c++-kinds flag is required as header files,
'\*.h', are interpreted as c++ files.

Note that using the default "projectile-regenerate-tags" in spacemacs (SPC p G) does not add tags for function
prototypes (function declarations) resulting in jumping to function declarations not being supported.

## Keybindings for Code Navigation

In spacemacs, the 'C-]' keybinding can be used to jump to function definition/declaration and struct
declarations when the cursor is positioned over the name of interest; the 'C-t' keybinding can then
be used to pop from the stack of jumps. Note that 'C-' above stands for holding the 'Ctrl' key. When
a single option is available, pressing 'C-]' will immediately jump to the tag. When multiple options
are available, a list will appear from which to choose.
