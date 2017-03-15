@page CodingStyle Coding Style

# use clang-format

Clang-format is a tool that automatically beautifies the source code. The settings are defined in `.clang-format` (and we can discuss them - have a look [at this](https://clangformat.com/) for some instructions). Use it in the command line with

~~~
clang-format -i somefile
~~~

However, every decent editor should support clang format, below are some instructions. 

### vim

https://github.com/rhysd/vim-clang-format

### clion

1) via the [ClangFormatIJ](https://plugins.jetbrains.com/plugin/8396-clangformatij) plugin (goto settings, plugins, install it)
2) by defining clang-format as an [external tool](http://stackoverflow.com/questions/34648255/using-clang-format-in-clion)

It is certainly convenient to add a keybinding for it (settings, keybinding, search for "clang" and assign a key

### QtCreator

(not much experience myself, but I found promising instructions [here](http://stackoverflow.com/a/40174996))


# other

- member variables start with `m`
- add Doxygen comments in headers, not in .cpp
- as `const` as possible
- as `private`/`protected` as possible
 
Example code: look in a file
 
In addition, [read this](http://www.se.rit.edu/~tabeec/RIT_441/Resources_files/How%20To%20Write%20Unmaintainable%20Code.pdf) and do exactly the opposite.


