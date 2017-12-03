@page BuildWithClang Building With Clang

Basically, all you need to do is tell CMake that you want to use a different
compiler, i.e.

```
cmake -DCMAKE_CXX_COMPILER=/usr/bin/clang++ path/to/nuto_src
```
