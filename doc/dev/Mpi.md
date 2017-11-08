@page Mpi Mpi concept

# General ideas

- Domain: part of the mesh that is available in one process. For 4 processes, the mesh is divided into 4 domains.
- apart from the `Solve`, each process runs like it would run in a serial simulation
    - defining constraings
    - build _local_ dof numbering that results in a `J/K` block structure
        - `J`: independent dof numbers - [0 .. `numIndependentDofs` - 1]. 
        - `K`: dependent dof numbers - [`numIndependentDofs` - `numDofs` - 1].
    - assemble _Gradient_(`R`) and _Hessians_(`K,C,M`)
    - depending on the time integration scheme, do calculations with _Gradient_ and _Hessian_, e.g. `Hessian = K + C * x + M * m`
    - apply constraint matrix
    - [here would be the solve, this is done differently, see blow]
    - postprocessing

# `MpiSolver`

The `MpiSolver` is only interested in the independent dofs. It maps each _local_ independent dof for a process to a _global_ dof number. Each process has a _local_ independent dof numbering, starting from 0. Neighboring domains share nodes. Shared nodes must be mapped to the same _global_ dof number. 

Example:

~~~
Simple 1D mesh

Domain 0:                                      ¦   Domain 1:
-----------------------------------------------¦-----------------------------------
- 3 independent dofs (0..2)                    ¦    - 3 independent dofs (0..2)
- 1 dependent dof (3)                          ¦
                                               ¦
                    ||>|-----|-----|-----|     ¦    |------|------| 
                                               ¦
local dof numbers      3     2     0     1     ¦    0      1      2
                                               ¦
global dof numbers           0     1     2     ¦    2      3      4     
~~~

Process 0 (domain 0) shares a node with process 1 (domain 1). In process 0, it has the _local_ dof number 1. In process 1, it has the _local_ dof number 0. The `MpiSolver` now has to provide a mapping, such that

~~~
GlobalDofNumber(process 0, 1) == GlobalDofNumber(process 1, 0) = 2;
~~~
We will call this mapping **LocalToGlobalMapping** and its calculation probably the main issue. 

# Calculate **LocalToGlobalMapping**

## Hardcode

- (This is what Bernd did in [here](https://github.com/nutofem/nuto/blob/2bb2e0ed78d6e8187b7cc0b65d15e6e0275f99a9/applications/custom/TestClasses/2D_Test.cpp))
- Build it manually with knowledge of your mesh files.
- mainly for test purposes.

## Use a mesh file

- The mesh file contains a node numbering for each domain. The nodes are numbered such that shared nodes have the same node number. Problems arise when constraints require a renumbering. This renumbering has to be communicated with other processes.
- This file should _not_ contain constraints.
- ...

## Calculate internally with 3D data structure

- Have something like a 3D subbox structure.
- ...


# Implementation

Say, we have this LocalToGlobalMapping. We can implement the `MpiSolver` according to the stuff Bernd did. No problem here.

- build a graph that contains the mapping
- convert the Eigen::SparseMatrix from the assembler to a Epetra_CsrMatrix with this mapping
- solve using some solver provided by trilinos 
- return the solution


