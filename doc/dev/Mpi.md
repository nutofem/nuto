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
We will call this mapping **LocalToGlobalMapping** and its calculation is the main issue. 

# Calculate **LocalToGlobalMapping**

__Assumption:__ The meshes we want to solve fit on one core. So on only one core, we can

- import/create a geometry mesh
- create as many "dof layers" as we want
- perform the numbering
- split the mesh into `nProc` smaller meshes (external mesh partitioning lib)
- via `boost::serialize` either
    - MPI communicate mesh to other processes
    - write meshes to disk and read them with other processes

Translated into a `main.cpp`, this could look like:

~~~cpp
MeshFem meshForThisRank;

if (rank == 0)
{
    MeshFem mesh = Mesh::ReadFromGmsh(pathToWholeGmshMesh);
    Mesh::AddDofInterpolation(&mesh, ...);
    Constraints dummyConstraintsForNumbering = Main::ApplyConstraints(mesh);
    DofNumbering::Build(&mesh, dummyConstraintsForNumbering);
    // Each node now has a _global_ dof number.
   
    std::vector<MeshFem> meshes = Mesh::Partition(mesh, nProcs);

    meshForThisRank = meshes[0];
    for (int iRank = 1; iRank < nProcs; ++iRank)
        boostMpiCommunicator.send(iRank, meshes[iRank]); // requires boost::serialize for MeshFem (!)
}
else
{
    boostMpiCommunicator.recieve(&meshForThisRank);
}

Constraints constraints = Main::ApplyConstraints(meshForThisRank);
auto localToGlobalMapping = DofNumbering::BuildMpi(&meshForThisRank, constraints);

auto cells = Main::CreateIntegrationCells(meshForThisRank);

auto domainMatrix = SimpleAssembler::BuildMatrix(cells, ...);
auto domainRhs = SimpleAssembler::BuildVector(cells, ...);
// apply constraint-matrix ...

MpiSolver solver(localToGlobalMapping);
auto domainSolution = solver.Solve(domainMatrix.JJ, domainRhs.J);
//...
~~~

## Some key functions

create the partitioned mesh

~~~cpp
std::vector<MeshFem> Mesh::Partition(MeshFem mesh, int nProcs)
{
    auto graph = Mesh::BuildGraph(mesh);
    std::vector<SomePartitionInfo> partitioningOutput = External::Partitioning(graph, nProcs);
    std::vector<Group<ElementFem>> partitionedElements = Mesh::PartitionPostprocess(partitioningOutput);
    return Mesh::SplitMesh(mesh, partitionedElements);
}
~~~


Building an mpi dof numbering

~~~cpp
std::vector<int> DofNumbering::BuildMpi(MeshFem* rMesh, Constraints constraints);
{
    // 1) extract global dof numbers (= current dof numbers at the nodes)
    // 2) perform DofNumbering::Build() to create a local numbering at the nodes
    // 3) compare both numberings to
    return localToGlobal;
}
~~~

mpi solver

~~~cpp
Eigen::VectorXd MpiSolver::Solve(Eigen::SparseMatrix<double> A, Eigen::VectorXd b)
{
    MpiHelper::AssembleLocalToGlobal(A, &mTpetraMatrix);
    MpiHelper::AssembleLocalToGlobal(b, &mTpetraVector);
    auto solution = mSomeTrilinosSolver.Solve(mTpetraMatrix, mTpetraVector);
    return MpiHelper::ToEigen(solution.GetContributionForThisRank());
}
~~~

