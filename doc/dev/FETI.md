@page FETI Finite element tearing and interconnecting

# MPI application

The FETI method uses MPI for the communication between subdomains.

```cpp
#include <mpi.h>
// other includes

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    // your code

    MPI_Finalize();
}
```    

# Implementation

```cpp
NuTo::StructureFETI structure(dim);
structure.SetNumTimeDerivatives(0);

const int interpolationTypeId = structure.InterpolationTypeCreate(eShapeType::QUAD2D);
structure.InterpolationTypeAdd(interpolationTypeId, eDof::COORDINATES,     eTypeOrder::EQUIDISTANT1);
structure.InterpolationTypeAdd(interpolationTypeId, eDof::DISPLACEMENTS,   eTypeOrder::EQUIDISTANT1);

structure.ImportMeshJson(meshFile,interpolationTypeId);
```    

# Assembly of the connectivity matrix

![alt text][logo]

[logo]: doc/images/FETI_B_matrix_assembly.png "Logo Title Text 2"
