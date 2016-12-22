// $Id$

#include "visualize/CellHexahedron.h"

// constructor
NuTo::CellHexahedron::CellHexahedron(const unsigned int *rPoints, const std::vector<VisualizeDataType>& rDataTypes):CellBase(rDataTypes)
{
    for (unsigned int count = 0; count < 8; count++)
    {
        this->mPoints[count] = rPoints[count];
    }
}

// number of cell points
unsigned int NuTo::CellHexahedron::GetNumPoints() const
{
    return 8;
}

// cell points
const unsigned int* NuTo::CellHexahedron::GetPoints() const
{
    return this->mPoints;
}

// Vtk cell type
unsigned int NuTo::CellHexahedron::GetVtkCellType() const
{
    return 12;
}
