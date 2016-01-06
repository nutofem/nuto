// $Id$

#include "nuto/visualize/CellVertex.h"

// constructor
NuTo::CellVertex::CellVertex(const unsigned int *rPoints, const std::vector<VisualizeDataType>& rDataTypes):CellBase(rDataTypes)
{
    this->mPoints[0] = rPoints[0];
}

// number of cell points
unsigned int NuTo::CellVertex::GetNumPoints() const
{
    return 1;
}

// cell points
const unsigned int* NuTo::CellVertex::GetPoints() const
{
    return this->mPoints;
}

// Vtk cell type
unsigned int NuTo::CellVertex::GetVtkCellType() const
{
    return 1;
}
