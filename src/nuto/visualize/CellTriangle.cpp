// $Id$

#include "nuto/visualize/CellTriangle.h"

// constructor
NuTo::CellTriangle::CellTriangle(const unsigned int *rPoints, const std::vector<VisualizeDataType>& rDataTypes):CellBase(rDataTypes)
{
    for (unsigned int count = 0; count < 3; count++)
    {
        this->mPoints[count] = rPoints[count];
    }
}

// number of cell points
unsigned int NuTo::CellTriangle::GetNumPoints() const
{
    return 3;
}

// cell points
const unsigned int* NuTo::CellTriangle::GetPoints() const
{
    return this->mPoints;
}

// Vtk cell type
unsigned int NuTo::CellTriangle::GetVtkCellType() const
{
    return 5;
}
