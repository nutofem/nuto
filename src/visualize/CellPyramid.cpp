// $Id$

#include "visualize/CellPyramid.h"

// constructor
NuTo::CellPyramid::CellPyramid(const unsigned int *rPoints, const std::vector<VisualizeDataType>& rDataTypes):CellBase(rDataTypes)
{
    for (unsigned int count = 0; count < 5; count++)
    {
        this->mPoints[count] = rPoints[count];
    }
}

// number of cell points
unsigned int NuTo::CellPyramid::GetNumPoints() const
{
    return 5;
}

// cell points
const unsigned int* NuTo::CellPyramid::GetPoints() const
{
    return this->mPoints;
}

// Vtk cell type
unsigned int NuTo::CellPyramid::GetVtkCellType() const
{
    return 14;
}
