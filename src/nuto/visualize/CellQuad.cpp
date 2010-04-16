// $Id$

#include "nuto/visualize/CellQuad.h"

// constructor
NuTo::CellQuad::CellQuad(const unsigned int *rPoints, const std::vector<VisualizeDataType>& rDataTypes):CellBase(rDataTypes)
{
    for (unsigned int count = 0; count < 4; count++)
    {
        this->mPoints[count] = rPoints[count];
    }
}

// number of cell points
unsigned int NuTo::CellQuad::GetNumPoints() const
{
    return 4;
}

// cell points
const unsigned int* NuTo::CellQuad::GetPoints() const
{
    return this->mPoints;
}

// Vtk cell type
unsigned int NuTo::CellQuad::GetVtkCellType() const
{
    return 9;
}
