// $Id$

#include "nuto/visualize/CellTetra.h"

// constructor
NuTo::CellTetra::CellTetra(const unsigned int *rPoints, const std::vector<VisualizeDataType>& rDataTypes):CellBase(rDataTypes)
{
    for (unsigned int count = 0; count < 4; count++)
    {
        this->mPoints[count] = rPoints[count];
    }
}

// number of cell points
unsigned int NuTo::CellTetra::GetNumPoints() const
{
    return 4;
}

// cell points
const unsigned int* NuTo::CellTetra::GetPoints() const
{
    return this->mPoints;
}

// Vtk cell type
unsigned int NuTo::CellTetra::GetVtkCellType() const
{
    return 10;
}
