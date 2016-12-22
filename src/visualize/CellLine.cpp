// $Id$

#include "visualize/CellLine.h"

// constructor
NuTo::CellLine::CellLine(const unsigned int *rPoints, const std::vector<VisualizeDataType>& rDataTypes):CellBase(rDataTypes)
{
    this->mPoints[0] = rPoints[0];
    this->mPoints[1] = rPoints[1];
}

// number of cell points
unsigned int NuTo::CellLine::GetNumPoints() const
{
    return 2;
}

// cell points
const unsigned int* NuTo::CellLine::GetPoints() const
{
    return this->mPoints;
}

// Vtk cell type
unsigned int NuTo::CellLine::GetVtkCellType() const
{
    return 3;
}
