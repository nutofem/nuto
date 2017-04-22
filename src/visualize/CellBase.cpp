// $Id$

#include "visualize/CellBase.h"
#include "visualize/VisualizeException.h"

int CellTypeToVtk(NuTo::eCellTypes cellType)
{
    switch (cellType)
    {
    case NuTo::eCellTypes::HEXAHEDRON:
        return 12;
    case NuTo::eCellTypes::LINE:
        return 3;
    case NuTo::eCellTypes::PYRAMID:
        return 14;
    case NuTo::eCellTypes::QUAD:
        return 9;
    case NuTo::eCellTypes::TETRAEDER:
        return 10;
    case NuTo::eCellTypes::TRIANGLE:
        return 5;
    case NuTo::eCellTypes::VERTEX:
        return 1;
    }
}

NuTo::CellBase::CellBase(std::vector<int> pointIds, int numData, eCellTypes cellType)
    : mPointIds(pointIds)
    , mVtkCellType(CellTypeToVtk(cellType))
{
    mData.resize(numData);
}

int NuTo::CellBase::GetNumPoints() const
{
    return mPointIds.size();
}

//! @brief ... return point id's
//! @return ... array of point id's
const std::vector<int>& NuTo::CellBase::GetPointIds() const
{
    return mPointIds;
}

//! @brief ... setter for mPointIds
void NuTo::CellBase::SetPointIds(std::vector<int> pointIds)
{
    mPointIds = pointIds;
}

//! @brief ... returns the corresponding Vtk cell type
//! @return ... Vtk cell type
int NuTo::CellBase::GetVtkCellType() const
{
    return mVtkCellType;
}

const Eigen::VectorXd& NuTo::CellBase::GetData(int dataIndex) const
{
    if (dataIndex >= this->mData.size())
        throw VisualizeException(__PRETTY_FUNCTION__, "invalid data index.");
    return mData[dataIndex];
}

void NuTo::CellBase::SetData(int dataIndex, Eigen::VectorXd data)
{
    if (dataIndex >= this->mData.size())
        throw NuTo::VisualizeException(__PRETTY_FUNCTION__, "invalid data index.");
    mData[dataIndex] = data;
}
