#include "visualize/CellBase.h"
#include "visualize/VisualizeException.h"

NuTo::CellBase::CellBase(std::vector<int> pointIds, eCellTypes cellType, int numData)
    : mPointIds(pointIds)
    , mCellType(cellType)
{
    mData.resize(numData);
}

int NuTo::CellBase::GetNumPoints() const
{
    return mPointIds.size();
}

const std::vector<int>& NuTo::CellBase::GetPointIds() const
{
    return mPointIds;
}

void NuTo::CellBase::SetPointIds(std::vector<int> pointIds)
{
    mPointIds = pointIds;
}

NuTo::eCellTypes NuTo::CellBase::GetCellType() const
{
    return mCellType;
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
