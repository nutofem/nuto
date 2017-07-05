#include "visualize/Cell.h"
#include "visualize/VisualizeException.h"

using namespace NuTo::Visualize;

Cell::Cell(std::vector<int> pointIds, eCellTypes cellType, int numData)
    : mPointIds(pointIds)
    , mCellType(cellType)
{
    mData.resize(numData);
}

int Cell::GetNumPoints() const
{
    return mPointIds.size();
}

const std::vector<int>& Cell::GetPointIds() const
{
    return mPointIds;
}

void Cell::SetPointIds(std::vector<int> pointIds)
{
    mPointIds = pointIds;
}

NuTo::eCellTypes Cell::GetCellType() const
{
    return mCellType;
}

const Eigen::VectorXd& Cell::GetData(int dataIndex) const
{
    if (dataIndex >= static_cast<int>(this->mData.size()))
        throw VisualizeException(__PRETTY_FUNCTION__, "invalid data index.");
    return mData[dataIndex];
}

void Cell::SetData(int dataIndex, Eigen::VectorXd data)
{
    if (dataIndex >= static_cast<int>(this->mData.size()))
        throw NuTo::VisualizeException(__PRETTY_FUNCTION__, "invalid data index.");
    mData[dataIndex] = data;
}
