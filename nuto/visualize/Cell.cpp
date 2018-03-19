#include "nuto/visualize/Cell.h"
#include "nuto/base/Exception.h"

using namespace NuTo::Visualize;

Cell::Cell(std::vector<int> pointIds, eCellTypes cellType)
    : mPointIds(pointIds)
    , mCellType(cellType)
{
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
    if (dataIndex >= static_cast<int>(mData.size()))
        throw Exception(__PRETTY_FUNCTION__, "invalid data index.");
    return mData[dataIndex];
}

void Cell::SetData(int dataIndex, Eigen::VectorXd data)
{
    if (dataIndex >= static_cast<int>(mData.size()))
        mData.resize(dataIndex + 1);
    mData[dataIndex] = data;
}
