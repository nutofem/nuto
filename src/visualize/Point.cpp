#include "visualize/Point.h"
#include "base/Exception.h"

using namespace NuTo::Visualize;

Point::Point(Eigen::Vector3d coordinates, int numData)
    : mCoordinates(coordinates)
{
    mData.resize(numData);
}

void Point::SetData(int dataIndex, Eigen::VectorXd data)
{
    if (dataIndex >= static_cast<int>(this->mData.size()))
        throw NuTo::Exception(__PRETTY_FUNCTION__, "invalid data index.");
    mData[dataIndex] = data;
}

const Eigen::Vector3d& Point::GetCoordinates() const
{
    return mCoordinates;
}

const Eigen::VectorXd& Point::GetData(int dataIndex) const
{
    return mData[dataIndex];
}
