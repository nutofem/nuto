#include "visualize/Point.h"
#include "visualize/VisualizeException.h"

using namespace NuTo::Visualize;

Point::Point(Eigen::Vector3d coordinates, int numData)
    : mCoordinates(coordinates)
{
    mData.resize(numData);
}

void Point::SetData(int dataIndex, Eigen::VectorXd data)
{
    if (dataIndex >= this->mData.size())
        throw NuTo::VisualizeException(__PRETTY_FUNCTION__, "invalid data index.");
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
