// $Id$

#include "visualize/Point.h"
#include "visualize/VisualizeException.h"

NuTo::Point::Point(Eigen::Vector3d coordinates, int numData)
    : mCoordinates(coordinates)
{
    mData.resize(numData);
}

void NuTo::Point::SetData(int dataIndex, Eigen::VectorXd data)
{
    if (dataIndex >= this->mData.size())
        throw NuTo::VisualizeException(__PRETTY_FUNCTION__, "invalid data index.");
    mData[dataIndex] = data;
}

const Eigen::Vector3d& NuTo::Point::GetCoordinates() const
{
    return mCoordinates;
}

const Eigen::VectorXd& NuTo::Point::GetData(unsigned int rDataIndex) const
{
    return mData[rDataIndex];
}
