#include "nuto/visualize/Point.h"
#include "nuto/base/Exception.h"
#include "nuto/math/EigenCompanion.h"

using namespace NuTo::Visualize;

Point::Point(Eigen::VectorXd coordinates)
    : mCoordinates(EigenCompanion::To3D(coordinates))
{
}

void Point::SetData(int dataIndex, Eigen::VectorXd data)
{
    if (dataIndex >= static_cast<int>(mData.size()))
        mData.resize(dataIndex + 1);
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
