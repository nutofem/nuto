// $Id$

#include "visualize/CellBase.h"
#include "visualize/VisualizeException.h"

NuTo::CellBase::CellBase(int num)
{
    mData.resize(num);
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

