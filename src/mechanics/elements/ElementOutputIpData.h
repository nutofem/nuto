#pragma once


#include "mechanics/elements/ElementOutputBase.h"

namespace NuTo
{
namespace IpData
{
enum class eIpStaticDataType;
} // namespace IpData

//! @author Joerg F. Unger
//! @date Apr 29, 2010
//! @brief ...
class ElementOutputIpData : public ElementOutputBase
{
public:
    ElementOutputIpData()
    {
    }

    ElementOutputIpData(IpData::eIpStaticDataType rIpDataType)
    {
        mIpData[rIpDataType].resize(0, 0);
    }


    ElementOutputIpData* Clone() const override
    {
        return new ElementOutputIpData(*this);
    }

    std::map<IpData::eIpStaticDataType, Eigen::MatrixXd>& GetIpDataMap()
    {
        return mIpData;
    }

    ElementOutputIpData& GetIpData() override
    {
        return *this;
    }

private:
    std::map<IpData::eIpStaticDataType, Eigen::MatrixXd> mIpData;
};
}
