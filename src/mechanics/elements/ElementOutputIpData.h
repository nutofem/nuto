#pragma once

#include "math/FullMatrix_Def.h"

#include "mechanics/elements/ElementOutputBase.h"

namespace NuTo
{
namespace IpData
{
    enum class eIpStaticDataType;
}// namespace IpData

//! @author Joerg F. Unger
//! @date Apr 29, 2010
//! @brief ...
class ElementOutputIpData : public ElementOutputBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    ElementOutputIpData() {}

    ElementOutputIpData(IpData::eIpStaticDataType rIpDataType)
    {
        mIpData[rIpDataType].resize(0,0);
    }


#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

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
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ElementOutputIpData)
#endif // ENABLE_SERIALIZATION
