#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/visualize/VisualizeComponentBase.h"

namespace NuTo
{
//! @author Philip Huschke
//! @date October 2015
//! @brief ...
class VisualizeComponentBondStress : public VisualizeComponentBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    VisualizeComponentBondStress();

    inline NuTo::VisualizeBase::eVisualizeWhat GetComponentEnum()const
    {
    	return NuTo::VisualizeBase::BOND_STRESS;
    }

    inline std::string GetComponentName()const
    {
    	return std::string("BondStress");
    }

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION
};
}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::VisualizeComponentBondStress)
#endif // ENABLE_SERIALIZATION
