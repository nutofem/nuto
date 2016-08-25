// $Id$ 
#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/visualize/VisualizeComponent.h"

namespace NuTo
{
class ElementBase;

//! @author Joerg F. Unger
//! @date Apr 27, 2010
//! @brief visualize the nonlocal weights for integration point mIp in element mElement
class VisualizeComponentNonlocalWeight : public VisualizeComponent
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    VisualizeComponentNonlocalWeight(int rElementId, int rIp);

    std::string GetComponentName(void) const override;

    int GetElementId(void) const override;

    int GetIp(void) const override;

    NuTo::eVisualizeWhat GetComponentEnum(void) const override;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

protected:
    //! @brief just for serialization
    VisualizeComponentNonlocalWeight(){};
    int mElementId;
    int	mIp;
};
}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::VisualizeComponentNonlocalWeight)
#endif // ENABLE_SERIALIZATION

