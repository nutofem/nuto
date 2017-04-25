// $Id$ 
#pragma once
#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include <string>
#include "visualize/VisualizeEnum.h"

namespace NuTo
{
namespace Visualize
{
class Component
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
    
    //! @brief just for serialization
    Component() = default;

#endif // ENABLE_SERIALIZATION
public:
    Component(NuTo::eVisualizeWhat visualizeComponent);

    NuTo::eVisualizeWhat GetComponentEnum() const;

    std::string GetComponentName() const;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

private:
    eVisualizeWhat mComponent;

};
} // namespace Visualize
} // namespace NuTo

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::Visualize::Component)
#endif // ENABLE_SERIALIZATION

