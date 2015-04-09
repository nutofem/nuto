// $Id$
#ifndef VISUALIZECOMPONENTWATERVOLUMEFRACTION_H_
#define VISUALIZECOMPONENTWATERVOLUMEFRACTION_H_

#include <string>

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/visualize/VisualizeComponentBase.h"


namespace NuTo
{
//! @author Volker Hirthammer
//! @date Apr 08, 2015
//! @brief ...
class VisualizeComponentWaterVolumeFraction : public VisualizeComponentBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    VisualizeComponentWaterVolumeFraction();

    inline NuTo::VisualizeBase::eVisualizeWhat GetComponentEnum()const
    {
        return NuTo::VisualizeBase::WATER_VOLUME_FRACTION;
    }

    inline std::string GetComponentName()const
    {
        return std::string("WaterVolumeFraction");
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
BOOST_CLASS_EXPORT_KEY(NuTo::VisualizeComponentWaterVolumeFraction)
#endif // ENABLE_SERIALIZATION
#endif /* VISUALIZECOMPONENTWATERVOLUMEFRACTION_H_ */
