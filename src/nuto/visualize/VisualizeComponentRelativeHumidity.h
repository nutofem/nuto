// $Id$
#ifndef VISUALIZECOMPONENTRELATIVEHUMIDITY_H_
#define VISUALIZECOMPONENTRELATIVEHUMIDITY_H_

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
class VisualizeComponentRelativeHumidity : public VisualizeComponentBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    VisualizeComponentRelativeHumidity();

    inline NuTo::VisualizeBase::eVisualizeWhat GetComponentEnum()const
    {
        return NuTo::VisualizeBase::RELATIVE_HUMIDITY;
    }

    inline std::string GetComponentName()const
    {
        return std::string("RelativeHumidity");
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
BOOST_CLASS_EXPORT_KEY(NuTo::VisualizeComponentRelativeHumidity)
#endif // ENABLE_SERIALIZATION
#endif /* VISUALIZECOMPONENTRELATIVEHUMIDITY_H_ */
