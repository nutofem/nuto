// $Id$
#ifndef VISUALIZECOMPONENTROTATION_H_
#define VISUALIZECOMPONENTROTATION_H_

#include <string>

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/visualize/VisualizeComponentBase.h"


namespace NuTo
{
//! @author Joerg F. Unger
//! @date Apr 27, 2010
//! @brief ...
class VisualizeComponentRotation : public VisualizeComponentBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
	VisualizeComponentRotation();

    inline NuTo::VisualizeBase::eVisualizeWhat GetComponentEnum()const
    {
    	return NuTo::VisualizeBase::ROTATION;
    }

    inline std::string GetComponentName()const
    {
    	return std::string("Rotations");
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
BOOST_CLASS_EXPORT_KEY(NuTo::VisualizeComponentRotation)
#endif // ENABLE_SERIALIZATION
#endif /* VISUALIZECOMPONENTROTATION_H_ */
