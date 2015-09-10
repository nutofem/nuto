// $Id$
#ifndef VISUALIZECOMPONENTLOCALEQSTRAIN_H_
#define VISUALIZECOMPONENTLOCALEQSTRAIN_H_

#include <string>

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/visualize/VisualizeComponentBase.h"


namespace NuTo
{
//! @author Thomas Titscher
//! @date Aug 26, 2015
//! @brief ...
class VisualizeComponentLocalEqStrain : public VisualizeComponentBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    VisualizeComponentLocalEqStrain();

    inline NuTo::VisualizeBase::eVisualizeWhat GetComponentEnum()const
    {
    	return NuTo::VisualizeBase::LOCAL_EQ_STRAIN;
    }

    inline std::string GetComponentName()const
    {
    	return std::string("LocalEqStrain");
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
BOOST_CLASS_EXPORT_KEY(NuTo::VisualizeComponentLocalEqStrain)
#endif // ENABLE_SERIALIZATION
#endif /* VISUALIZECOMPONENTLOCALEQSTRAIN_H_ */
