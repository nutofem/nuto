// $Id$
#ifndef VISUALIZECOMPONENTNONLOCALEQSTRAIN_H_
#define VISUALIZECOMPONENTNONLOCALEQSTRAIN_H_

#include <string>

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/visualize/VisualizeComponentBase.h"


namespace NuTo
{
//! @author Thomas Titscher
//! @date Dez 16, 2014
//! @brief ...
class VisualizeComponentNonlocalEqStrain : public VisualizeComponentBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
	VisualizeComponentNonlocalEqStrain();

    inline NuTo::VisualizeBase::eVisualizeWhat GetComponentEnum()const
    {
    	return NuTo::VisualizeBase::NONLOCAL_EQ_STRAIN;
    }

    inline std::string GetComponentName()const
    {
    	return std::string("NonlocalEqStrain");
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
BOOST_CLASS_EXPORT_KEY(NuTo::VisualizeComponentNonlocalEqStrain)
#endif // ENABLE_SERIALIZATION
#endif /* VISUALIZECOMPONENTNONLOCALEQSTRAIN_H_ */
