/*
 * VisualizeComponentEngineeringPlasticStrain.h
 *
 *  Created on: Jun 30, 2010
 *      Author: unger3
 */

#ifndef VISUALIZECOMPONENTENGINEERINGPLASTICSTRAIN_H_
#define VISUALIZECOMPONENTENGINEERINGPLASTICSTRAIN_H_

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/visualize/VisualizeComponentBase.h"

namespace NuTo
{
//! @author Joerg F. Unger
//! @date Apr 27, 2010
//! @brief visualize the nonlocal weights for integration point mIp in element mElement
class VisualizeComponentEngineeringPlasticStrain : public VisualizeComponentBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
	VisualizeComponentEngineeringPlasticStrain();

    inline NuTo::VisualizeBase::eVisualizeWhat GetComponentEnum()const
    {
    	return NuTo::VisualizeBase::ENGINEERING_PLASTIC_STRAIN;
    }

    inline std::string GetComponentName()const
    {
    	return std::string("EngineeringPlasticStrain");
    }
#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

protected:
};
}

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::VisualizeComponentEngineeringPlasticStrain)
#endif // ENABLE_SERIALIZATION

#endif /* VISUALIZECOMPONENTENGINEERINGPLASTICSTRAIN_H_ */
