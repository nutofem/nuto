/*
 * VisualizeComponentTotalInelasticEqStrain.h
 *
 *  Created on: Jun 30, 2010
 *      Author: unger3
 */

#ifndef VISUALIZECOMPONENTTOTALINELASTICEQSTRAIN_H_
#define VISUALIZECOMPONENTTOTALINELASTICEQSTRAIN_H_

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/visualize/VisualizeComponentBase.h"

namespace NuTo
{
//! @author V. Kindrachuk
//! @date Oct, 2015
//! @brief visualize total inelastic equivalent strain
class VisualizeComponentTotalInelasticEqStrain : public VisualizeComponentBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
	VisualizeComponentTotalInelasticEqStrain();

    inline NuTo::VisualizeBase::eVisualizeWhat GetComponentEnum()const
    {
    	return NuTo::VisualizeBase::TOTAL_INELASTIC_EQ_STRAIN;
    }

    inline std::string GetComponentName()const
    {
    	return std::string("TotalInelasticEqStrain");
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
BOOST_CLASS_EXPORT_KEY(NuTo::VisualizeComponentTotalInelasticEqStrain)
#endif // ENABLE_SERIALIZATION

#endif /* VISUALIZECOMPONENTTOTALINELASTICEQSTRAIN_H_ */
