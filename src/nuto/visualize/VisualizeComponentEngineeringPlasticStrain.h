/*
 * VisualizeComponentEngineeringPlasticStrain.h
 *
 *  Created on: Jun 30, 2010
 *      Author: unger3
 */

#ifndef VISUALIZECOMPONENTENGINEERINGPLASTICSTRAIN_H_
#define VISUALIZECOMPONENTENGINEERINGPLASTICSTRAIN_H_
#include "nuto/visualize/VisualizeComponentBase.h"

namespace NuTo
{
//! @author Joerg F. Unger
//! @date Apr 27, 2010
//! @brief visualize the nonlocal weights for integration point mIp in element mElement
class VisualizeComponentEngineeringPlasticStrain : public VisualizeComponentBase
{
public:
	VisualizeComponentEngineeringPlasticStrain() : VisualizeComponentBase::VisualizeComponentBase()
	{}

    inline NuTo::VisualizeBase::eVisualizeWhat GetComponentEnum()const
    {
    	return NuTo::VisualizeBase::ENGINEERING_PLASTIC_STRAIN;
    }

    inline std::string GetComponentName()const
    {
    	return std::string("ENGINEERING_PLASTIC_STRAIN");
    }

protected:
};
}



#endif /* VISUALIZECOMPONENTENGINEERINGPLASTICSTRAIN_H_ */
