/*
 * VisualizeComponentConstitutive.h
 *
 *  Created on: Jun 30, 2010
 *      Author: unger3
 */

#ifndef VISUALIZECOMPONENTCONSTITUTIVE_H_
#define VISUALIZECOMPONENTCONSTITUTIVE_H_

#include "nuto/visualize/VisualizeComponentBase.h"

namespace NuTo
{
//! @author Joerg F. Unger
//! @date Apr 27, 2010
//! @brief visualize the nonlocal weights for integration point mIp in element mElement
class VisualizeComponentConstitutive : public VisualizeComponentBase
{
public:
	VisualizeComponentConstitutive() : VisualizeComponentBase::VisualizeComponentBase()
	{}

    inline NuTo::VisualizeBase::eVisualizeWhat GetComponentEnum()const
    {
    	return NuTo::VisualizeBase::CONSTITUTIVE;
    }

    inline std::string GetComponentName()const
    {
    	return std::string("ConstitutiveModel");
    }

protected:
};
}

#endif /* VISUALIZECOMPONENTCONSTITUTIVE_H_ */
