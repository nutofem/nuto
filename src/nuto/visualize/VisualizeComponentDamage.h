/*
 * VisualizeComponentDamage.h
 *
 *  Created on: Jun 30, 2010
 *      Author: unger3
 */

#ifndef VISUALIZECOMPONENTDAMAGE_H_
#define VISUALIZECOMPONENTDAMAGE_H_

#include "nuto/visualize/VisualizeComponentBase.h"

namespace NuTo
{
//! @author Joerg F. Unger
//! @date Apr 27, 2010
//! @brief visualize the nonlocal weights for integration point mIp in element mElement
class VisualizeComponentDamage : public VisualizeComponentBase
{
public:
	VisualizeComponentDamage() : VisualizeComponentBase::VisualizeComponentBase()
	{}

    inline NuTo::VisualizeBase::eVisualizeWhat GetComponentEnum()const
    {
    	return NuTo::VisualizeBase::DAMAGE;
    }

    inline std::string GetComponentName()const
    {
    	return std::string("DAMAGE");
    }

protected:
};
}

#endif /* VISUALIZECOMPONENTDAMAGE_H_ */
