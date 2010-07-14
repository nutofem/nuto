/*
 * VisualizeComponentSection.h
 *
 *  Created on: Jun 30, 2010
 *      Author: unger3
 */

#ifndef VISUALIZECOMPONENTSECTION_H_
#define VISUALIZECOMPONENTSECTION_H_

#include "nuto/visualize/VisualizeComponentBase.h"

namespace NuTo
{
//! @author Joerg F. Unger
//! @date Apr 27, 2010
//! @brief visualize the nonlocal weights for integration point mIp in element mElement
class VisualizeComponentSection : public VisualizeComponentBase
{
public:
	VisualizeComponentSection() : VisualizeComponentBase::VisualizeComponentBase()
	{}

    inline NuTo::VisualizeBase::eVisualizeWhat GetComponentEnum()const
    {
    	return NuTo::VisualizeBase::SECTION;
    }

    inline std::string GetComponentName()const
    {
    	return std::string("Section");
    }

protected:
};
}

#endif /* VISUALIZECOMPONENTSECTION_H_ */
