// $ld: $ 
#ifndef VISUALIZECOMPONENTENGINEERINGSTRAIN_H_
#define VISUALIZECOMPONENTENGINEERINGSTRAIN_H_

#include <boost/serialization/vector.hpp>

#include "nuto/visualize/VisualizeComponentBase.h"

namespace NuTo
{
//! @author Joerg F. Unger
//! @date Apr 27, 2010
//! @brief ...
class VisualizeComponentEngineeringStrain : public VisualizeComponentBase
{
public:
	VisualizeComponentEngineeringStrain();

    inline NuTo::VisualizeBase::eVisualizeWhat GetComponentEnum()const
    {
    	return NuTo::VisualizeBase::ENGINEERING_STRAIN;
    }

    inline std::string GetComponentName()const
    {
    	return std::string("EngineeringStrain");
    }
};
}
#endif /* VISUALIZECOMPONENTENGINEERINGSTRAIN_H_ */
