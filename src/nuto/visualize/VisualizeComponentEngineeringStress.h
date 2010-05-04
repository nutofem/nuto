// $ld: $ 
#ifndef VISUALIZECOMPONENTENGINEERINGSTRESS_H_
#define VISUALIZECOMPONENTENGINEERINGSTRESS_H_

#include <boost/serialization/vector.hpp>

#include "nuto/visualize/VisualizeComponentBase.h"

namespace NuTo
{
//! @author Joerg F. Unger
//! @date Apr 27, 2010
//! @brief ...
class VisualizeComponentEngineeringStress : public VisualizeComponentBase
{
public:
	VisualizeComponentEngineeringStress();

    inline NuTo::VisualizeBase::eVisualizeWhat GetComponentEnum()const
    {
    	return NuTo::VisualizeBase::ENGINEERING_STRESS;
    }

    inline std::string GetComponentName()const
    {
    	return std::string("ENGINEERING_STRESS");
    }
};
}
#endif /* VISUALIZECOMPONENTENGINEERINGSTRESS_H_ */
