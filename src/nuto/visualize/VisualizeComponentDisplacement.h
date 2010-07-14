// $ld: $ 
#ifndef VISUALIZECOMPONENTDISPLACEMENT_H_
#define VISUALIZECOMPONENTDISPLACEMENT_H_

#include <boost/serialization/vector.hpp>
#include <string>

#include "nuto/visualize/VisualizeComponentBase.h"
#include "nuto/visualize/VisualizeComponentBase.h"

namespace NuTo
{
//! @author Joerg F. Unger
//! @date Apr 27, 2010
//! @brief ...
class VisualizeComponentDisplacement : public VisualizeComponentBase
{
public:
	VisualizeComponentDisplacement();

    inline NuTo::VisualizeBase::eVisualizeWhat GetComponentEnum()const
    {
    	return NuTo::VisualizeBase::DISPLACEMENTS;
    }

    inline std::string GetComponentName()const
    {
    	return std::string("Displacements");
    }

};
}
#endif /* VISUALIZECOMPONENTDISPLACEMENT_H_ */
