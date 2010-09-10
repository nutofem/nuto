// $ld: $ 
#ifndef VISUALIZECOMPONENTPRINCIPALENGINEERINGSTRESS_H_
#define VISUALIZECOMPONENTPRINCIPALENGINEERINGSTRESS_H_

#include <boost/serialization/vector.hpp>

#include "nuto/visualize/VisualizeComponentBase.h"

namespace NuTo
{
//! @author Joerg F. Unger
//! @date Apr 27, 2010
//! @brief ...this routine is normally not needed, since the stress tensor is export
//! unfortunately, the extraction of the principal stresses in paraview was not straigthforward
//! that is why this class is implemented
class VisualizeComponentPrincipalEngineeringStress : public VisualizeComponentBase
{
public:
    VisualizeComponentPrincipalEngineeringStress();

    inline NuTo::VisualizeBase::eVisualizeWhat GetComponentEnum()const
    {
    	return NuTo::VisualizeBase::PRINCIPAL_ENGINEERING_STRESS;
    }

    inline std::string GetComponentName()const
    {
    	return std::string("PrincipalEngineeringStress");
    }
};
}
#endif /* VISUALIZECOMPONENTPRINCIPALENGINEERINGSTRESS_H_ */
