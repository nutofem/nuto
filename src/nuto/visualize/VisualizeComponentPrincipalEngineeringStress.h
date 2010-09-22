// $ld: $ 
#ifndef VISUALIZECOMPONENTPRINCIPALENGINEERINGSTRESS_H_
#define VISUALIZECOMPONENTPRINCIPALENGINEERINGSTRESS_H_

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

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
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
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

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION
};
}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::VisualizeComponentPrincipalEngineeringStress)
#endif // ENABLE_SERIALIZATION
#endif /* VISUALIZECOMPONENTPRINCIPALENGINEERINGSTRESS_H_ */
