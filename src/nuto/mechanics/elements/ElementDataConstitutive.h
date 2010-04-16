#ifndef ELEMENT_DATA_CONSTITUTIVE_H
#define ELEMENT_DATA_CONSTITUTIVE_H

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/elements/ElementDataBase.h"

namespace NuTo
{
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... class for elements with a single material per element and no static data
class ElementDataConstitutive : public ElementDataBase
{

#ifdef ENABLE_SERIALIZATION
	friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    ElementDataConstitutive(const NuTo::IntegrationTypeBase* rIntegrationType);
    
    void SetConstitutiveLaw(const ElementBase* rElement, NuTo::ConstitutiveBase* rConstitutiveLaw);

    void SetConstitutiveLaw(const ElementBase* rElement, int rIp, NuTo::ConstitutiveBase* rConstitutiveLaw);

    ConstitutiveStaticDataBase* GetStaticData(int rIp);

    const ConstitutiveStaticDataBase* GetStaticData(int rIp)const;

    ConstitutiveBase* GetConstitutiveLaw(int rIp);

    const ConstitutiveBase* GetConstitutiveLaw(int rIp)const;

    //! @brief update the information related to a modification of the integration type, e.g. reallocation of the static data
    void UpdateForModifiedIntegrationType(const ElementBase* rElement);

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {    
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ElementDataBase)
           & BOOST_SERIALIZATION_NVP(mConstitutiveLaw);
    }
#endif  // ENABLE_SERIALIZATION
    
protected:
    ConstitutiveBase* mConstitutiveLaw;
};
}//namespace NuTo
#endif //ELEMENT_DATA_CONSTITUTIVE_H
