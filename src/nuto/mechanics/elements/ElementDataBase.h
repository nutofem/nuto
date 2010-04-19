#ifndef ELEMENT_DATA_BASE_H
#define ELEMENT_DATA_BASE_H

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/MechanicsException.h"

namespace NuTo
{
class ConstitutiveStaticDataBase;
class ConstitutiveBase;
class IntegrationTypeBase;
class ElementBase;

//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... standard abstract class for all possible element data classes
class ElementDataBase
{
#ifdef ENABLE_SERIALIZATION
	friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION

public:
    enum eElementDataType
    {
        CONSTITUTIVELAWELEMENT_NOSTATICDATA=0,   //!< constitutive law stored for at element level, no static data at ip
        CONSTITUTIVELAWELEMENT_STATICDATA,     //!< constitutive law stored for at element level, static data at ip
        CONSTITUTIVELAWIP_NOSTATICDATA,        //!< constitutive law stored for at integration point level, no static data at ip
        CONSTITUTIVELAWIP_STATICDATA,          //!< constitutive law stored for at integration point level, static data at ip
    };

    //! @brief constructor
    ElementDataBase(const NuTo::IntegrationTypeBase* rIntegrationType);

    //! @brief constructor
    virtual ~ElementDataBase();
    
    //! @brief sets the constitutive law for all integration points of the element
    //! @param rConstitutiveLaw constitutive law
    virtual void SetConstitutiveLaw(const ElementBase* rElement, NuTo::ConstitutiveBase* rConstitutiveLaw);

    //! @brief sets the constitutive law for a single integration point of the element
    //! @param rConstitutiveLaw constitutive law
    //! @param rIp integration point
   virtual void SetConstitutiveLaw(const ElementBase* rElement, int rIp, NuTo::ConstitutiveBase* rConstitutiveLaw);

    //! @brief returns the static data of an integration point
    //! @param rIp integration point
    //! @return static data
    virtual ConstitutiveStaticDataBase* GetStaticData(int rIp);

    //! @brief returns the static data of an integration point
    //! @param rIp integration point
    //! @return static data
   virtual const ConstitutiveStaticDataBase* GetStaticData(int rIp)const;

    //! @brief returns the constitutive law of an integration point
    //! @param rIp integration point
    //! @return constitutive law
    virtual  ConstitutiveBase* GetConstitutiveLaw(int rIp);

    //! @brief returns the constitutive law of an integration point
    //! @param rIp integration point
    //! @return constitutive law
    virtual  const ConstitutiveBase* GetConstitutiveLaw(int rIp)const;

    //! @brief sets the integration type of an element
    //! implemented with an exception for all elements, reimplementation required for those elements
    //! which actually need an integration type
    //! @param rElement pointer to element
    //! @param rIntegrationType pointer to integration type
    void SetIntegrationType(const ElementBase* rElement, const NuTo::IntegrationTypeBase* rIntegrationType);

    //! @brief returns a pointer to the integration type of an element
    //! implemented with an exception for all elements, reimplementation required for those elements
    //! which actually need an integration type
    //! @return pointer to integration type
    const IntegrationTypeBase* GetIntegrationType()const;

    //! @brief update the information related to a modification of the integration type, e.g. reallocation of the static data
    //! @param rElement pointer to element
    //! @param rIp integration point
    virtual void UpdateForModifiedIntegrationType(const ElementBase* rElement)=0;

    //! @brief update the information related to a modification of the constitutive law, e.g. reallocation of the static data, ip data
    //! @param rElement pointer to element
    //! @param rIp integration point
    virtual void UpdateForModifiedConstitutiveLaw(const ElementBase* rElement)=0;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {    
        ar & BOOST_SERIALIZATION_NVP(mIntegrationType);
    }
#endif  // ENABLE_SERIALIZATION
    
protected:
    //ElementDataBase(){};
    const IntegrationTypeBase *mIntegrationType;
};
}//namespace NuTo
#endif //ELEMENT_DATA_BASE_H
