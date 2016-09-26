// $Id$
#pragma once

#include "nuto/mechanics/elements/ElementDataBase.h"


namespace NuTo
{
//! @author Andrea Keszler
//! @date October 2014
//! @brief ... class for elements with a variable materials per element and no static data
class ElementDataVariableConstitutiveBase : public virtual ElementDataBase
{

#ifdef ENABLE_SERIALIZATION
	friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    ElementDataVariableConstitutiveBase();

    virtual ~ElementDataVariableConstitutiveBase();

    virtual void SetConstitutiveLaw(const ElementBase* rElement, NuTo::ConstitutiveBase* rConstitutiveLaw);

    //! @brief sets the constitutive law for a single integration point of the element
    //! @param rConstitutiveLaw constitutive law
    //! @param rIp integration point
    virtual void SetConstitutiveLaw(const ElementBase* rElement, int rIp, NuTo::ConstitutiveBase* rConstitutiveLaw);

    bool HasConstitutiveLawAssigned(int rIp) const override;

    //! @brief returns the constitutive law of an integration point
    //! @param rIp integration point
    //! @return constitutive law
    virtual  ConstitutiveBase* GetConstitutiveLaw(int rIp);

    //! @brief returns the constitutive law of an integration point
    //! @param rIp integration point
    //! @return constitutive law
    virtual  const ConstitutiveBase* GetConstitutiveLaw(int rIp)const;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

protected:
	std::vector<ConstitutiveBase*>  mVarConstitutiveLaw;
};
}//namespace NuTo

