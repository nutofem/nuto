/*
 * ElementDataConstitutiveBase.h
 *
 *  Created on: Apr 16, 2010
 *      Author: Joerg F. Unger
 */

#ifndef ELEMENTDATACONSTITUTIVEBASE_H_
#define ELEMENTDATACONSTITUTIVEBASE_H_

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
class ElementDataConstitutiveBase : public virtual ElementDataBase
{

#ifdef ENABLE_SERIALIZATION
	friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    ElementDataConstitutiveBase(const NuTo::IntegrationTypeBase* rIntegrationType);

    virtual void SetConstitutiveLaw(const ElementBase* rElement, NuTo::ConstitutiveBase* rConstitutiveLaw);

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

#endif /* ELEMENTDATACONSTITUTIVEBASE_H_ */
