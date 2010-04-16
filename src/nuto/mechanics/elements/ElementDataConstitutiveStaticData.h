// $Id: $
#ifndef ElementDataConstitutiveStaticData_H
#define ElementDataConstitutiveStaticData_H

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include <boost/ptr_container/ptr_vector.hpp>
#include "nuto/mechanics/elements/ElementDataBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveStaticDataBase.h"

namespace NuTo
{
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... class for elements with a single material per element static data for each integration point
class ElementDataConstitutiveStaticData : public ElementDataBase
{

#ifdef ENABLE_SERIALIZATION
	friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION

public:
    //! @brief constructor
	ElementDataConstitutiveStaticData(const NuTo::IntegrationTypeBase* rIntegrationType);

	~ElementDataConstitutiveStaticData();

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
           & BOOST_SERIALIZATION_NVP(mConstitutiveLaw)
           & BOOST_SERIALIZATION_NVP(mStaticData);
    }
#endif  // ENABLE_SERIALIZATION

protected:
    ConstitutiveBase* mConstitutiveLaw;
    // if the class is extended (e.g. for nonlocal data, plot data, stresses or strains at the integration point
    // please use a vector<IP_data >, where IP_data is a struct comprising all the required data
    // the size of mStaticData corresponds to the number of integration points, so for each integration point a static data object ist stored
    // no ptr_vector is used since it does not (in general) support null pointer
    std::vector<ConstitutiveStaticDataBase*> mStaticData;
};
}

#endif //ElementDataConstitutiveStaticData_H
