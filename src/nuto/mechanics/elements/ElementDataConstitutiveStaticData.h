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
#include "nuto/mechanics/elements/ElementDataConstitutiveBase.h"
#include "nuto/mechanics/elements/ElementDataStaticDataBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveStaticDataBase.h"

namespace NuTo
{
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... class for elements with a single material per element static data for each integration point
class ElementDataConstitutiveStaticData : public ElementDataConstitutiveBase, public ElementDataStaticDataBase
{

#ifdef ENABLE_SERIALIZATION
	friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION

public:
    //! @brief constructor
	ElementDataConstitutiveStaticData(const NuTo::IntegrationTypeBase* rIntegrationType);

	~ElementDataConstitutiveStaticData();

    //! @brief update the information related to a modification of the integration type, e.g. reallocation of the static data
    void UpdateForModifiedIntegrationType(const ElementBase* rElement);

    //! @brief update the information related to a modification of the constitutive law, e.g. reallocation of the static data
    void UpdateForModifiedConstitutiveLaw(const ElementBase* rElement);

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ElementDataConstitutiveBase)
           & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ElementDataStaticDataBase);
    }
#endif  // ENABLE_SERIALIZATION

protected:
};
}

#endif //ElementDataConstitutiveStaticData_H
