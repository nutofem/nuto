// $Id: $
#ifndef ElementDataConstitutiveIpData_H
#define ElementDataConstitutiveIpData_H

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
#include "nuto/mechanics/elements/ElementDataIpBase.h"
#include "nuto/mechanics/elements/IpDataEnum.h"

namespace NuTo
{
class IntegrationTypeBase;
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... class for elements with a single material per element static data for each integration point
class ElementDataConstitutiveIp : public ElementDataConstitutiveBase, public ElementDataIpBase
{
#ifdef ENABLE_SERIALIZATION
	friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION

public:
    //! @brief constructor
	ElementDataConstitutiveIp(const ElementBase *rElement, const NuTo::IntegrationTypeBase* rIntegrationType, NuTo::IpData::eIpDataType rIpDataType);

	virtual ~ElementDataConstitutiveIp();

    //! @brief updates the data related to changes of the constitutive model (e.g. reallocation of static data, nonlocal weights etc.)
    //! @param rElement element
   virtual void InitializeUpdatedConstitutiveLaw(const ElementBase* rElement);

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ElementDataConstitutiveBase)
           & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ElementDataIpBase);
    }
#endif  // ENABLE_SERIALIZATION

protected:
};
}

#endif //ElementDataConstitutiveIpData_H
