// $ld: $ 
#ifndef ELEMENTDATACONSTITUTIVEIPNONLOCAL_H_
#define ELEMENTDATACONSTITUTIVEIPNONLOCAL_H_

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/elements/ElementDataConstitutiveBase.h"
#include "nuto/mechanics/elements/ElementDataNonlocalBase.h"
#include "nuto/mechanics/elements/ElementDataIpBase.h"
#include "nuto/mechanics/elements/IpDataBase.h"

namespace NuTo
{

//! @author Joerg F. Unger
//! @date Apr 23, 2010
//! @brief ...
class ElementDataConstitutiveIpNonlocal : public ElementDataConstitutiveBase, public ElementDataNonlocalBase, public ElementDataIpBase
{

#ifdef ENABLE_SERIALIZATION
	friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION

public:
	//! @brief constructor
	ElementDataConstitutiveIpDataNonlocal(const ElementWithDataBase *rElement, IpDataBase::eIpDataType rIpDataType, const NuTo::IntegrationTypeBase* rIntegrationType);

	virtual ~ElementDataConstitutiveIpDataStaticDataNonlocalWeights();

#ifdef ENABLE_SERIALIZATION
	//! @brief serializes the class
	//! @param ar         archive
	//! @param version    version
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
		ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ElementDataConstitutiveBase)
		   & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ElementDataIpDataBase)
		   & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ElementDataNonlocalBase);
	}
#endif  // ENABLE_SERIALIZATION

protected:

};
}
#endif /* ELEMENTDATACONSTITUTIVEIPNONLOCAL_H_ */
