// $Id$
#ifndef ElementDataVariableConstitutiveIp_H
#define ElementDataVariableConstitutiveIp_H

#include "nuto/mechanics/elements/ElementDataVariableConstitutiveBase.h"
#include "nuto/mechanics/elements/ElementDataIpBase.h"
#include "nuto/mechanics/elements/IpDataEnum.h"

namespace NuTo
{
class IntegrationTypeBase;
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... class for elements with a single material per element static data for each integration point
class ElementDataVariableConstitutiveIp : public ElementDataVariableConstitutiveBase, public ElementDataIpBase
{
#ifdef ENABLE_SERIALIZATION
	friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION

public:
    //! @brief constructor
	ElementDataVariableConstitutiveIp(const ElementBase *rElement, const NuTo::IntegrationTypeBase* rIntegrationType, NuTo::IpData::eIpDataType rIpDataType);

	//! @brief constructor
	ElementDataVariableConstitutiveIp(const ElementBase *rElement, int rNumIp, NuTo::IpData::eIpDataType rIpDataType);

	virtual ~ElementDataVariableConstitutiveIp();

   //! @brief updates the data related to changes of the constitutive model (e.g. reallocation of static data, nonlocal weights etc.)
   //! @param rElement element
   virtual void InitializeUpdatedConstitutiveLaw(const ElementBase* rElement)
   {
 	   throw MechanicsException("[NuTo::ElementDataVariableConstitutiveIp::InitializeUpdatedConstitutiveLaw] need ip for constitutive law assignment.");
   }

    //! @brief updates the data related to changes of the constitutive model (e.g. reallocation of static data, nonlocal weights etc.)
    //! @param rElement element
   virtual void InitializeUpdatedConstitutiveLaw(const ElementBase* rElement, int rNumIp);

   //! @brief returns the enum of element data type
   //! @return enum of ElementDataType
   virtual const NuTo::ElementData::eElementDataType GetElementDataType()const;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

protected:
    //! @brief ...just for serialization
    ElementDataVariableConstitutiveIp(){}
};
}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ElementDataVariableConstitutiveIp)
namespace boost{
template<>
struct is_virtual_base_of<NuTo::ElementDataConstitutiveBase, NuTo::ElementDataVariableConstitutiveIp>: public mpl::true_ {};
template<>
struct is_virtual_base_of<NuTo::ElementDataIpBase, NuTo::ElementDataVariableConstitutiveIp>: public mpl::true_ {};
}
#endif // ENABLE_SERIALIZATION

#endif //ElementDataVariableConstitutiveIp_H
