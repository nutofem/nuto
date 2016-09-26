// $Id$
#pragma once

#include "nuto/mechanics/elements/ElementDataConstitutiveBase.h"
#include "nuto/mechanics/elements/ElementDataCrackBase.h"
#include "nuto/mechanics/elements/ElementDataIpBase.h"

namespace NuTo
{

//! @author Daniel Arnold
//! @date October 2010
//! @brief Standart class for ElementData IP, Constitutive and Cracks
class ElementDataConstitutiveIpCrack : public ElementDataConstitutiveBase, public ElementDataCrackBase, public ElementDataIpBase
{

public:
	//! @brief constructor
	ElementDataConstitutiveIpCrack(const ElementBase *rElement, const NuTo::IntegrationTypeBase* rIntegrationType, NuTo::IpData::eIpDataType rIpDataType);

	//! @brief constructor
	ElementDataConstitutiveIpCrack(const ElementBase *rElement, int rNumIp, NuTo::IpData::eIpDataType rIpDataType);

	//! @brief destructor
	virtual ~ElementDataConstitutiveIpCrack();

	//! @brief updates the data related to changes of the constitutive model (e.g. reallocation of static data, nonlocal weights etc.)
    //! @param rElement element
    virtual void InitializeUpdatedConstitutiveLaw(const ElementBase* rElement);

    //! @brief updates the data related to changes of the constitutive model (e.g. reallocation of static data, nonlocal weights etc.)
    //! @param rElement element
    //! @param rIp Ip
    virtual void InitializeUpdatedConstitutiveLaw(const ElementBase* rElement,int rIp);
		
    //! @brief returns the enum of element data type
    //! @return enum of ElementDataType
    const NuTo::ElementData::eElementDataType GetElementDataType()const;

protected:

};
}


