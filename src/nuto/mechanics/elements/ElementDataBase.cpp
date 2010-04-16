#include <assert.h>
#include "nuto/mechanics/elements/ElementDataBase.h"

//! @brief sets the integration type of an element
//! implemented with an exception for all elements, reimplementation required for those elements
//! which actually need an integration type
//! @param rIntegrationType pointer to integration type
void NuTo::ElementDataBase::SetIntegrationType(const ElementBase* rElement, const NuTo::IntegrationTypeBase* rIntegrationType)
{
	mIntegrationType = rIntegrationType;
	UpdateForModifiedIntegrationType(rElement);
}

//! @brief returns a pointer to the integration type of an element
//! implemented with an exception for all elements, reimplementation required for those elements
//! which actually need an integration type
//! @return pointer to integration type
const NuTo::IntegrationTypeBase* NuTo::ElementDataBase::GetIntegrationType()const
{
	assert(mIntegrationType!=0);
	return mIntegrationType;
}
