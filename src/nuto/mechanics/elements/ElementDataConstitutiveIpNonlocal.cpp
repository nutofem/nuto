// $Id: $
#include "nuto/mechanics/elements/ElementDataConstitutiveIpNonlocal.h"
#include "nuto/mechanics/constitutive/ConstitutiveBase.h"
#include <assert.h>


NuTo::ElementDataConstitutiveIpNonlocal::ElementDataConstitutiveIpNonlocal(const ElementWithDataBase *rElement,
		const NuTo::IntegrationTypeBase* rIntegrationType) :
   NuTo::ElementDataBase::ElementDataBase(), ElementDataConstitutiveBase(), ElementDataNonlocalBase() , ElementDataIpDataBase(rElement,rIntegrationType)
{
}

NuTo::ElementDataConstitutiveIpNonlocal::~ElementDataConstitutiveStaticDataNonlocal()
{
	//std::cout << "NuTo::ElementDataConstitutiveStaticDataNonlocal::~ElementDataConstitutiveStaticDataNonlocal()" << std::endl;
}
