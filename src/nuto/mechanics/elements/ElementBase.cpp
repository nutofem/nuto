#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/sections/SectionBase.h"

//! @brief sets the section of an element
//! implemented with an exception for all elements, reimplementation required for those elements
//! which actually need a section
//! @param rSection pointer to section
//! @return pointer to constitutive law
void NuTo::ElementBase::SetSection(const SectionBase* rSection)
{
    throw MechanicsException("[NuTo::ElementBase::SetSection] Section for this type of elements not required.");
}

//! @brief returns a pointer to the section of an element
//! implemented with an exception for all elements, reimplementation required for those elements
//! which actually need a section
//! @return pointer to section
const NuTo::SectionBase* NuTo::ElementBase::GetSection()const
{
    throw MechanicsException("[NuTo::ElementBase::GetSection] Section for this type of elements not required.");
}

//! @brief sets the integration type of an element
//! implemented with an exception for all elements, reimplementation required for those elements
//! which actually need an integration type
//! @param rIntegrationType pointer to integration type
void NuTo::ElementBase::SetIntegrationType(const IntegrationTypeBase* rIntegrationType)
{
    throw MechanicsException("[NuTo::ElementBase::SetIntegrationType] Integration type for is type of elements not required.");
}

//! @brief returns a pointer to the integration type of an element
//! implemented with an exception for all elements, reimplementation required for those elements
//! which actually need an integration type
//! @return pointer to integration type
const NuTo::IntegrationTypeBase* NuTo::ElementBase::GetIntegrationType()const
{
    throw MechanicsException("[NuTo::ElementBase::GetIntegrationType] Integration type for is type of elements not required.");
}


void NuTo::ElementBase::InterpolateCoordinatesFrom1D(double rLocalCoordinates, double rGlobalCoordinates[3]) const
{
    throw MechanicsException("[NuTo::ElementBase::InterpolateCoordinatesFrom1D] 1D geometry interpolation routine not implemented.");
}

void NuTo::ElementBase::InterpolateCoordinatesFrom2D(double rLocalCoordinates[2], double rGlobalCoordinates[3]) const
{
    throw MechanicsException("[NuTo::ElementBase::InterpolateCoordinatesFrom2D] 2D geometry interpolation routine not implemented.");
}

void NuTo::ElementBase::InterpolateCoordinatesFrom3D(double rLocalCoordinates[3], double rGlobalCoordinates[3]) const
{
    throw MechanicsException("[NuTo::ElementBase::InterpolateCoordinatesFrom3D] 3D geometry interpolation routine not implemented.");
}

void NuTo::ElementBase::InterpolateDisplacementsFrom1D(double rLocalCoordinates, double rGlobalDisplacements[3]) const
{
    throw MechanicsException("[NuTo::ElementBase::InterpolateDisplacementsFrom1D] 1D displacement interpolation routine not implemented.");
}

void NuTo::ElementBase::InterpolateDisplacementsFrom2D(double rLocalCoordinates[2], double rGlobalDisplacements[3]) const
{
    throw MechanicsException("[NuTo::ElementBase::InterpolateDisplacementsFrom2D] 2D displacement interpolation routine not implemented.");
}

void NuTo::ElementBase::InterpolateDisplacementsFrom3D(double rLocalCoordinates[3], double rGlobalDisplacements[3]) const
{
    throw MechanicsException("[NuTo::ElementBase::InterpolateDisplacementsFrom3D] 3D displacement interpolation routine not implemented.");
}
