// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/base/Logger.h"
#include "nuto/mechanics/constitutive/ConstitutiveInputBase.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient1D.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient2D.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient3D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress1D.h"
#include "nuto/mechanics/constitutive/mechanics/NonlocalEqPlasticStrain.h"
#include "nuto/mechanics/constitutive/thermal/Temperature.h"
#include "nuto/mechanics/constitutive/thermal/TemperatureGradient3D.h"

#include "nuto/mechanics/MechanicsException.h"

// constructor
NuTo::ConstitutiveInputBase::ConstitutiveInputBase()
{
}

const NuTo::DeformationGradient1D& NuTo::ConstitutiveInputBase::GetDeformationGradient1D()const
{
	throw MechanicsException("[NuTo::ConstitutiveInputBase::GetDeformationGradient1D] not implemented for this input object.");
}

const NuTo::DeformationGradient2D& NuTo::ConstitutiveInputBase::GetDeformationGradient2D()const
{
	throw MechanicsException("[NuTo::ConstitutiveInputBase::GetDeformationGradient2D] not implemented for this input object.");
}

const NuTo::DeformationGradient3D& NuTo::ConstitutiveInputBase::GetDeformationGradient3D()const
{
	throw MechanicsException("[NuTo::ConstitutiveInputBase::GetDeformationGradient3D] not implemented for this input object.");
}

double NuTo::ConstitutiveInputBase::GetTemperature()const
{
	throw MechanicsException("[NuTo::ConstitutiveInputBase::GetTemperature] not implemented for this input object.");
}

const NuTo::TemperatureGradient3D& NuTo::ConstitutiveInputBase::GetTemperatureGradient3D()const
{
	throw MechanicsException("[NuTo::ConstitutiveInputBase::GetTemperatureGradient3D] not implemented for this input object.");
}

const NuTo::NonlocalEqPlasticStrain& NuTo::ConstitutiveInputBase::GetNonlocalEqPlasticStrain()const
{
	throw MechanicsException("[NuTo::ConstitutiveInputBase::GetNonlocalEqPlasticStrain] not implemented for this input object.");
}

const NuTo::NonlocalEqStrain& NuTo::ConstitutiveInputBase::GetNonlocalEqStrain()const
{
    throw MechanicsException("[NuTo::ConstitutiveInputBase::GetNonlocalEqStrain] not implemented for this input object.");
}

const NuTo::EngineeringStrain1D& NuTo::ConstitutiveInputBase::GetEngineeringStrain1D()const
{
	throw MechanicsException("[NuTo::ConstitutiveInputBase::GetEngineeringStrain1D] not implemented for this input object.");
}

const NuTo::EngineeringStrain2D& NuTo::ConstitutiveInputBase::GetEngineeringStrain2D()const
{
	throw MechanicsException("[NuTo::ConstitutiveInputBase::GetEngineeringStrain2D] not implemented for this input object.");
}

const NuTo::EngineeringStrain3D& NuTo::ConstitutiveInputBase::GetEngineeringStrain3D()const
{
	throw MechanicsException("[NuTo::ConstitutiveInputBase::GetEngineeringStrain3D] not implemented for this input object.");
}

const NuTo::EngineeringStress1D& NuTo::ConstitutiveInputBase::GetEngineeringStress1D()const
{
	throw MechanicsException("[NuTo::ConstitutiveInputBase::GetEngineeringStress1D] not implemented for this input object.");
}

const NuTo::RelativeHumidity& NuTo::ConstitutiveInputBase::GetRelativeHumidity()const
{
    throw MechanicsException("[NuTo::ConstitutiveInputBase::GetRelativeHumidity] not implemented for this input object.");
}

const NuTo::RelativeHumidityGradient2D &NuTo::ConstitutiveInputBase::GetRelativeHumidityGradient2D() const
{
    throw MechanicsException("[NuTo::ConstitutiveInputBase::GetRelativeHumidityGradient2D] not implemented for this input object.");
}

const NuTo::RelativeHumidityGradient3D &NuTo::ConstitutiveInputBase::GetRelativeHumidityGradient3D() const
{
    throw MechanicsException("[NuTo::ConstitutiveInputBase::GetRelativeHumidityGradient3D] not implemented for this input object.");
}

const NuTo::WaterVolumeFraction& NuTo::ConstitutiveInputBase::GetWaterVolumeFraction()const
{
    throw MechanicsException("[NuTo::ConstitutiveInputBase::GetWaterVolumeFraction] not implemented for this input object.");
}

const NuTo::WaterVolumeFractionGradient2D &NuTo::ConstitutiveInputBase::GetWaterVolumeFractionGradient2D() const
{
    throw MechanicsException("[NuTo::ConstitutiveInputBase::GetWaterVolumeFractionGradient2D] not implemented for this input object.");
}

const NuTo::WaterVolumeFractionGradient3D &NuTo::ConstitutiveInputBase::GetWaterVolumeFractionGradient3D() const
{
    throw MechanicsException("[NuTo::ConstitutiveInputBase::GetWaterVolumeFractionGradient3D] not implemented for this input object.");
}

const NuTo::InterfaceSlip &NuTo::ConstitutiveInputBase::GetInterfaceSlip() const
{
    throw MechanicsException("[NuTo::ConstitutiveInputBase::GetInterfaceSlip] not implemented for this input object.");
}


#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::ConstitutiveInputBase::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveInputBase::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveInputBase::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveInputBase::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveInputBase::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveInputBase::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstitutiveInputBase::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize constitutiveInputBase" << std::endl;
#endif
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstitutiveInputBase" << std::endl;
#endif
}
#endif // ENABLE_SERIALIZATION
