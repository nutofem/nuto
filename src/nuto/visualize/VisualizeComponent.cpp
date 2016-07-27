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

#include "nuto/visualize/VisualizeComponent.h"
#include "nuto/visualize/VisualizeException.h"


NuTo::VisualizeComponent::VisualizeComponent(NuTo::VisualizeBase::eVisualizeWhat rVisualizeComponent) :
        mVisualizeComponent(rVisualizeComponent)
{
}



int NuTo::VisualizeComponent::GetElementId()const
{
    throw VisualizeException(__PRETTY_FUNCTION__, "Visualization component has no ElementId.");
}


int NuTo::VisualizeComponent::GetIp()const
{
    throw VisualizeException(__PRETTY_FUNCTION__, "Visualization component has no Integration point.");
}


NuTo::VisualizeBase::eVisualizeWhat NuTo::VisualizeComponent::GetComponentEnum() const
{
    return mVisualizeComponent;
}

std::string NuTo::VisualizeComponent::GetComponentName() const
{
    switch (mVisualizeComponent)
    {
    case VisualizeBase::ACCELERATION:
        return "Accelerations";
    case VisualizeBase::ANGULAR_ACCELERATION:
        return "AngularAccelerations";
    case VisualizeBase::ANGULAR_VELOCITY:
        return "AngularVelocities";
    case VisualizeBase::BOND_STRESS:
        return "BondStress";
    case VisualizeBase::CONSTITUTIVE:
        return "ConstitutiveModel";
    case VisualizeBase::CRACK:
        return "Crack";
    case VisualizeBase::DAMAGE:
        return "Damage";
    case VisualizeBase::DAMAGE_PHASE_FIELD:
        return "Damage_Phase_Field";
    case VisualizeBase::DISPLACEMENTS:
        return "Displacements";
    case VisualizeBase::ELEMENT:
        return "Element";
    case VisualizeBase::ENGINEERING_PLASTIC_STRAIN:
        return "EngineeringPlasticStrain";
    case VisualizeBase::ENGINEERING_STRAIN:
        return "EngineeringStrain";
    case VisualizeBase::ENGINEERING_STRESS:
        return "EngineeringStress";
    case VisualizeBase::HEAT_FLUX:
        return "HeatFlux";
    case VisualizeBase::LATTICE_STRAIN:
        return "LatticeStrain";
    case VisualizeBase::LATTICE_STRESS:
        return "LatticeStress";
    case VisualizeBase::LOCAL_EQ_STRAIN:
        return "LocalEqStrain";
    case VisualizeBase::NONLOCAL_EQ_STRAIN:
        return "NonlocalEqStrain";
    case VisualizeBase::PARTICLE_RADIUS:
        return "ParticleRadius";
    case VisualizeBase::PRINCIPAL_ENGINEERING_STRESS:
        return "PrincipalEngineeringStress";
    case VisualizeBase::RELATIVE_HUMIDITY:
        return "RelativeHumidity";
    case VisualizeBase::ROTATION:
        return "Rotations";
    case VisualizeBase::SECTION:
        return "Section";
    case VisualizeBase::SHRINKAGE_STRAIN:
        return "ShrinkageStrains";
    case VisualizeBase::SLIP:
        return "Slip";
    case VisualizeBase::TEMPERATURE:
        return "Temperature";
    case VisualizeBase::THERMAL_STRAIN:
        return "ThermalStrain";
    case VisualizeBase::TOTAL_INELASTIC_EQ_STRAIN:
        return "TotalInelasticEqStrain";
    case VisualizeBase::VELOCITY:
        return "Velocities";
    case VisualizeBase::WATER_VOLUME_FRACTION:
        return "WaterVolumeFraction";
    default:
        throw VisualizeException(__PRETTY_FUNCTION__, "Visualization component not implemented.");
    }
}



#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::VisualizeComponent::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::VisualizeComponent::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::VisualizeComponent::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::VisualizeComponent::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::VisualizeComponent::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::VisualizeComponent::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::VisualizeComponent::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize VisualizeComponent" << std::endl;
#endif
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize VisualizeComponent" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::VisualizeComponent)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::VisualizeComponent)
#endif // ENABLE_SERIALIZATION
