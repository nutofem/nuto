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
#include "nuto/visualize/VisualizeEnum.h"
#include "nuto/visualize/VisualizeException.h"


NuTo::VisualizeComponent::VisualizeComponent(NuTo::eVisualizeWhat rVisualizeComponent) :
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


NuTo::eVisualizeWhat NuTo::VisualizeComponent::GetComponentEnum() const
{
    return mVisualizeComponent;
}

std::string NuTo::VisualizeComponent::GetComponentName() const
{
    switch (mVisualizeComponent)
    {
    case eVisualizeWhat::ACCELERATION:
        return "Accelerations";
    case eVisualizeWhat::ANGULAR_ACCELERATION:
        return "AngularAccelerations";
    case eVisualizeWhat::ANGULAR_VELOCITY:
        return "AngularVelocities";
    case eVisualizeWhat::BOND_STRESS:
        return "BondStress";
    case eVisualizeWhat::CONSTITUTIVE:
        return "ConstitutiveModel";
    case eVisualizeWhat::DAMAGE:
        return "Damage";
    case eVisualizeWhat::CRACK_PHASE_FIELD:
        return "CrackPhaseField";
    case eVisualizeWhat::DISPLACEMENTS:
        return "Displacements";
    case eVisualizeWhat::ELEMENT:
        return "Element";
    case eVisualizeWhat::ENGINEERING_PLASTIC_STRAIN:
        return "EngineeringPlasticStrain";
    case eVisualizeWhat::ENGINEERING_STRAIN:
        return "EngineeringStrain";
    case eVisualizeWhat::ENGINEERING_STRESS:
        return "EngineeringStress";
    case eVisualizeWhat::HEAT_FLUX:
        return "HeatFlux";
    case eVisualizeWhat::LATTICE_STRAIN:
        return "LatticeStrain";
    case eVisualizeWhat::LATTICE_STRESS:
        return "LatticeStress";
    case eVisualizeWhat::LOCAL_EQ_STRAIN:
        return "LocalEqStrain";
    case eVisualizeWhat::NONLOCAL_EQ_STRAIN:
        return "NonlocalEqStrain";
    case eVisualizeWhat::PARTICLE_RADIUS:
        return "ParticleRadius";
    case eVisualizeWhat::PRINCIPAL_ENGINEERING_STRESS:
        return "PrincipalEngineeringStress";
    case eVisualizeWhat::RELATIVE_HUMIDITY:
        return "RelativeHumidity";
    case eVisualizeWhat::ROTATION:
        return "Rotations";
    case eVisualizeWhat::SECTION:
        return "Section";
    case eVisualizeWhat::SHRINKAGE_STRAIN:
        return "ShrinkageStrains";
    case eVisualizeWhat::SLIP:
        return "Slip";
    case eVisualizeWhat::TEMPERATURE:
        return "Temperature";
    case eVisualizeWhat::THERMAL_STRAIN:
        return "ThermalStrain";
    case eVisualizeWhat::TOTAL_INELASTIC_EQ_STRAIN:
        return "TotalInelasticEqStrain";
    case eVisualizeWhat::VELOCITY:
        return "Velocities";
    case eVisualizeWhat::WATER_VOLUME_FRACTION:
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
