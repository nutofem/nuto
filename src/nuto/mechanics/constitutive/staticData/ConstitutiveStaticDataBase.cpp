// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#ifdef ENABLE_VISUALIZE
#include "nuto/visualize/VisualizeBase.h"
#include "nuto/visualize/VisualizeUnstructuredGrid.h"
#include <boost/ptr_container/ptr_list.hpp>
#endif // ENABLE_VISUALIZE

#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/constitutive/staticData/ConstitutiveStaticDataBase.h"


#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::ConstitutiveStaticDataBase::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataBase::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataBase::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataBase::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataBase::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataBase::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstitutiveStaticDataBase::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize constitutive static data base" << std::endl;
    std::cout << "finish serialize constitutive static data base" << std::endl;
#endif
}
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::ConstitutiveStaticDataBase)
#endif // ENABLE_SERIALIZATION

//!@ brief reinterpret as nonlocal damage2d static data
NuTo::ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain* NuTo::ConstitutiveStaticDataBase::AsNonlocalDamagePlasticity2DPlaneStrain()
{
    throw NuTo::MechanicsException(__PRETTY_FUNCTION__,"Static data is not of type NonlocalDamagePlasticity2DPlaneStrain.");
}

//!@ brief reinterpret as nonlocal damage2d static data
const NuTo::ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain* NuTo::ConstitutiveStaticDataBase::AsNonlocalDamagePlasticity2DPlaneStrain()const
{
    throw NuTo::MechanicsException(__PRETTY_FUNCTION__,"Static data is not of type NonlocalDamagePlasticity2DPlaneStrain.");
}

//!@ brief reinterpret as gradient damage 1d static data
NuTo::ConstitutiveStaticDataGradientDamagePlasticity1D* NuTo::ConstitutiveStaticDataBase::AsGradientDamagePlasticity1D()
{
    throw NuTo::MechanicsException(__PRETTY_FUNCTION__,"Static data is not of type GradientDamagePlasticity1D.");
}

//!@ brief reinterpret as nonlocal damage 1d static data
const NuTo::ConstitutiveStaticDataGradientDamagePlasticity1D* NuTo::ConstitutiveStaticDataBase::AsGradientDamagePlasticity1D()const
{
    throw NuTo::MechanicsException(__PRETTY_FUNCTION__,"Static data is not of type GradientDamagePlasticity1D.");
}

//!@ brief reinterpret as gradient damage 1d static data
NuTo::ConstitutiveStaticDataGradientDamage* NuTo::ConstitutiveStaticDataBase::AsGradientDamage()
{
    throw NuTo::MechanicsException(__PRETTY_FUNCTION__,"Static data is not of type GradientDamage.");
}

//!@ brief reinterpret as nonlocal damage 1d static data
const NuTo::ConstitutiveStaticDataGradientDamage* NuTo::ConstitutiveStaticDataBase::AsGradientDamage()const
{
    throw NuTo::MechanicsException(__PRETTY_FUNCTION__,"Static data is not of type GradientDamage.");
}

//!@ brief reinterpret as gradient damage 1d static data
NuTo::ConstitutiveStaticDataGradientDamage1DFatigue* NuTo::ConstitutiveStaticDataBase::AsGradientDamage1DFatigue()
{
    throw NuTo::MechanicsException(__PRETTY_FUNCTION__,"Static data is not of type GradientDamage1DFatigue.");
}

//!@ brief reinterpret as nonlocal damage 1d static data
const NuTo::ConstitutiveStaticDataGradientDamage1DFatigue* NuTo::ConstitutiveStaticDataBase::AsGradientDamage1DFatigue()const
{
    throw NuTo::MechanicsException(__PRETTY_FUNCTION__,"Static data is not of type GradientDamage1D.");
}

//!@ brief reinterpret as gradient damage 2d static data
NuTo::ConstitutiveStaticDataGradientDamage2DFatigue* NuTo::ConstitutiveStaticDataBase::AsGradientDamage2DFatigue()
{
    throw NuTo::MechanicsException(__PRETTY_FUNCTION__,"Static data is not of type GradientDamage2DFatigue.");
}

//!@ brief reinterpret as nonlocal damage 2d static data
const NuTo::ConstitutiveStaticDataGradientDamage2DFatigue* NuTo::ConstitutiveStaticDataBase::AsGradientDamage2DFatigue()const
{
    throw NuTo::MechanicsException(__PRETTY_FUNCTION__,"Static data is not of type GradientDamage2DFatigue.");
}

//!@ brief reinterpret as bond stress slip static data
NuTo::ConstitutiveStaticDataBondStressSlip* NuTo::ConstitutiveStaticDataBase::AsBondStressSlip()
{
    throw NuTo::MechanicsException(__PRETTY_FUNCTION__,"Static data is not of type BondStressSlip.");
}

//!@ brief reinterpret as bond stress slip static data
const NuTo::ConstitutiveStaticDataBondStressSlip* NuTo::ConstitutiveStaticDataBase::AsBondStressSlip()const
{
    throw NuTo::MechanicsException(__PRETTY_FUNCTION__,"Static data is not of type BondStressSlip.");
}

//!@ brief reinterpret as multiscale2d static data
NuTo::ConstitutiveStaticDataMultiscale2DPlaneStrain* NuTo::ConstitutiveStaticDataBase::AsMultiscale2DPlaneStrain()
{
    throw NuTo::MechanicsException(__PRETTY_FUNCTION__,"Static data is not of type NonlocalDamagePlasticity2DPlaneStrain.");
}

//!@ brief reinterpret as multiscale2d static data
const NuTo::ConstitutiveStaticDataMultiscale2DPlaneStrain* NuTo::ConstitutiveStaticDataBase::AsMultiscale2DPlaneStrain()const
{
    throw NuTo::MechanicsException(__PRETTY_FUNCTION__,"Static data is not of type NonlocalDamagePlasticity2DPlaneStrain.");
}

//!@ brief reinterpret as lattice concrete 2D static data
NuTo::ConstitutiveStaticDataLatticeConcrete2D* NuTo::ConstitutiveStaticDataBase::AsConstitutiveStaticDataLatticeConcrete2D()
{
    throw NuTo::MechanicsException(__PRETTY_FUNCTION__,"Static data is not of type ConstitutiveStaticDataLatticeConcrete2D.");
}

//!@ brief reinterpret as lattice concrete 2D static data
const NuTo::ConstitutiveStaticDataLatticeConcrete2D* NuTo::ConstitutiveStaticDataBase::AsConstitutiveStaticDataLatticeConcrete2D()const
{
    throw NuTo::MechanicsException(__PRETTY_FUNCTION__,"Static data is not of type ConstitutiveStaticDataLatticeConcrete2D.");
}

//!@ brief reinterpret as strain gradient damage1d static data
NuTo::ConstitutiveStaticDataStrainGradientDamagePlasticity1D* NuTo::ConstitutiveStaticDataBase::AsStrainGradientDamagePlasticity1D()
{
    throw NuTo::MechanicsException(__PRETTY_FUNCTION__,"Static data is not of type StrainGradientDamagePlasticity1D.");
}

//!@ brief reinterpret as strain gradient damage plasticity static data
const NuTo::ConstitutiveStaticDataStrainGradientDamagePlasticity1D* NuTo::ConstitutiveStaticDataBase::AsStrainGradientDamagePlasticity1D()const
{
    throw NuTo::MechanicsException(__PRETTY_FUNCTION__,"Static data is not of type GradientDamagePlasticity1D.");
}

//!@ brief reinterpret as damage viscoplasticity static data
NuTo::ConstitutiveStaticDataDamageViscoPlasticity3D* NuTo::ConstitutiveStaticDataBase::AsDamageViscoPlasticity3D()
{
    throw NuTo::MechanicsException(__PRETTY_FUNCTION__,"Static data is not of type ConstitutiveStaticDataDamageViscoPlasticity3D.");
}

//!@ brief reinterpret as damage viscoplasticity static data
const NuTo::ConstitutiveStaticDataDamageViscoPlasticity3D* NuTo::ConstitutiveStaticDataBase::AsDamageViscoPlasticity3D()const
{
    throw NuTo::MechanicsException(__PRETTY_FUNCTION__,"Static data is not of type ConstitutiveStaticDataDamageViscoPlasticity3D.");
}

//!@ brief reinterpret as damage viscoplasticity static data with damage
NuTo::ConstitutiveStaticDataDamageViscoPlasticity3DFatigue* NuTo::ConstitutiveStaticDataBase::AsDamageViscoPlasticity3DFatigue()
{
    throw NuTo::MechanicsException(__PRETTY_FUNCTION__,"Static data is not of type ConstitutiveStaticDataDamageViscoPlasticity3DFatigue.");
}

//!@ brief reinterpret as damage viscoplasticity static data with damage
const NuTo::ConstitutiveStaticDataDamageViscoPlasticity3DFatigue* NuTo::ConstitutiveStaticDataBase::AsDamageViscoPlasticity3DFatigue()const
{
    throw NuTo::MechanicsException(__PRETTY_FUNCTION__,"Static data is not of type ConstitutiveStaticDataDamageViscoPlasticity3DFatigue.");
}

//!@ brief reinterpret as moisture transport
NuTo::ConstitutiveStaticDataMoistureTransport* NuTo::ConstitutiveStaticDataBase::AsMoistureTransport()
{
    throw NuTo::MechanicsException(__PRETTY_FUNCTION__,"Static data is not of type ConstitutiveStaticDataMoistureTransport.");
}

//!@ brief reinterpret as moisture transport
const NuTo::ConstitutiveStaticDataMoistureTransport* NuTo::ConstitutiveStaticDataBase::AsMoistureTransport()const
{
    throw NuTo::MechanicsException(__PRETTY_FUNCTION__,"Static data is not of type ConstitutiveStaticDataMoistureTransport.");
}

//!@ brief reinterpret as multi physics
NuTo::ConstitutiveStaticDataMultipleConstitutiveLaws *NuTo::ConstitutiveStaticDataBase::AsMultipleConstitutiveLaws()
{
    throw NuTo::MechanicsException(__PRETTY_FUNCTION__,"Static data is not of type ConstitutiveStaticDataMultipleConstitutiveLaws.");
}

//!@ brief reinterpret as multi physics
const NuTo::ConstitutiveStaticDataMultipleConstitutiveLaws *NuTo::ConstitutiveStaticDataBase::AsMultipleConstitutiveLaws() const
{
    throw NuTo::MechanicsException(__PRETTY_FUNCTION__,"Static data is not of type ConstitutiveStaticDataMultipleConstitutiveLaws.");
}


void NuTo::ConstitutiveStaticDataBase::SetFineScaleModel(std::string rFileName, double rMacroLength, double rCenter[2], std::string rIPName)
{
    throw NuTo::MechanicsException(__PRETTY_FUNCTION__,"Static data has no fine scale model");
}


//! @brief sets the fine scale model parameters
void NuTo::ConstitutiveStaticDataBase::SetFineScaleParameter(const std::string& rName, double rParameter)
{
    throw NuTo::MechanicsException(__PRETTY_FUNCTION__,"Static data has no fine scale model");
}

//! @brief sets the fine scale model parameters
void NuTo::ConstitutiveStaticDataBase::SetFineScaleParameter(const std::string& rName, std::string rParameter)
{
    throw NuTo::MechanicsException(__PRETTY_FUNCTION__,"Static data has no fine scale model");
}

#ifdef ENABLE_VISUALIZE
//Visualize for all integration points the fine scale structure
void NuTo::ConstitutiveStaticDataBase::VisualizeIpMultiscale(VisualizeUnstructuredGrid& rVisualize,
		const boost::ptr_list<NuTo::VisualizeComponentBase>& rWhat, bool rVisualizeDamage)const
{
    //do nothing if not multiscale static data
}
#endif



