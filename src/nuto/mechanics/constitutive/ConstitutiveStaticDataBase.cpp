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
#include "nuto/mechanics/constitutive/ConstitutiveStaticDataBase.h"

NuTo::ConstitutiveStaticDataBase& NuTo::ConstitutiveStaticDataBase::operator= (NuTo::ConstitutiveStaticDataBase const& rOther)
{
    return (*this);
}


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
    throw NuTo::MechanicsException("[NuTo::ConstitutiveStaticDataBase::AsNonlocalDamagePlasticity2DPlaneStrain] Static data is not of type NonlocalDamagePlasticity2DPlaneStrain.");
}

//!@ brief reinterpret as nonlocal damage2d static data
const NuTo::ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain* NuTo::ConstitutiveStaticDataBase::AsNonlocalDamagePlasticity2DPlaneStrain()const
{
    throw NuTo::MechanicsException("[NuTo::ConstitutiveStaticDataBase::AsNonlocalDamagePlasticity2DPlaneStrain] Static data is not of type NonlocalDamagePlasticity2DPlaneStrain.");
}

//!@ brief reinterpret as multiscale2d static data
NuTo::ConstitutiveStaticDataMultiscale2DPlaneStrain* NuTo::ConstitutiveStaticDataBase::AsMultiscale2DPlaneStrain()
{
    throw NuTo::MechanicsException("[NuTo::ConstitutiveStaticDataBase::ConstitutiveStaticDataMultiscale2DPlaneStrain] Static data is not of type NonlocalDamagePlasticity2DPlaneStrain.");
}

//!@ brief reinterpret as multiscale2d static data
const NuTo::ConstitutiveStaticDataMultiscale2DPlaneStrain* NuTo::ConstitutiveStaticDataBase::AsMultiscale2DPlaneStrain()const
{
    throw NuTo::MechanicsException("[NuTo::ConstitutiveStaticDataBase::ConstitutiveStaticDataMultiscale2DPlaneStrain] Static data is not of type NonlocalDamagePlasticity2DPlaneStrain.");
}

//!@ brief reinterpret as lattice concrete 2D static data
NuTo::ConstitutiveStaticDataLatticeConcrete2D* NuTo::ConstitutiveStaticDataBase::AsConstitutiveStaticDataLatticeConcrete2D()
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveStaticDataBase::ConstitutiveStaticDataLatticeConcrete2D] Static data is not of type ConstitutiveStaticDataLatticeConcrete2D.");
}

//!@ brief reinterpret as lattice concrete 2D static data
const NuTo::ConstitutiveStaticDataLatticeConcrete2D* NuTo::ConstitutiveStaticDataBase::AsConstitutiveStaticDataLatticeConcrete2D()const
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveStaticDataBase::ConstitutiveStaticDataLatticeConcrete2D] Static data is not of type ConstitutiveStaticDataLatticeConcrete2D.");
}

void NuTo::ConstitutiveStaticDataBase::SetFineScaleModel(std::string rFileName, double rMacroLength, double rCenter[2], std::string rIPName)
{
    throw NuTo::MechanicsException("[NuTo::ConstitutiveStaticDataBase::SetFineScaleModel] Static data has no fine scale model");
}

//! @brief sets the fine scale model parameters
void NuTo::ConstitutiveStaticDataBase::SetFineScaleParameter(const std::string& rName, double rParameter)
{
    throw NuTo::MechanicsException("[NuTo::ConstitutiveStaticDataBase::SetFineScaleParameter] Static data has no fine scale model");
}

//! @brief sets the fine scale model parameters
void NuTo::ConstitutiveStaticDataBase::SetFineScaleParameter(const std::string& rName, std::string rParameter)
{
    throw NuTo::MechanicsException("[NuTo::ConstitutiveStaticDataBase::SetFineScaleParameter] Static data has no fine scale model");
}

#ifdef ENABLE_VISUALIZE
//Visualize for all integration points the fine scale structure
void NuTo::ConstitutiveStaticDataBase::VisualizeIpMultiscale(VisualizeUnstructuredGrid& rVisualize,
		const boost::ptr_list<NuTo::VisualizeComponentBase>& rWhat, bool rVisualizeDamage)const
{
    //do nothing if not multiscale static data
}
#endif



