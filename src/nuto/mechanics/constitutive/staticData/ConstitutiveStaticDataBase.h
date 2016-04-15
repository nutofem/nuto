// $Id$

#ifndef CONSTITUTIVESTATICDATABASE_H_
#define CONSTITUTIVESTATICDATABASE_H_

#include <string>
#include "nuto/math/FullMatrix_Def.h"
#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/elements/ElementEnum.h"

#ifdef ENABLE_VISUALIZE
#include <boost/ptr_container/ptr_list.hpp>
#endif // ENABLE_VISUALIZE

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

//! @brief ... base class, storing the static data (history variables) of a constitutive relationship
//! @author Stefan Eckardt, ISM
//! @date November 2009
namespace NuTo
{
class ConstitutiveStaticDataLatticeConcrete2D;
class ConstitutiveStaticDataGradientDamagePlasticity1D;
class ConstitutiveStaticDataGradientDamage;
class ConstitutiveStaticDataGradientDamage1DFatigue;
class ConstitutiveStaticDataGradientDamage2DFatigue;
class ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain;
template <int TDim> class ConstitutiveStaticDataMisesPlasticity;
class ConstitutiveStaticDataMultiscale2DPlaneStrain;
class ConstitutiveStaticDataStrainGradientDamagePlasticity1D;
class ConstitutiveStaticDataDamageViscoPlasticity3D;
class ConstitutiveStaticDataDamageViscoPlasticity3DFatigue;
class ConstitutiveStaticDataMoistureTransport;
class ConstitutiveStaticDataMultipleConstitutiveLaws;
class ConstitutiveStaticDataBondStressSlip;
class ElementBase;
class VisualizeUnstructuredGrid;
class VisualizeComponentBase;

class ConstitutiveStaticDataBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
    template<class Archive> void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION
public:
    //! @brief constructor
    ConstitutiveStaticDataBase()
    {}

    //! @brief destructor (virtual, in order to make the class a polymorphic type)
    virtual ~ConstitutiveStaticDataBase()
    {};

    //! @brief constructor
    virtual ConstitutiveStaticDataBase* Clone()const=0;

    NuTo::ConstitutiveStaticDataBase& operator= (NuTo::ConstitutiveStaticDataBase const& rOther) = default;

    //! @brief check, if the static data is compatible with a given element and a given constitutive model
    virtual bool CheckConstitutiveCompatibility(NuTo::Constitutive::eConstitutiveType constitutiveType, NuTo::Element::eElementType elementType)const=0;

    //! @brief sets the fine scale model (deserialization from a binary file)
    virtual void SetFineScaleModel(std::string rFileName, double rMacroLength, double rCenter[2], std::string rIPName);

    //! @brief sets the fine scale model (deserialization from a binary file)
    virtual void SetFineScaleParameter(const std::string& rName, double rParameter);

    //! @brief sets the fine scale model (deserialization from a binary file)
    virtual void SetFineScaleParameter(const std::string& rName, std::string rParameter);

    //!@ brief reinterpret as nonlocal damage2d static data
    virtual ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain* AsNonlocalDamagePlasticity2DPlaneStrain();

    //!@ brief reinterpret as nonlocal damage2d static data
    virtual const ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain* AsNonlocalDamagePlasticity2DPlaneStrain()const;

    //!@ brief reinterpret as gradient damage 1d static data
    virtual ConstitutiveStaticDataGradientDamagePlasticity1D* AsGradientDamagePlasticity1D();

    //!@ brief reinterpret as gradient damage1d static data
    virtual const ConstitutiveStaticDataGradientDamagePlasticity1D* AsGradientDamagePlasticity1D()const;

    //!@ brief reinterpret as gradient damage 1d static data
    virtual ConstitutiveStaticDataGradientDamage* AsGradientDamage();

    //!@ brief reinterpret as gradient damage 1d static data
    virtual const ConstitutiveStaticDataGradientDamage* AsGradientDamage() const;

    //!@ brief reinterpret as gradient damage 1d static data with fatigue
    virtual ConstitutiveStaticDataGradientDamage1DFatigue* AsGradientDamage1DFatigue();

    //!@ brief reinterpret as gradient damage 1d static data with fatigue
    virtual const ConstitutiveStaticDataGradientDamage1DFatigue* AsGradientDamage1DFatigue() const;

    //!@ brief reinterpret as gradient damage 2d static data with fatigue
    virtual ConstitutiveStaticDataGradientDamage2DFatigue* AsGradientDamage2DFatigue();

    //!@ brief reinterpret as gradient damage 2d static data with fatigue
    virtual const ConstitutiveStaticDataGradientDamage2DFatigue* AsGradientDamage2DFatigue() const;

    //!@ brief reinterpret as bond stress slip static data
    virtual ConstitutiveStaticDataBondStressSlip* AsBondStressSlip();

    //!@ brief reinterpret as bond stress slip static data
    virtual const ConstitutiveStaticDataBondStressSlip* AsBondStressSlip() const;

    //!@ brief reinterpret as lattice concrete 2D static data
    virtual NuTo::ConstitutiveStaticDataLatticeConcrete2D* AsConstitutiveStaticDataLatticeConcrete2D();

    //!@ brief reinterpret as lattice concrete 2D static data
    virtual const NuTo::ConstitutiveStaticDataLatticeConcrete2D* AsConstitutiveStaticDataLatticeConcrete2D()const;

    //!@ brief reinterpret as nonlocal damage2d static data
    virtual ConstitutiveStaticDataMultiscale2DPlaneStrain* AsMultiscale2DPlaneStrain();

    //!@ brief reinterpret as nonlocal damage2d static data
    virtual const ConstitutiveStaticDataMultiscale2DPlaneStrain* AsMultiscale2DPlaneStrain()const;

    virtual ConstitutiveStaticDataMisesPlasticity<1>* AsConstitutiveStaticDataMisesPlasticity1D()
    {throw NuTo::MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] Static data is not of type ConstitutiveStaticDataMisesPlasticity 1D.");}

    virtual const ConstitutiveStaticDataMisesPlasticity<1>* AsConstitutiveStaticDataMisesPlasticity1D() const
    {throw NuTo::MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] Static data is not of type ConstitutiveStaticDataMisesPlasticity 1D.");}

    virtual ConstitutiveStaticDataMisesPlasticity<2>* AsConstitutiveStaticDataMisesPlasticity2D()
    {throw NuTo::MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] Static data is not of type ConstitutiveStaticDataMisesPlasticity 2D.");}

    virtual const ConstitutiveStaticDataMisesPlasticity<2>* AsConstitutiveStaticDataMisesPlasticity2D() const
    {throw NuTo::MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] Static data is not of type ConstitutiveStaticDataMisesPlasticity 2D.");}

    virtual ConstitutiveStaticDataMisesPlasticity<3>* AsConstitutiveStaticDataMisesPlasticity3D()
    {throw NuTo::MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] Static data is not of type ConstitutiveStaticDataMisesPlasticity 3D.");}

    virtual const ConstitutiveStaticDataMisesPlasticity<3>* AsConstitutiveStaticDataMisesPlasticity3D() const
    {throw NuTo::MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] Static data is not of type ConstitutiveStaticDataMisesPlasticity 3D.");}


    //!@ brief reinterpret as strain gradient damage 1d static data
    virtual ConstitutiveStaticDataStrainGradientDamagePlasticity1D* AsStrainGradientDamagePlasticity1D();

    //!@ brief reinterpret as strain gradient damage1d static data
    virtual const ConstitutiveStaticDataStrainGradientDamagePlasticity1D* AsStrainGradientDamagePlasticity1D()const;

    //!@ brief reinterpret as damage viscoplasticity static data
    virtual ConstitutiveStaticDataDamageViscoPlasticity3D* AsDamageViscoPlasticity3D();

    //!@ brief reinterpret as damage viscoplasticity static data
    virtual const ConstitutiveStaticDataDamageViscoPlasticity3D* AsDamageViscoPlasticity3D()const;

    //!@ brief reinterpret as damage viscoplasticity static data with fatigue
    virtual ConstitutiveStaticDataDamageViscoPlasticity3DFatigue* AsDamageViscoPlasticity3DFatigue();

    //!@ brief reinterpret as damage viscoplasticity static data
    virtual const ConstitutiveStaticDataDamageViscoPlasticity3DFatigue* AsDamageViscoPlasticity3DFatigue()const;

    //!@ brief reinterpret as moisture transport
    virtual ConstitutiveStaticDataMoistureTransport* AsMoistureTransport();

    //!@ brief reinterpret as moisture transport
    virtual const ConstitutiveStaticDataMoistureTransport* AsMoistureTransport()const;

    //!@ brief reinterpret as multi physics
    virtual ConstitutiveStaticDataMultipleConstitutiveLaws* AsMultipleConstitutiveLaws();

    //!@ brief reinterpret as multi physics
    virtual const ConstitutiveStaticDataMultipleConstitutiveLaws* AsMultipleConstitutiveLaws()const;

#ifdef ENABLE_VISUALIZE
    //Visualize for all integration points the fine scale structure
    virtual void VisualizeIpMultiscale(VisualizeUnstructuredGrid& rVisualize,
    		const boost::ptr_list<NuTo::VisualizeComponentBase>& rWhat, bool rVisualizeDamage)const;
#endif


};

}

#endif // CONSTITUTIVESTATICDATABASE_H_ 
