// $Id$

#ifndef CONSTITUTIVEENGINEERINGSTRESSSTRAIN_H_
#define CONSTITUTIVEENGINEERINGSTRESSSTRAIN_H_

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"

namespace NuTo
{
// forward declarations
class ElementBase;
class ConstitutiveTangentBase;
class DeformationGradient1D;
class DeformationGradient2D;
class DeformationGradient3D;
class EngineeringStress1D;
class EngineeringStress2D;
class EngineeringStress3D;
class EngineeringStrain1D;
class EngineeringStrain2D;
class EngineeringStrain3D;

//! @brief ... base class for the mechanical constitutive relationship using engineering stress and strains
//! @author Jörg F. Unger, ISM
//! @date December 2009
class ConstitutiveEngineeringStressStrain : public virtual ConstitutiveBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:

    //! @brief ... constructor
    ConstitutiveEngineeringStressStrain();

    //  Engineering strain /////////////////////////////////////
    //! @brief ... calculate engineering strain from deformation gradient in 1D (truss is assumed to be plane stress)
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rEngineeringStrain ... engineering strain
    Error::eError GetEngineeringStrain(const ElementBase* rElement, int rIp,
                                      const DeformationGradient1D& rDeformationGradient, EngineeringStrain3D& rEngineeringStrain) const;

    //  Engineering strain /////////////////////////////////////
    //! @brief ... calculate engineering strain from deformation gradient in 2D (e.g. plane stress)
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rEngineeringStrain ... engineering strain
    Error::eError GetEngineeringStrain(const ElementBase* rElement, int rIp,
                                      const DeformationGradient2D& rDeformationGradient, EngineeringStrain3D& rEngineeringStrain) const;

    //  Engineering strain /////////////////////////////////////
    //! @brief ... calculate engineering strain from deformation gradient in 3D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rEngineeringStrain ... engineering strain
    Error::eError GetEngineeringStrain(const ElementBase* rElement, int rIp,
                                      const DeformationGradient3D& rDeformationGradient, EngineeringStrain3D& rEngineeringStrain) const;

    //  Engineering strain /////////////////////////////////////
    //! @brief ... convert an engineering strain from 2D into 3D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rEngineeringStrain2D ... engineering strain 2D (input)
    //! @param rEngineeringStrain3D ... engineering strain 3D (output)
    virtual Error::eError GetEngineeringStrainFromEngineeringStrain(const ElementBase* rElement, int rIp,
                                      const EngineeringStrain2D& rEngineeringStrain2D, EngineeringStrain3D& rEngineeringStrain3D) const;

    //  Engineering strain /////////////////////////////////////
    //! @brief ... calculate engineering plastic strain from deformation gradient in 3D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rEngineeringStrain ... engineering strain
    virtual Error::eError GetEngineeringPlasticStrain(const ElementBase* rElement, int rIp,
                                      const DeformationGradient1D& rDeformationGradient, EngineeringStrain3D& rEngineeringPlasticStrain) const=0;

    //  Engineering strain /////////////////////////////////////
    //! @brief ... calculate engineering plastic strain from deformation gradient in 3D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rEngineeringStrain ... engineering strain
    virtual Error::eError GetEngineeringPlasticStrain(const ElementBase* rElement, int rIp,
                                      const DeformationGradient2D& rDeformationGradient, EngineeringStrain3D& rEngineeringPlasticStrain) const=0;

    //  Engineering strain /////////////////////////////////////
    //! @brief ... calculate engineering plastic strain from deformation gradient in 3D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rEngineeringStrain ... engineering strain
    virtual Error::eError GetEngineeringPlasticStrain(const ElementBase* rElement, int rIp,
                                      const DeformationGradient3D& rDeformationGradient, EngineeringStrain3D& rEngineeringPlasticStrain) const=0;

    // Engineering stress - Engineering strain /////////////////////////////////////
    //! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rEngineeringStress ... Engineering stress
    virtual Error::eError GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient1D& rDeformationGradient, EngineeringStress1D& rEngineeringStress) const=0;

    // Engineering stress - Engineering strain /////////////////////////////////////
    //! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rEngineeringStress ... Engineering stress
    virtual Error::eError GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient1D& rDeformationGradient, EngineeringStress3D& rEngineeringStress) const=0;

    // Engineering stress - Engineering strain /////////////////////////////////////
    //! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rEngineeringStress ... Engineering stress
    virtual Error::eError GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient2D& rDeformationGradient, EngineeringStress2D& rEngineeringStress) const=0;

    // Engineering stress - Engineering strain /////////////////////////////////////
    //! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rEngineeringStress ... Engineering stress
    virtual Error::eError GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient2D& rDeformationGradient, EngineeringStress3D& rEngineeringStress) const=0;

    // Engineering stress - Engineering strain /////////////////////////////////////
    //! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rCauchyStress ... Cauchy stress
    virtual Error::eError GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient3D& rDeformationGradient, EngineeringStress3D& rEngineeringStress) const=0;

    //  Damage /////////////////////////////////////
    //! @brief ... calculate isotropic damage from deformation gradient in 1D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rDamage ... damage variable
    virtual Error::eError GetDamage(const ElementBase* rElement, int rIp,
                                      const DeformationGradient1D& rDeformationGradient, double& rDamage) const=0;

    //  Damage /////////////////////////////////////
    //! @brief ... calculate isotropic damage from deformation gradient in 2D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rDamage ... damage variable
    virtual Error::eError GetDamage(const ElementBase* rElement, int rIp,
                                      const DeformationGradient2D& rDeformationGradient, double& rDamage) const=0;

    //  Damage /////////////////////////////////////
    //! @brief ... calculate isotropic damage from deformation gradient in 3D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rDamage ... damage variable
    virtual Error::eError GetDamage(const ElementBase* rElement, int rIp,
                                      const DeformationGradient3D& rDeformationGradient, double& rDamage) const=0;

    //! @brief ... calculate the tangent (derivative of the Engineering stresses with respect to the engineering strains) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rTangent ... tangent
    virtual Error::eError GetTangent_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient1D& rDeformationGradient,
            ConstitutiveTangentBase* rTangent) const=0;

    //! @brief ... calculate the tangent (derivative of the Engineering stresses with respect to the engineering strains) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rTangent ... tangent
    virtual Error::eError GetTangent_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient2D& rDeformationGradient,
            ConstitutiveTangentBase* rTangent) const=0;

    //! @brief ... calculate the tangent (derivative of the Engineering stresses with respect to the engineering strains) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rTangent ... tangent
    virtual Error::eError GetTangent_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient3D& rDeformationGradient,
            ConstitutiveTangentBase* rTangent) const=0;

    //! @brief ... update static data (history variables) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    virtual Error::eError UpdateStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
            const DeformationGradient1D& rDeformationGradient) const=0;

    //! @brief ... update static data (history variables) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    virtual Error::eError UpdateStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
            const DeformationGradient2D& rDeformationGradient) const=0;

    //! @brief ... update static data (history variables) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    virtual Error::eError UpdateStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
            const DeformationGradient3D& rDeformationGradient) const=0;

    //! @brief ... update tmp static data (history variables) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    virtual Error::eError UpdateTmpStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
            const DeformationGradient1D& rDeformationGradient) const=0;

    //! @brief ... update tmp static data (history variables) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    virtual Error::eError UpdateTmpStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
            const DeformationGradient2D& rDeformationGradient) const=0;

    //! @brief ... update tmp static data (history variables) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    virtual Error::eError UpdateTmpStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
            const DeformationGradient3D& rDeformationGradient) const=0;

    //! @brief ... checks, if a model has to be switched from linear to nonlinear, and then performs the adaption
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    virtual Error::eError MultiscaleSwitchToNonlinear(ElementBase* rElement, int rIp,
            const DeformationGradient2D& rDeformationGradient)const;

    //! @brief ... create new static data object for an integration point
    //! @return ... pointer to static data object
    virtual ConstitutiveStaticDataBase* AllocateStaticDataEngineeringStress_EngineeringStrain1D(const ElementBase* rElement) const=0;

    //! @brief ... create new static data object for an integration point
    //! @return ... pointer to static data object
    virtual ConstitutiveStaticDataBase* AllocateStaticDataEngineeringStress_EngineeringStrain2D(const ElementBase* rElement) const=0;

    //! @brief ... create new static data object for an integration point
    //! @return ... pointer to static data object
    virtual ConstitutiveStaticDataBase* AllocateStaticDataEngineeringStress_EngineeringStrain3D(const ElementBase* rElement) const=0;

    //! @brief ... calculate the total energy density
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    virtual Error::eError GetInternalEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient1D& rDeformationGradient, double& energy) const;

    //! @brief ... calculate the total energy density
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    virtual Error::eError GetInternalEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient2D& rDeformationGradient, double& energy) const;

    //! @brief ... calculate the total energy density
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    virtual Error::eError GetInternalEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient3D& rDeformationGradient, double& energy) const;

    //! @brief ... calculate the elastic energy density
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    virtual Error::eError GetElasticEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient1D& rDeformationGradient, double& energy) const;

    //! @brief ... calculate the elastic energy density
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    virtual Error::eError GetElasticEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient2D& rDeformationGradient, double& energy) const;

    //! @brief ... calculate the elastic energy density
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    virtual Error::eError GetElasticEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient3D& rDeformationGradient, double& energy) const;

    //! @brief ... calculates the difference of the elastic strain between the current state and the previous update
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rDeltaElasticEngineeringStrain ... delta elastic engineering strain (return value)
    virtual void GetDeltaElasticEngineeringStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient1D& rDeformationGradient, EngineeringStrain1D& rDeltaElasticEngineeringStrain) const=0;

    //! @brief ... calculates the difference of the elastic strain between the current state and the previous update
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rDeltaElasticEngineeringStrain ... delta elastic engineering strain (return value)
    virtual void GetDeltaElasticEngineeringStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient2D& rDeformationGradient, EngineeringStrain2D& rDeltaElasticEngineeringStrain) const=0;

    //! @brief ... calculates the difference of the elastic strain between the current state and the previous update
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rDeltaElasticEngineeringStrain ... delta elastic engineering strain (return value)
    virtual void GetDeltaElasticEngineeringStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient3D& rDeformationGradient, EngineeringStrain3D& rDeltaElasticEngineeringStrain) const=0;

    //! @brief ... avoid dynamic cast
    //! @return ... see brief explanation
    virtual const ConstitutiveEngineeringStressStrain* AsConstitutiveEngineeringStressStrain()const;

    //! @brief ... avoid dynamic cast
    //! @return ... see brief explanation
    virtual ConstitutiveEngineeringStressStrain* AsConstitutiveEngineeringStressStrain();

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION
protected:
    //! @brief ... if true, the static data should be derived from StaticDataPrevStressStrain
    bool mEnergyFlag;

};

}

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ConstitutiveEngineeringStressStrain)
#endif // ENABLE_SERIALIZATION


#endif // CONSTITUTIVEENGINEERINGSTRESSSTRAIN_H_
