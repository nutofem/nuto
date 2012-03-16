// $Id$

#ifndef CONSTITUTIVEPIOLAKIRCHHOFFIIGREENLAGRANGE_H_
#define CONSTITUTIVEPIOLAKIRCHHOFFIIGREENLAGRANGE_H_

#include "nuto/mechanics/constitutive/ConstitutiveBase.h"

namespace NuTo
{
// forward declarations
class StructureBase;
class ElementBase;
class ConstitutiveStaticDataBase;
class ConstitutiveTangentLocal1x1;
class ConstitutiveTangentLocal3x3;
class ConstitutiveTangentLocal6x6;
class DeformationGradient1D;
class DeformationGradient2D;
class DeformationGradient3D;
class SecondPiolaKirchhoffStress1D;
class SecondPiolaKirchhoffStress2D;
class SecondPiolaKirchhoffStress3D;

//! @brief ... base class for the mechanical constitutive relationship using Piola-Kirch II stresses and Green Lagrange strains
//! @author Jörg F. Unger, ISM
//! @date December 2009
class ConstitutivePiolaKirchhoffIIGreenLagrange : virtual public ConstitutiveBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    //! @brief ... constructor
    ConstitutivePiolaKirchhoffIIGreenLagrange();

    // Second PiolaKirchhoff Stress - Green Lagrange Strains /////////////////////////////////////
    //! @brief ... calculate second Piola-Kirchhoff stress from GreenLagrangeStrains (which is calculated from the deformation gradient)
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rSecondPiolaKirchhoffStress ... Second Piola Kirchhoff Stress
    virtual Error::eError GetSecondPiolaKirchhoffStressFromGreenLagrangeStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient1D& rDeformationGradient, SecondPiolaKirchhoffStress1D& rSecondPiolaKirchhoffStress) const=0;

    // Second PiolaKirchhoff Stress - Green Lagrange Strains /////////////////////////////////////
    //! @brief ... calculate second Piola-Kirchhoff stress from GreenLagrangeStrains (which is calculated from the deformation gradient)
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rSecondPiolaKirchhoffStress ... Second Piola Kirchhoff Stress
    virtual Error::eError GetSecondPiolaKirchhoffStressFromGreenLagrangeStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient2D& rDeformationGradient, SecondPiolaKirchhoffStress2D& rSecondPiolaKirchhoffStress) const=0;

    // Second PiolaKirchhoff Stress - Green Lagrange Strains /////////////////////////////////////
    //! @brief ... calculate second Piola-Kirchhoff stress from GreenLagrangeStrains (which is calculated from the deformation gradient)
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rSecondPiolaKirchhoffStress ... Second Piola Kirchhoff Stress
    virtual Error::eError GetSecondPiolaKirchhoffStressFromGreenLagrangeStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient3D& rDeformationGradient, SecondPiolaKirchhoffStress3D& rSecondPiolaKirchhoffStress) const=0;

    //! @brief ... calculate the tangent (derivative of the Engineering stresses with respect to the engineering strains) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rTangent ... tangent
    virtual Error::eError GetTangent_SecondPiolaKirchhoffStress_GreenLagrangeStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient1D& rDeformationGradient,
            ConstitutiveTangentLocal1x1& rTangent) const=0;

    //! @brief ... calculate the tangent (derivative of the Engineering stresses with respect to the engineering strains) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rTangent ... tangent
    virtual Error::eError GetTangent_SecondPiolaKirchhoffStress_GreenLagrangeStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient2D& rDeformationGradient,
            ConstitutiveTangentLocal3x3& rTangent) const=0;

    //! @brief ... calculate the tangent (derivative of the Engineering stresses with respect to the engineering strains) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rTangent ... tangent
    virtual Error::eError GetTangent_SecondPiolaKirchhoffStress_GreenLagrangeStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient3D& rDeformationGradient,
            ConstitutiveTangentLocal6x6& rTangent) const=0;

    //! @brief ... update static data (history variables) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    virtual Error::eError UpdateStaticData_SecondPiolaKirchhoffStress_GreenLagrangeStrain(ElementBase* rElement, int rIp,
            const DeformationGradient1D& rDeformationGradient) const=0;

    //! @brief ... update static data (history variables) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    virtual Error::eError UpdateStaticData_SecondPiolaKirchhoffStress_GreenLagrangeStrain(ElementBase* rElement, int rIp,
            const DeformationGradient2D& rDeformationGradient) const=0;

    //! @brief ... update static data (history variables) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    virtual Error::eError UpdateStaticData_SecondPiolaKirchhoffStress_GreenLagrangeStrain(ElementBase* rElement, int rIp,
            const DeformationGradient3D& rDeformationGradient) const=0;

    //! @brief ... create new static data object for an integration point
    //! @return ... pointer to static data object
    virtual ConstitutiveStaticDataBase* AllocateStaticDataSecondPiolaKirchhoffStress_GreenLagrangeStrain1D(const ElementBase* rElement) const=0;

    //! @brief ... create new static data object for an integration point
    //! @return ... pointer to static data object
    virtual ConstitutiveStaticDataBase* AllocateStaticDataSecondPiolaKirchhoffStress_GreenLagrangeStrain2D(const ElementBase* rElement) const=0;

    //! @brief ... create new static data object for an integration point
    //! @return ... pointer to static data object
    virtual ConstitutiveStaticDataBase* AllocateStaticDataSecondPiolaKirchhoffStress_GreenLagrangeStrain3D(const ElementBase* rElement) const=0;

    //! @brief ... calculate the total energy density
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    virtual Error::eError GetInternalEnergy_SecondPiolaKirchhoffStress_GreenLagrangeStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient1D& rDeformationGradient, double& rEnergy) const=0;

    //! @brief ... calculate the total energy density
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    virtual Error::eError GetInternalEnergy_SecondPiolaKirchhoffStress_GreenLagrangeStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient2D& rDeformationGradient, double& rEnergy) const=0;

    //! @brief ... calculate the total energy density
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    virtual Error::eError GetInternalEnergy_SecondPiolaKirchhoffStress_GreenLagrangeStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient3D& rDeformationGradient, double& rEnergy) const=0;

    //! @brief ... calculate the elastic energy density
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    virtual Error::eError GetElasticEnergy_SecondPiolaKirchhoffStress_GreenLagrangeStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient1D& rDeformationGradient, double& rEnergy) const=0;

    //! @brief ... calculate the elastic energy density
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    virtual Error::eError GetElasticEnergy_SecondPiolaKirchhoffStress_GreenLagrangeStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient2D& rDeformationGradient, double& rEnergy) const=0;

    //! @brief ... calculate the elastic energy density
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    virtual Error::eError GetElasticEnergy_SecondPiolaKirchhoffStress_GreenLagrangeStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient3D& rDeformationGradient, double& rEnergy) const=0;
    ///////////////////////////////////////////////////////////////////////////

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION
protected:

};

}

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ConstitutivePiolaKirchhoffIIGreenLagrange)
#endif // ENABLE_SERIALIZATION

#endif // CONSTITUTIVEPIOLAKIRCHHOFFIIGREENLAGRANGE_H_
