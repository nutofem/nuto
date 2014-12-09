// $Id$

#ifndef SECTIONBASE_H
#define SECTIONBASE_H

#include <iostream>

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/sections/SectionEnum.h"

namespace NuTo
{
//! @author Stefan Eckardt, ISM
//! @date October 2009
//! @brief ... standard abstract base class for sections
class SectionBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    //! @brief ... constructor
    SectionBase();

    //! @brief ... destructor
    virtual ~SectionBase(){};

    //! @brief ... get the section type
    //! @return ... section type
    virtual Section::eSectionType GetType() const = 0;

    //! @brief ... returns the number of inputs in the constitutive model
    int GetNumInputConstitutive()const;

    //! @brief ... returns the number of inputs in the constitutive model
    int GetNumDofs()const;

    //! @brief... get if displacements are dof
    //! @return ... true, if displacements are dofs
    bool GetIsDisplacementDof()const;

    //! @brief... set if displacements are dofs
    //! @param rFlag ... true, if displacements are dofs
    void SetIsDisplacementDof(bool rFlag);

    //! @brief... get if rotations are dofs
    //! @return ... true, if displacements are dofs
    bool GetIsRotationDof()const;

    //! @brief... set if rotations are dofs
    //! @param rFlag ... true, if displacements are dofs
    void SetIsRotationDof(bool rFlag);

    //! @brief... get if temperatures are dof
    //! @return ... true, if temperatures are dofs
    bool GetIsTemperatureDof()const;

    //! @brief... set if temperatures are dofs
    //! @param rFlag ... true, if temperatures are dofs
    void SetIsTemperatureDof(bool rFlag);

    //! @brief... get if damage are dof
    //! @return ... true, if damage are dofs
    bool GetIsNonlocalEqPlasticStrainDof()const;

    //! @brief... set if damage are dofs
    //! @param rFlag ... true, if damage are dofs
    void SetIsNonlocalEqPlasticStrainDof(bool rFlag);

    //! @brief... get if damage are dof
    //! @return ... true, if damage are dofs
    bool GetIsNonlocalTotalStrainDof()const;

    //! @brief... set if damage are dofs
    //! @param rFlag ... true, if damage are dofs
    void SetIsNonlocalTotalStrainDof(bool rFlag);

    //! @brief... get if nonlocal eq strain is dofs
    //! @return ... true, if nonlocal eq strain is dofs
    bool GetIsNonlocalEqStrainDof()const;

    //! @brief... set if nonlocal eq strain is dofs
    //! @param rFlag ... true, if nonlocal eq strain is dofs
    void SetIsNonlocalEqStrainDof(bool rFlag);

    //! @brief... get if water phase fraction is dofs
    //! @return ... true, if water phase fraction is dofs
    bool GetIsWaterPhaseFractionDof()const;

    //! @brief... set if water phase fraction is dofs
    //! @param rFlag ... true, if water phase fraction is dofs
    void SetIsWaterPhaseFractionDof(bool rFlag);

    //! @brief... get if relative humidity is dofs
    //! @return ... true, if relative humidity is dofs
    bool GetIsRelativeHumidityDof()const;

    //! @brief... set if relative humidity is dofs
    //! @param rFlag ... true, if relative humidity is dofs
    void SetIsRelativeHumidityDof(bool rFlag);

    //! @brief... get if temperatures are to be used as input to the constitutive model
    //! @return ... true, if temperatures are to be used as input to the constitutive model
    bool GetInputConstitutiveIsTemperature()const;

    //! @brief... set if temperatures are to be used as input to the constitutive model
    //! @param rFlag ... true, if temperatures are to be used as input to the constitutive model
    void SetInputConstitutiveIsTemperature(bool rFlag);

    //! @brief... get if damage is to be used as input to the constitutive model (gradient damage formulation)
    //! @return ... true, if damage is to be used as input to the constitutive model
    bool GetInputConstitutiveIsNonlocalEqPlasticStrain()const;

    //! @brief... set if damage is to be used as input to the constitutive model (gradient damage formulation)
    //! @param rFlag ... true, if damage is to be used as input to the constitutive model
    void SetInputConstitutiveIsNonlocalEqPlasticStrain(bool rFlag);

    //! @brief... get if damage is to be used as input to the constitutive model (gradient damage formulation)
    //! @return ... true, if damage is to be used as input to the constitutive model
    bool GetInputConstitutiveIsNonlocalTotalStrain()const;

    //! @brief... set if damage is to be used as input to the constitutive model (gradient damage formulation)
    //! @param rFlag ... true, if damage is to be used as input to the constitutive model
    void SetInputConstitutiveIsNonlocalTotalStrain(bool rFlag);

    //! @brief... get if damage is to be used as input to the constitutive model (gradient damage formulation)
    //! @return ... true, if damage is to be used as input to the constitutive model
    bool GetInputConstitutiveIsNonlocalEqStrain()const;

    //! @brief... set if damage is to be used as input to the constitutive model (gradient damage formulation)
    //! @param rFlag ... true, if damage is to be used as input to the constitutive model
    void SetInputConstitutiveIsNonlocalEqStrain(bool rFlag);

    //! @brief... get if temperature gradients are to be used as input to the constitutive model
    //! @return ... true, if temperature gradients are to be used as input to the constitutive model
    bool GetInputConstitutiveIsTemperatureGradient()const;

    //! @brief... set if temperature gradients are to be used as input to the constitutive model
    //! @param rFlag ... true, if temperature gradients are to be used as input to the constitutive model
    void SetInputConstitutiveIsTemperatureGradient(bool rFlag);

    //! @brief... get if deformation Gradient is to be used as input to the constitutive model
    //! @return ... true, if deformation Gradient are to be used as input to the constitutive model
    bool GetInputConstitutiveIsDeformationGradient()const;

    //! @brief... set if deformation Gradient are to be used as input to the constitutive model
    //! @param rFlag ... true, if deformation Gradient are to be used as input to the constitutive model
    void SetInputConstitutiveIsDeformationGradient(bool rFlag);

    //! @brief... get if water phase fraction is to be used as input to the constitutive model
    //! @return ... true, if water phase fraction are to be used as input to the constitutive model
    bool GetInputConstitutiveIsWaterPhaseFraction() const;

    //! @brief... set if water phase fraction are to be used as input to the constitutive model
    //! @param rFlag ... true, if water phase fraction are to be used as input to the constitutive model
    void SetInputConstitutiveIsWaterPhaseFraction(bool rFlag);

    //! @brief... get if relative humidity is to be used as input to the constitutive model
    //! @return ... true, if relative humidity are to be used as input to the constitutive model
    bool GetInputConstitutiveIsRelativeHumidity() const;

    //! @brief... set if relative humidity are to be used as input to the constitutive model
    //! @param rFlag ... true, if relative humidity are to be used as input to the constitutive model
    void SetInputConstitutiveIsRelativeHumidity(bool rFlag);

    //! @brief ... get the cross-section area of the section
    //! @return ... section cross-section area
    virtual double GetArea() const;

    //! @brief ... set the cross-section area of the section
    //! @param rArea ... cross-section area
    virtual void SetArea(double rArea);

    //! @brief ... get the thickness of the section
    //! @return ... section thickness
    virtual double GetThickness() const;

    //! @brief ... get the thickness of the section
    //! @return ... section thickness
    virtual void SetThickness(double rThickness);

    //! @brief ... print information about the section
    //! @param rVerboseLevel ... verbosity of the information
    virtual void Info(unsigned short rVerboseLevel) const;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

#endif // ENABLE_SERIALIZATION
private:
    bool mInputConstitutiveIsTemperature;
    bool mInputConstitutiveIsNonlocalEqPlasticStrain;
    bool mInputConstitutiveIsNonlocalTotalStrain;
    bool mInputConstitutiveIsNonlocalEqStrain;
    bool mInputConstitutiveIsTemperatureGradient;
    bool mInputConstitutiveIsDeformationGradient;
    bool mInputConstitutiveIsWaterPhaseFraction = false;
    bool mInputConstitutiveIsRelativeHumidity = false;


    bool mDisplacementDof;
    bool mTemperatureDof;
    bool mNonlocalEqPlasticStrainDof;
    bool mNonlocalTotalStrainDof;
    bool mNonlocalEqStrainDof;
    bool mWaterPhaseFractionDof = false;
    bool mRelativeHumidityDoF = false;
};

}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::SectionBase)
#endif // ENABLE_SERIALIZATION

#endif // SECTIONBASE_H
