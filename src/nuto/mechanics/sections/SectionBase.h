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

    //! @brief... get if temperatures are dof
    //! @return ... true, if temperatures are dofs
    bool GetIsTemperatureDof()const;

    //! @brief... set if temperatures are dofs
    //! @param rFlag ... true, if temperatures are dofs
    void SetIsTemperatureDof(bool rFlag);

    //! @brief... get if temperatures are to be used as input to the constitutive model
    //! @return ... true, if temperatures are to be used as input to the constitutive model
    bool GetInputConstitutiveIsTemperature()const;

    //! @brief... set if temperatures are to be used as input to the constitutive model
    //! @param rFlag ... true, if temperatures are to be used as input to the constitutive model
    void SetInputConstitutiveIsTemperature(bool rFlag);

    //! @brief... get if temperature gradients are to be used as input to the constitutive model
    //! @return ... true, if temperature gradients are to be used as input to the constitutive model
    bool GetInputConstitutiveIsTemperatureGradient()const;

    //! @brief... set if temperature gradients are to be used as input to the constitutive model
    //! @param rFlag ... true, if temperature gradients are to be used as input to the constitutive model
    void SetInputConstitutiveIsTemperatureGradient(bool rFlag);

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

    //! @brief... get if deformation Gradient is to be used as input to the constitutive model
    //! @return ... true, if deformation Gradient are to be used as input to the constitutive model
    bool GetInputConstitutiveIsDeformationGradient()const;

    //! @brief... set if deformation Gradient are to be used as input to the constitutive model
    //! @param rFlag ... true, if deformation Gradient are to be used as input to the constitutive model
    void SetInputConstitutiveIsDeformationGradient(bool rFlag);

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
    bool mInputConstitutiveIsTemperatureGradient;
    bool mInputConstitutiveIsDeformationGradient;

    bool mDisplacementDof;
    bool mTemperatureDof;
};

}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::SectionBase)
#endif // ENABLE_SERIALIZATION

#endif // SECTIONBASE_H
