// $Id$

#ifndef CONSTITUTIVEINPUTBASE_H_
#define CONSTITUTIVEINPUTBASE_H_

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include <string>

#include "nuto/base/ErrorEnum.h"
#include "nuto/math/FullMatrix_Def.h"
#include "nuto/mechanics/elements/ElementEnum.h"
#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"

namespace NuTo
{
// forward declarations
class ConstitutiveHeatFluxTemperatureGradient;
class ConstitutiveEngineeringStressStrain;
class ConstitutiveLatticeStressStrain;
class DeformationGradient1D;
class DeformationGradient2D;
class DeformationGradient3D;
class ElementBase;
class EngineeringStrain3D;
class EngineeringStress1D;
class EngineeringStress2D;
class EngineeringStress3D;
class Temperature;
class TemperatureGradient3D;
class Logger;
class SecondPiolaKirchhoffStress3D;
class StructureBase;

//! @brief ... base class for the constitutive relationship, e.g. material laws
//! @author Jörg F. Unger, BAM
//! @date July 2012
class ConstitutiveInputBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    //! @brief ... constructor
    ConstitutiveInputBase();

    //! @brief ... constructor
    virtual ~ConstitutiveInputBase()
    {}

    virtual const DeformationGradient1D& GetDeformationGradient1D()const;
    virtual const DeformationGradient2D& GetDeformationGradient2D()const;
    virtual const DeformationGradient3D& GetDeformationGradient3D()const;
    virtual double GetTemperature()const;
    virtual const TemperatureGradient3D& GetTemperatureGradient3D()const;
    virtual double GetDamage()const;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

};
}
#endif // CONSTITUTIVEINPUTBASE_H_

