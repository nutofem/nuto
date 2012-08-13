// $Id

#ifndef LINEARHEATFLUX_H_
#define LINEARHEATFLUX_H_

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
class TemperatureGradient1D;
class TemperatureGradient2D;
class TemperatureGradient3D;
class HeatFlux1D;
class HeatFlux2D;
class HeatFlux3D;

//! @brief ... base class for the constitutive relationship between heat flux and temperature gradient
//! @author Jörg F. Unger, ISM
//! @date December 2009
class LinearHeatFlux : public virtual ConstitutiveBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:

    //! @brief ... constructor
    LinearHeatFlux();

    //! @brief ... evaluate the constitutive relation in 1D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
    //! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
    NuTo::Error::eError Evaluate3D(ElementBase* rElement, int rIp,
    		const std::map<NuTo::Constitutive::eInput, const ConstitutiveInputBase*>& rConstitutiveInput,
    		std::map<NuTo::Constitutive::eOutput, ConstitutiveOutputBase*>& rConstitutiveOutput) override;

    //! @brief ... allocate the correct static data
    //! @return ... see brief explanation
    NuTo::ConstitutiveStaticDataBase* AllocateStaticDataEngineeringStress_EngineeringStrain3D(ElementBase* rElement)const;

    //! @brief ... get type of constitutive relationship
    //! @return ... type of constitutive relationship
    //! @sa eConstitutiveType
    NuTo::Constitutive::eConstitutiveType GetType() const;

    //! @brief ... check compatibility between element type and type of constitutive relationship
    //! @param rElementType ... element type
    //! @return ... <B>true</B> if the element is compatible with the constitutive relationship, <B>false</B> otherwise.
    bool CheckElementCompatibility(NuTo::Element::eElementType rElementType) const;

    // check parameters
    void CheckParameters()const;

    //! @brief ... returns true, if a material model has tmp static data (which has to be updated before stress or stiffness are calculated)
    //! @return ... see brief explanation
    bool HaveTmpStaticData() const
    {
    	return false;
    }

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION
protected:
    //! @brief ... get HeatCapacity
    //! @return ... HeatCapacity
    double GetHeatCapacity() const;

    //! @brief ... set HeatCapacity
    //! @param rE ... HeatCapacity
    void SetHeatCapacity(double rHeatCapacity);

    //! @brief ... get ThermalConductivity
    //! @return ... ThermalConductivity
    double GetThermalConductivity() const;

    //! @brief ... set density
    //! @param rRho ... density
    void SetThermalConductivity(double rThermalConductivity);

    //! @brief ... check if heat capacity is positive
    //! @param rHeatCapacity ... heat capacity
    void CheckHeatCapacity(double rHeatCapacity) const;

    //! @brief ... check if thermal conductivity is positive
    //! @param rThermalConductivity ... thermal conductivity
    void CheckThermalConductivity(double rThermalConductivity) const;

    //! @brief ... heat capacity (per volume)
    double mHeatCapacity;

    //! @brief ... thermal conductivity
    double mThermalConductivity;


};

}

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::LinearHeatFlux)
#endif // ENABLE_SERIALIZATION


#endif // LINEARHEATFLUX_H_
