#pragma once

#include "nuto/mechanics/constitutive/ConstitutiveBase.h"

namespace NuTo
{
class ConstitutiveIOBase;
class InterpolationType;

class HeatConduction: public ConstitutiveBase
{

#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    HeatConduction();

    //! @brief Determines the constitutive inputs needed to evaluate the constitutive outputs.
    //! @param rConstitutiveOutput Desired constitutive outputs
    //! @param rInterpolationType Interpolation type to determine additional inputs
    //! @return constitutive inputs needed for the evaluation
    ConstitutiveInputMap GetConstitutiveInputs(
            const ConstitutiveOutputMap& rConstitutiveOutput,
            const InterpolationType& rInterpolationType) const override;

    //! @brief Evaluate the constitutive relation in 1D.
    //! @param rElement Element
    //! @param rIp Integration point
    //! @param rConstitutiveInput Input to the constitutive law
    //! @param rConstitutiveOutput Output to the constitutive law
    template<int TDim>
    NuTo::Error::eError Evaluate(ElementBase* rElement, int rIp,
            const ConstitutiveInputMap& rConstitutiveInput,
            const ConstitutiveOutputMap& rConstitutiveOutput);

    virtual NuTo::Error::eError Evaluate1D(ElementBase* rElement, int rIp,
            const ConstitutiveInputMap& rConstitutiveInput,
            const ConstitutiveOutputMap& rConstitutiveOutput) override
    {
        return Evaluate<1>(rElement, rIp, rConstitutiveInput, rConstitutiveOutput);
    }

    //! @brief ... evaluate the constitutive relation in 2D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
    //! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
    virtual NuTo::Error::eError Evaluate2D(ElementBase* rElement, int rIp,
            const ConstitutiveInputMap& rConstitutiveInput,
            const ConstitutiveOutputMap& rConstitutiveOutput) override
    {
        return Evaluate<2>(rElement, rIp, rConstitutiveInput, rConstitutiveOutput);
    }

    //! @brief ... evaluate the constitutive relation in 3D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
    //! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
    virtual NuTo::Error::eError Evaluate3D(ElementBase* rElement, int rIp,
            const ConstitutiveInputMap& rConstitutiveInput,
            const ConstitutiveOutputMap& rConstitutiveOutput) override
    {
        return Evaluate<3>(rElement, rIp, rConstitutiveInput, rConstitutiveOutput);
    }

    //! @brief Create new static data object for an integration point.
    //! @return Pointer to static data object
    ConstitutiveStaticDataBase* AllocateStaticData1D(const ElementBase* rElement) const override {return nullptr;}

    //! @brief Create new static data object for an integration point.
    //! @return Pointer to static data object
    ConstitutiveStaticDataBase* AllocateStaticData2D(const ElementBase* rElement) const override  {return nullptr;}

    //! @brief Create new static data object for an integration point.
    //! @return Pointer to static data object
    ConstitutiveStaticDataBase* AllocateStaticData3D(const ElementBase* rElement) const override  {return nullptr;}

    //! @brief ... determines which submatrices of a multi-doftype problem can be solved by the constitutive law
    //! @param rDofRow ... row dof
    //! @param rDofCol ... column dof
    //! @param rTimeDerivative ... time derivative
    virtual bool CheckDofCombinationComputable(Node::eDof rDofRow,
                                                Node::eDof rDofCol,
                                                int rTimeDerivative) const override;

    //! @brief Checks if the constitutive law has a specific parameter.
    //! @param rIdentifier Enum to identify the requested parameter
    virtual bool CheckHaveParameter(Constitutive::eConstitutiveParameter rIdentifier) const override;

    //! @brief Gets a parameter of the constitutive law which is selected by an enum.
    //! @param rIdentifier Enum to identify the requested parameter
    //! @return Value of the requested variable
    virtual double GetParameterDouble(Constitutive::eConstitutiveParameter rIdentifier) const override;

    //! @brief Sets a parameter of the constitutive law which is selected by an enum.
    //! @param rIdentifier Enum to identify the requested parameter
    //! @param rValue New value for requested variable
    virtual void SetParameterDouble(Constitutive::eConstitutiveParameter rIdentifier, double rValue) override;

    //! @brief Gets a set of all constitutive output enums that are compatible with the constitutive law.
    //! @return Set of all constitutive output enums that are compatible with the constitutive law
    virtual bool CheckOutputTypeCompatibility(Constitutive::Output::eOutput rOutputEnum) const override;

    //! @brief Get type of constitutive relationship.
    //! @return Type of constitutive relationship
    //! @sa eConstitutiveType
    Constitutive::eConstitutiveType GetType() const override;

    //! @brief Check parameters of the constitutive relationship.
    void CheckParameters()const;

    //! @brief Check compatibility between element type and type of constitutive relationship.
    //! @param rElementType Element type
    //! @return `true` if the element is compatible with the constitutive relationship, `false` otherwise.
    bool CheckElementCompatibility(Element::eElementType rElementType) const override;

    //! @brief Print information about the object.
    //! @param rVerboseLevel Verbosity of the information
    //! @param rLogger Stream for the output
    void Info(unsigned short rVerboseLevel, Logger& rLogger) const override;

    //! @brief Returns true, if a material model has tmp static data (which has to be updated before stress or stiffness are calculated).
    bool HaveTmpStaticData() const override
    {
    	return false;
    }

#ifdef ENABLE_SERIALIZATION
    //! @brief Serializes the class.
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION


protected:
    //! @brief Thermal conduction coefficient \f$ k \f$
    double mK;

    //! @brief Specific heat capacity \f$ c_T \f$
    double mCt;

    //! @brief Density \f$ \rho \f$
    double mRho;

    template <int TDim>
    struct InputData
    {
        Eigen::Matrix<double, TDim, 1> mTemperatureGradient;
        double mTemperatureChange;
    };

};

}

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::HeatConduction)
#endif //ENABLE_SERIALIZATION
