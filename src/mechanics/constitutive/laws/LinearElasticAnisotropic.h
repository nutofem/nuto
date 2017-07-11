#pragma once


#include "mechanics/constitutive/ConstitutiveBase.h"


namespace NuTo
{
namespace Constitutive
{
class IPConstitutiveLawBase;
} // namespace Constitutive
class InterpolationType;

//! @brief ... linear elastic material model fully anisotropic
/*!
  */
//!
//!
class LinearElasticAnisotropic : public ConstitutiveBase
{

public:
    LinearElasticAnisotropic();

    std::unique_ptr<Constitutive::IPConstitutiveLawBase> CreateIPLaw() override;


    //! @brief ... evaluate the constitutive relation in 1D
    //! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
    //! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
    template <int TDim>
    void Evaluate(const ConstitutiveInputMap& rConstitutiveInput, const ConstitutiveOutputMap& rConstitutiveOutput);


    ConstitutiveInputMap GetConstitutiveInputs(const ConstitutiveOutputMap& rConstitutiveOutput) const override;

    //! @brief ... determines which submatrices of a multi-doftype problem can be solved by the constitutive law
    //! @param rDofRow ... row dof
    //! @param rDofCol ... column dof
    //! @param rTimeDerivative ... time derivative
    virtual bool CheckDofCombinationComputable(Node::eDof rDofRow, Node::eDof rDofCol,
                                               int rTimeDerivative) const override;

    // parameters /////////////////////////////////////////////////////////////

    //! @brief ... checks if the constitutive law has a specific parameter
    //! @param rIdentifier ... Enum to identify the requested parameter
    //! @return ... true/false
    virtual bool CheckHaveParameter(Constitutive::eConstitutiveParameter rIdentifier) const override;

    //! @brief ... gets a parameter of the constitutive law which is selected by an enum
    //! @param rIdentifier ... Enum to identify the requested parameter
    //! @return ... value of the requested variable
    virtual double GetParameterDouble(Constitutive::eConstitutiveParameter rIdentifier) const override;

    //! @brief ... sets a parameter of the constitutive law which is selected by an enum
    //! @param rIdentifier ... Enum to identify the requested parameter
    //! @param rValue ... new value for requested variable
    virtual void SetParameterDouble(Constitutive::eConstitutiveParameter rIdentifier, double rValue) override;

    //! @brief ... gets a parameter of the constitutive law which is selected by an enum
    //! @param rIdentifier ... Enum to identify the requested parameter
    //! @return ... value of the requested variable
    virtual Eigen::MatrixXd GetParameterMatrixDouble(Constitutive::eConstitutiveParameter rIdentifier) const override;

    //! @brief ... sets a parameter of the constitutive law which is selected by an enum
    //! @param rIdentifier ... Enum to identify the requested parameter
    //! @param rValue ... new value for requested variable
    virtual void SetParameterMatrixDouble(Constitutive::eConstitutiveParameter rIdentifier,
                                          Eigen::MatrixXd rValue) override;

    ///////////////////////////////////////////////////////////////////////////


    //! @brief ... gets a set of all constitutive output enums that are compatible with the constitutive law
    //! @return ... set of all constitutive output enums that are compatible with the constitutive law
    virtual bool CheckOutputTypeCompatibility(Constitutive::eOutput rOutputEnum) const override;

    //! @brief ... get type of constitutive relationship
    //! @return ... type of constitutive relationship
    //! @sa eConstitutiveType
    Constitutive::eConstitutiveType GetType() const override;

    //! @brief ... check parameters of the constitutive relationship
    void CheckParameters() const override;

    //! @brief ... print information about the object
    //! @param rVerboseLevel ... verbosity of the information
    //! @param rLogger stream for the output
    void Info(unsigned short rVerboseLevel, Logger& rLogger) const override;

    //! @brief ... returns true, if a material model has tmp static data (which has to be updated before stress or
    //! stiffness are calculated)
    //! @return ... see brief explanation
    bool HaveTmpStaticData() const override
    {
        return false;
    }


protected:
    //! @brief ... Stiffness tensor components in Voigt notation
    Eigen::Matrix<double, 6, 6> mStiffness;

    //! @brief ... density \f$ \rho \f$
    double mRho;
};
}

//#ifdef ENABLE_SERIALIZATION
// BOOST_CLASS_EXPORT_KEY(NuTo::LinearElasticAnisotropic)
//#endif //ENABLE_SERIALIZATION
