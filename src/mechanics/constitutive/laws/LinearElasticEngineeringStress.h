#pragma once


#include "mechanics/constitutive/ConstitutiveBase.h"


namespace NuTo
{
namespace Constitutive
{
class IPConstitutiveLawBase;
} // namespace Constitutive

//! @brief ... linear elastic material model
/*!
 * Assuming linear elastic material behavior, the one-dimensional constitutive relationship reads
 * \f{align*}{
 *   \sigma_{xx} = E \varepsilon_{xx}
 * \f}
 * where \f$ E \f$ is the Young's modulus, \f$ \sigma_{xx} \f$ is the Cauchy stress
 * and \f$ \varepsilon_{xx} \f$ is the Engineering strain vector.
 * It is to be noted, that this relationship holds for Green strains and the second Piola-Kirchhoff stress as well.
 * Assuming linear elastic isotropic material behavior, the general three-dimensional constitutive relationship reads
 * \f{align*}{
 *  \begin{bmatrix}
 *    \sigma_{xx}\\
 *    \sigma_{yy}\\
 *    \sigma_{zz}\\
 *    \sigma_{xy}\\
 *    \sigma_{yz}\\
 *    \sigma_{zx}
 *  \end{bmatrix} = \dfrac{E (1 - \nu)}{(1+\nu)(1-2\nu)} \begin{bmatrix}
 *    1 & \dfrac{\nu}{1-\nu} & \dfrac{\nu}{1-\nu} & 0 & 0 & 0\\
 *    \dfrac{\nu}{1-\nu} & 1 & \dfrac{\nu}{1-\nu} & 0 & 0 & 0\\
 *    \dfrac{\nu}{1-\nu} & \dfrac{\nu}{1-\nu} & 1 & 0 & 0 & 0\\
 *    0 & 0 & 0 & \dfrac{1-2\nu}{2(1-\nu)} & 0 & 0\\
 *    0 & 0 & 0 & 0 & \dfrac{1-2\nu}{2(1-\nu)} & 0\\
 *    0 & 0 & 0 & 0 & 0 & \dfrac{1-2\nu}{2(1-\nu)}
 *  \end{bmatrix}
 *  \begin{bmatrix}
 *    \varepsilon_{xx}\\
 *    \varepsilon_{yy}\\
 *    \varepsilon_{zz}\\
 *    \gamma_{xy}\\
 *    \gamma_{yz}\\
 *    \gamma_{zx}
 *  \end{bmatrix},
 * \f}
 * two-dimensional plain stress constitutive relationship reads
 * \f{align*}{
 *  \begin{bmatrix}
 *    \sigma_{xx}\\
 *    \sigma_{yy}\\
 *    \sigma_{xy}
 *  \end{bmatrix} = \dfrac{E}{(1+\nu^2)} \begin{bmatrix}
 *    1 & \nu & 0\\
 *    \nu & 1 & 0\\
 *    0 & 0 \dfrac{1-\nu}{2}
 *  \end{bmatrix}
 *  \begin{bmatrix}
 *    \varepsilon_{xx}\\
 *    \varepsilon_{yy}\\
 *    \gamma_{xy}
 *  \end{bmatrix},
 * \f}
 * where \f$ E \f$ is the Young's modulus, \f$ \nu \f$ is the Poisson's ratio,
 * \f$ \boldsymbol{\sigma} \f$ are the components of the Cauchy stress vector,
 * and \f$ \boldsymbol{\varepsilon} \f$ are the components of the Engineering strain vector.
 * It is to be noted, that this relationship holds for Green strains and the second Piola-Kirchhoff stress as well.
 */
//! @author Jörg F. Unger, ISM
//! @date July 2012
class LinearElasticEngineeringStress : public ConstitutiveBase
{

public:
    LinearElasticEngineeringStress();

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
    //! @brief ... Young's modulus \f$ E \f$
    double mE;

    //! @brief ... Poisson's ratio \f$ \nu \f$
    double mNu;

    //! @brief ... density \f$ \rho \f$
    double mRho;
};
}
