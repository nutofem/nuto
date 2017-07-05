//============================================================================
// Name        : FibreMatrixBondStressSlip.cpp
// Author      : Philip Huschke
// Version     : 26 Aug 2015
// Copyright   :
// Description : Constitutive law for the interface between fibre and matrix
//============================================================================

#pragma once

#include "mechanics/constitutive/staticData/IPConstitutiveLaw.h"
#include "mechanics/constitutive/ConstitutiveBase.h"

namespace NuTo
{
template <int TRows, int TCols> class ConstitutiveMatrix;
template <int TRows> class ConstitutiveVector;
class ConstitutiveScalar;
class ConstitutiveStaticDataBondStressSlip;

class FibreMatrixBondStressSlip: public ConstitutiveBase
{

public:
    typedef double StaticDataType;
    using Data = typename Constitutive::StaticData::DataContainer<double>;

    std::unique_ptr<Constitutive::IPConstitutiveLawBase> CreateIPLaw() override
    {
        return std::make_unique<Constitutive::IPConstitutiveLaw<FibreMatrixBondStressSlip>>(*this, 0.);
    }

    //! @brief Create constitutive law for fibre-matrix interface.
    //! @param dimension Global dimension of the structure the fibres are embedded in.
    FibreMatrixBondStressSlip(int dimension);

    //! @brief ... determines the constitutive inputs needed to evaluate the constitutive outputs
    //! @param rConstitutiveOutput ... desired constitutive outputs
    //! @param rInterpolationType ... interpolation type to determine additional inputs
    //! @return constitutive inputs needed for the evaluation
    ConstitutiveInputMap GetConstitutiveInputs(
        const ConstitutiveOutputMap& rConstitutiveOutput,
        const InterpolationType& rInterpolationType) const override;

    //! @brief Evaluate the constitutive relation.
    //! @param rConstitutiveInput Input to the constitutive law (strain, temp gradient etc.).
    //! @param rConstitutiveOutput Output to the constitutive law (stress, stiffness, heat flux etc.).
    //! @param rStaticData static data container.
    template <int TDim>
    void Evaluate(
            const ConstitutiveInputMap& rConstitutiveInput,
            const ConstitutiveOutputMap& rConstitutiveOutput,
            Data& rStaticData);

    //! @brief Calculates the current static data based on the given CALCULATE_STATIC_DATA input.
    //! @param rStaticData Static data passed to the law.
    //! @param rConstitutiveInput Input to the constitutive law (strain, temp gradient etc.).
    //! @return Updated current static data.
    double GetCurrentStaticData(Data& rStaticData, const ConstitutiveInputMap& rConstitutiveInput) const;

    //! @brief ... gets a variable of the constitutive law which is selected by an enum
    //! @param rIdentifier ... Enum to identify the requested variable
    //! @return ... value of the requested variable
    virtual double GetParameterDouble(Constitutive::eConstitutiveParameter rIdentifier) const override;

    //! @brief ... sets a variable of the constitutive law which is selected by an enum
    //! @param rIdentifier ... Enum to identify the requested variable
    //! @param rValue ... new value for requested variable
    virtual void SetParameterDouble(Constitutive::eConstitutiveParameter rIdentifier, double rValue) override;

    //! @brief ... get type of constitutive relationship
    //! @return ... type of constitutive relationship
    //! @sa eConstitutiveType
    Constitutive::eConstitutiveType GetType() const override;

    //! @brief ... check parameters of the constitutive relationship
    void CheckParameters() const override;

    //! @brief ... determines which submatrices of a multi-doftype problem can be solved by the constitutive law
    //! @param rDofRow ... row dof
    //! @param rDofCol ... column dof
    //! @param rTimeDerivative ... time derivative
    virtual bool CheckDofCombinationComputable(Node::eDof rDofRow,
                                                Node::eDof rDofCol,
                                                int rTimeDerivative) const override;

    //! @brief ... print information about the object
    //! @param rVerboseLevel ... verbosity of the information
    //! @param rLogger stream for the output
    void Info(unsigned short rVerboseLevel, Logger& rLogger) const override;

    //! @brief ... returns true, if a material model has tmp static data (which has to be updated before stress or stiffness are calculated)
    //! @return ... see brief explanation
    bool HaveTmpStaticData() const override
    {
        return false;
    }

private:

    //! @brief ... maximum bond stress
    double mMaxBondStress;

    //! @brief ... residual bond stress
    double mResidualBondStress;

    //! @brief ... slip at maximum bond stress
    double mSlipAtMaxBondStress;

    //! @brief ... slip at residual bond stress (bond stress stays constant for greater slips)
    double mSlipAtResidualBondStress;

    //! @brief ... parameter that controls the pre peak shape of the bond stress-slip curve
    double mAlpha;

    //! @brief ... normal penalty stiffness to avoid penetration
    double mNormalStiffness;

    //! @brief Global dimension of the structure the fibres are embedded in.
    unsigned int mGlobalDimension;
};

}

