#pragma once

#include "nuto/mechanics/constitutive/staticData/IPConstitutiveLawBase.h"

namespace NuTo
{
namespace Constitutive
{

template<typename TLaw>
class IPConstitutiveLawWithoutData : public IPConstitutiveLawBase
{
public:

    //! @brief constructor
    //! @param rLaw underlying constitutive law
    IPConstitutiveLawWithoutData(TLaw& rLaw) : mLaw(rLaw) {}

    //! @brief default copy constructor
    IPConstitutiveLawWithoutData(const IPConstitutiveLawWithoutData&) = default;

    //! @brief default move constructor
    IPConstitutiveLawWithoutData(      IPConstitutiveLawWithoutData&&) = default;

    //! @brief default destuctor
    ~IPConstitutiveLawWithoutData() = default;

    TLaw& GetConstitutiveLaw() const
    {
        return mLaw;
    }

protected:
    //! @brief Evaluate the constitutive relation in 1D
    //! @param rConstitutiveInput Input to the constitutive law
    //! @param rConstitutiveOutput Output of the constitutive law
    eError Evaluate1D(const ConstitutiveInputMap& rConstitutiveInput,
                      const ConstitutiveOutputMap& rConstitutiveOutput) override
    {
        return mLaw.template Evaluate<1>(rConstitutiveInput, rConstitutiveOutput);
    }

    //! @brief Evaluate the constitutive relation in 1D
    //! @param rConstitutiveInput Input to the constitutive law
    //! @param rConstitutiveOutput Output of the constitutive law
    eError Evaluate2D(const ConstitutiveInputMap& rConstitutiveInput,
                      const ConstitutiveOutputMap& rConstitutiveOutput) override
    {
        return mLaw.template Evaluate<2>(rConstitutiveInput, rConstitutiveOutput);
    }

    //! @brief Evaluate the constitutive relation in 1D
    //! @param rConstitutiveInput Input to the constitutive law
    //! @param rConstitutiveOutput Output of the constitutive law
    eError Evaluate3D(const ConstitutiveInputMap& rConstitutiveInput,
                      const ConstitutiveOutputMap& rConstitutiveOutput) override
    {
        return mLaw.template Evaluate<3>(rConstitutiveInput, rConstitutiveOutput);
    }

private:

    TLaw& mLaw;
};


} // namespace Constitutive
} // namespace NuTo
