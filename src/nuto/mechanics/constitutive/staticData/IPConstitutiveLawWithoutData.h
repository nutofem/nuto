#pragma once

#include "nuto/mechanics/constitutive/staticData/IPConstitutiveLawBase.h"
#include <type_traits>

namespace NuTo
{
namespace Constitutive
{

template<typename TLaw>
class IPConstitutiveLawWithoutData : public IPConstitutiveLawBase
{
public:

    static_assert(std::is_base_of<ConstitutiveBase, TLaw>::value,"TLaw must be derived from NuTo::ConstitutiveBase");

    //! @brief constructor
    //! @param rLaw underlying constitutive law
    IPConstitutiveLawWithoutData(TLaw& rLaw) : mLaw(rLaw) {}

    ConstitutiveBase& GetConstitutiveLaw() const
    {
        return static_cast<ConstitutiveBase&>(mLaw);
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
