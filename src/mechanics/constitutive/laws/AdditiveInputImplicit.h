#pragma once

#include "mechanics/constitutive/laws/AdditiveBase.h"

#include "mechanics/constitutive/staticData/IPAdditiveInputImplicit.h"
#include "mechanics/constitutive/staticData/DataAdditiveInputImplicit.h"


#include <set>
#include <vector>

namespace NuTo
{

class AdditiveInputImplicit : public AdditiveBase
{
public:
    typedef Constitutive::StaticData::DataAdditiveInputImplicit StaticDataType;
    using Data = typename Constitutive::StaticData::DataContainer<StaticDataType>;

    //! @brief constructor
    AdditiveInputImplicit(const int& rNumTimeDerivatives)
        : AdditiveBase(rNumTimeDerivatives)
    {
    }

    //! @brief creates corresponding IPConstitutiveLaw
    std::unique_ptr<Constitutive::IPConstitutiveLawBase> CreateIPLaw() override
    {
        return std::make_unique<Constitutive::IPAdditiveInputImplicit>(*this, StaticDataType());
    }

    //! @brief Evaluate the constitutive relation.
    //! @param rConstitutiveInput Input to the constitutive law (strain, temp gradient etc.).
    //! @param rConstitutiveOutput Output to the constitutive law (stress, stiffness, heat flux etc.).
    //! @param staticData Pointer to the static data.
    template <int TDim>
    void Evaluate(const ConstitutiveInputMap&, const ConstitutiveOutputMap&)
    {
        throw MechanicsException(
                __PRETTY_FUNCTION__,
                "Additive Law cannot be evaluated. Their IPAdditiveInputExplicit should be evaluated instead.");
    }


    //! @brief Determines the constitutive inputs needed to evaluate the constitutive outputs.
    //! @param rConstitutiveOutput Desired constitutive outputs.
    //! @param rInterpolationType Interpolation type to determine additional inputs.
    //! @return Constitutive inputs needed for the evaluation.
    virtual ConstitutiveInputMap GetConstitutiveInputs(const ConstitutiveOutputMap& rConstitutiveOutput,
                                                       const InterpolationType& rInterpolationType) const override;

    //! @brief Get the type of the constitutive relationship.
    //! @return Type of the constitutive relationship.
    //! @see eConstitutiveType
    virtual Constitutive::eConstitutiveType GetType() const override;


private:
    //! @brief Adds all calculable dof combinations of an attached constitutive law to an internal storage
    //! @param rConstitutiveLaw ... constitutive law
    void AddCalculableDofCombinations(NuTo::ConstitutiveBase& rConstitutiveLaw);

    std::vector<std::set<std::pair<Node::eDof, Node::eDof>>> mComputableDofCombinations;
};
}
