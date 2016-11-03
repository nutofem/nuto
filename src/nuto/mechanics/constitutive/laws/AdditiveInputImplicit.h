#pragma once

#include "nuto/mechanics/constitutive/laws/AdditiveBase.h"
#include "nuto/mechanics/constitutive/staticData/IPConstitutiveLaw.h"

#include "nuto/mechanics/constitutive/staticData/DataMoistureTransport.h" // TODO ----> WEG!

#include <set>
#include <vector>

namespace NuTo
{

class AdditiveInputImplicit : public AdditiveBase
{
public:

    typedef Constitutive::StaticData::DataMoistureTransport StaticDataType;
    using Data = typename Constitutive::StaticData::DataContainer<StaticDataType>;

    //! @brief constructor
    AdditiveInputImplicit() : AdditiveBase()
    {
        //VHIRTHAMTODO ---> Get number time derivatives during construction (as parameter)
        mComputableDofCombinations.resize(2);
    }

//    // has no ip static data itself
//    std::unique_ptr<Constitutive::IPConstitutiveLawBase> CreateIPLaw()
//    {
//        return std::make_unique<Constitutive::IPConstitutiveLawWithoutData<AdditiveInputImplicit>>(*this);
//    }

    //! @brief creates corresponding IPConstitutiveLaw
    std::unique_ptr<Constitutive::IPConstitutiveLawBase> CreateIPLaw() override
    {
        return std::make_unique<Constitutive::IPConstitutiveLaw<AdditiveInputImplicit>>(*this, StaticDataType());
    }

    //! @brief Evaluate the constitutive relation.
    //! @param rConstitutiveInput Input to the constitutive law (strain, temp gradient etc.).
    //! @param rConstitutiveOutput Output to the constitutive law (stress, stiffness, heat flux etc.).
    //! @param staticData Pointer to the static data.
    template <int TDim>
    NuTo::eError Evaluate(
        const ConstitutiveInputMap& rConstitutiveInput,
        const ConstitutiveOutputMap& rConstitutiveOutput,
        Data& rStaticData);


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

    std::vector<std::set<std::pair<Node::eDof,Node::eDof>>> mComputableDofCombinations;

    //! @brief debug variable to avoid that a constitutive law can be attached after allocation of static data.
    mutable bool mStaticDataAllocated = false;
};

}
