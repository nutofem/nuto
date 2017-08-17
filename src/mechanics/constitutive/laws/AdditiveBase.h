#pragma once

#include <set>
#include <vector>

#include "mechanics/constitutive/ConstitutiveBase.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"

namespace NuTo
{
namespace Constitutive
{
class IPAdditiveInputExplicit;
class IPAdditiveInputImplicit;
class IPAdditiveOutput;
}
class AdditiveBase : public ConstitutiveBase
{
public:
    friend class Constitutive::IPAdditiveInputExplicit;
    friend class Constitutive::IPAdditiveInputImplicit;
    friend class Constitutive::IPAdditiveOutput;

    //! @brief ctor
    AdditiveBase(const int& rNumTimeDerivatives);

    virtual ~AdditiveBase() = default;
    AdditiveBase(const AdditiveBase& rOther) = default;
    AdditiveBase(AdditiveBase&& rOther) = default;

    AdditiveBase& operator=(const AdditiveBase& rOther) = default;
    AdditiveBase& operator=(AdditiveBase&& rOther) = default;


    //! @brief Adds a constitutive law to a model that combines multiple constitutive laws (additive, parallel)
    //! @param rConstitutiveLaw Constitutive law to be added.
    //! @param rModiesInput Enum which defines which input is modified by a constitutive law.
    virtual void AddConstitutiveLaw(NuTo::ConstitutiveBase& rConstitutiveLaw,
                                    Constitutive::eInput rModiesInput = Constitutive::eInput::NONE);

    //! @brief Check parameters of the constitutive relationship.
    virtual void CheckParameters() const override;

    //! @brief Checks whether a material model has tmp static data, which has to be updated before stress or
    //! stiffness are calculated.
    virtual bool HaveTmpStaticData() const override;

    //! @brief returns the sublaw with index rInted
    //! @param rIndex ... index
    //! @return reference to ip law
    ConstitutiveBase& GetSublaw(int rIndex);

    virtual ConstitutiveInputMap GetConstitutiveInputs(const ConstitutiveOutputMap& rConstitutiveOutput) const override;

    virtual bool CheckDofCombinationComputable(Node::eDof rDofRow, Node::eDof rDofCol,
                                               int rTimeDerivative) const override;

protected:
    //! @brief Adds all calculable DOF combinations of an attached constitutive law to an internal storage.
    //! @param rConstitutiveLaw Constitutive law whose DOF combinations are added to additive law.
    void AddCalculableDofCombinations(NuTo::ConstitutiveBase& rConstitutiveLaw);

    //! @brief Debug variable to avoid that a constitutive law can be attached after allocation of static data.
    mutable bool mStaticDataAllocated = false;

    //! @brief Reference to the variable that stores the number of time derivatives (Original variable should be stored
    //! at the structure)
    const int& mNumTimeDerivatives;


    //! @brief Vector storing the IPConstitutiveBase sublaws.
    std::vector<NuTo::ConstitutiveBase*> mSublaws;

    //! @brief Vector of all the computable DOF combinations.
    std::vector<std::set<std::pair<Node::eDof, Node::eDof>>> mComputableDofCombinations;
};
}
