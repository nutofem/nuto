#include "nuto/mechanics/constitutive/laws/AdditiveBase.h"
#include "nuto/mechanics/nodes/NodeEnum.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOBase.h"


void NuTo::AdditiveBase::AddConstitutiveLaw(NuTo::ConstitutiveBase& rConstitutiveLaw,
        Constitutive::eInput rModiesInput)
{
    if(rConstitutiveLaw.HaveTmpStaticData())
        throw MechanicsException(__PRETTY_FUNCTION__,
            "Constitutive law has tmp static data! The HaveTmpStaticData function is only called on construction of "
            "the AdditiveInputExplicit law, but at this time, no constitutive law is attached. Therefore it does not "
            "know if it will have tmpstatic data or not and returns false by default. Find a way to update this "
            "information at the necessary code sections if a law with tmpstatic data is attached.");

    if(mStaticDataAllocated)
        throw MechanicsException(__PRETTY_FUNCTION__,
                "All constitutive laws have to be attached before static data is allocated!");

    mSublaws.push_back(&rConstitutiveLaw);
    AddCalculableDofCombinations(rConstitutiveLaw);
}


void NuTo::AdditiveBase::AddCalculableDofCombinations(NuTo::ConstitutiveBase& rConstitutiveLaw)
{
    std::set<Node::eDof> allDofs = Node::GetDofSet();
    for (unsigned int i = 0;  i < mComputableDofCombinations.size(); ++i)
        for (auto itRow : allDofs)
            for (auto itCol : allDofs)
            {
                if (rConstitutiveLaw.CheckDofCombinationComputable(itRow,itCol,i))
                        mComputableDofCombinations[i].emplace(itRow,itCol);
            }
}


bool NuTo::AdditiveBase::CheckElementCompatibility(Element::eElementType rElementType) const
{
    for (auto sublaw : mSublaws)
    {
        if(!sublaw->CheckElementCompatibility(rElementType))
            return false;
    }
    return true;
}


void NuTo::AdditiveBase::CheckParameters() const
{
    for (auto sublaw : mSublaws)
    {
        sublaw->CheckParameters();
    }
}


bool NuTo::AdditiveBase::HaveTmpStaticData() const 
{
    for (auto sublaw : mSublaws)
    {
        if(sublaw->HaveTmpStaticData())
            return true;
    }
    return false;
}

bool NuTo::AdditiveBase::CheckDofCombinationComputable(NuTo::Node::eDof rDofRow, NuTo::Node::eDof rDofCol, 
        int rTimeDerivative) const
{
    if (mComputableDofCombinations[rTimeDerivative].find(std::pair<Node::eDof,Node::eDof>(rDofRow,rDofCol))
            != mComputableDofCombinations[rTimeDerivative].end())
        return true;
    else
        return false;
}


NuTo::ConstitutiveInputMap NuTo::AdditiveBase::GetConstitutiveInputs(
        const NuTo::ConstitutiveOutputMap &rConstitutiveOutput, const NuTo::InterpolationType &rInterpolationType) const
{
    ConstitutiveInputMap constitutiveInputMap;

    for (auto sublaw : mSublaws)
    {
        ConstitutiveInputMap singleLawInputMap = sublaw->GetConstitutiveInputs(rConstitutiveOutput, rInterpolationType);
        constitutiveInputMap.Merge(singleLawInputMap);
    }
    return constitutiveInputMap;
}

template NuTo::Constitutive::StaticData::Component* NuTo::AdditiveBase::AllocateStaticData<1>(const NuTo::ElementBase *rElement) const;
template NuTo::Constitutive::StaticData::Component* NuTo::AdditiveBase::AllocateStaticData<2>(const NuTo::ElementBase *rElement) const;
template NuTo::Constitutive::StaticData::Component* NuTo::AdditiveBase::AllocateStaticData<3>(const NuTo::ElementBase *rElement) const;
