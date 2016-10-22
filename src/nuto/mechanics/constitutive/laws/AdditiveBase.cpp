#include "nuto/mechanics/constitutive/laws/AdditiveBase.h"
#include "nuto/mechanics/nodes/NodeEnum.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOBase.h"



NuTo::AdditiveBase::AdditiveBase()
{
    mComputableDofCombinations.resize(3);
}

NuTo::AdditiveBase::AdditiveBase(const AdditiveBase& rOther)
{
    *this = rOther;
}

NuTo::AdditiveBase::AdditiveBase(AdditiveBase&& rOther)
{
    *this = std::move(rOther);
}

NuTo::AdditiveBase& NuTo::AdditiveBase::operator =(const AdditiveBase& rOther)
{
    mComputableDofCombinations = rOther.mComputableDofCombinations;
    mStaticDataAllocated = rOther.mStaticDataAllocated;
    mSublaws.reserve(rOther.mSublaws.size());
    for (auto& ipLaw : rOther.mSublaws)
        mSublaws.push_back(ipLaw->GetConstitutiveLaw().CreateIPLaw());
    return *this;
}
NuTo::AdditiveBase& NuTo::AdditiveBase::operator =(AdditiveBase&& rOther)
{
    mComputableDofCombinations = std::move(rOther.mComputableDofCombinations);
    mStaticDataAllocated = rOther.mStaticDataAllocated; // just a bool
    mSublaws.reserve(rOther.mSublaws.size());
    for (auto& ipLaw : rOther.mSublaws)
        mSublaws.push_back(std::move(ipLaw));
    return *this;
}


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

    mSublaws.push_back(rConstitutiveLaw.CreateIPLaw());
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
    for (auto& sublaw : mSublaws)
    {
        if(!sublaw->GetConstitutiveLaw().CheckElementCompatibility(rElementType))
            return false;
    }
    return true;
}


void NuTo::AdditiveBase::CheckParameters() const
{
    for (auto& sublaw : mSublaws)
    {
        sublaw->GetConstitutiveLaw().CheckParameters();
    }
}


bool NuTo::AdditiveBase::HaveTmpStaticData() const 
{
    for (auto& sublaw : mSublaws)
    {
        if(sublaw->GetConstitutiveLaw().HaveTmpStaticData())
            return true;
    }
    return false;
}

bool NuTo::AdditiveBase::CheckDofCombinationComputable(NuTo::Node::eDof rDofRow, NuTo::Node::eDof rDofCol, 
        int rTimeDerivative) const
{
    return mComputableDofCombinations[rTimeDerivative].find(std::pair<Node::eDof, Node::eDof>(rDofRow, rDofCol))
        != mComputableDofCombinations[rTimeDerivative].end();
}


NuTo::ConstitutiveInputMap NuTo::AdditiveBase::GetConstitutiveInputs(
        const NuTo::ConstitutiveOutputMap &rConstitutiveOutput, const NuTo::InterpolationType &rInterpolationType) const
{
    ConstitutiveInputMap constitutiveInputMap;

    for (auto& sublaw : mSublaws)
    {
        ConstitutiveInputMap singleLawInputMap = sublaw->GetConstitutiveLaw().GetConstitutiveInputs(rConstitutiveOutput, rInterpolationType);
        constitutiveInputMap.Merge(singleLawInputMap);
    }
    return constitutiveInputMap;
}
