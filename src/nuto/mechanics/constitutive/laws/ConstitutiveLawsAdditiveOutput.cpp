#include "ConstitutiveLawsAdditiveOutput.h"










//VHIRTHAMTODO --- maybe store possible dof combinations one time in a set when constitutive law is added instead of looping over each law every time
bool NuTo::ConstitutiveLawsAdditiveOutput::CheckDofCombinationComputeable(NuTo::Node::eDof rDofRow, NuTo::Node::eDof rDofCol, int rTimeDerivative) const
{
    if(mComputableDofCombinations[rTimeDerivative].find(std::pair<Node::eDof,Node::eDof>(rDofRow,rDofCol)) != mComputableDofCombinations[rTimeDerivative].end())
        return true;
    return false;
}

NuTo::Error::eError NuTo::ConstitutiveLawsAdditiveOutput::Evaluate1D(NuTo::ElementBase *rElement,
                                                                   int rIp,
                                                                   const NuTo::ConstitutiveInputMap &rConstitutiveInput,
                                                                   const NuTo::ConstitutiveOutputMap &rConstitutiveOutput)
{
    Error::eError error = Error::SUCCESSFUL;
    try
    {
        for(unsigned int i=0; i<mConstitutiveLaws.size(); ++i)
        {
            error = mConstitutiveLaws[i]->Evaluate1D(rElement,
                                                     rIp,
                                                     rConstitutiveInput,
                                                     rConstitutiveOutput);
            if(error != Error::SUCCESSFUL)
            {
                throw MechanicsException(__PRETTY_FUNCTION__,"Attached constitutive law returned an error code. Can't handle this");
            }
        }
    }
    catch(Exception e)
    {
        e.AddMessage(__PRETTY_FUNCTION__,"Exception while evaluating constitutive law attached to an additive link.");
        throw e;
    }
    return error;
}











NuTo::Error::eError NuTo::ConstitutiveLawsAdditiveOutput::Evaluate2D(NuTo::ElementBase *rElement,
                                                                   int rIp,
                                                                   const NuTo::ConstitutiveInputMap &rConstitutiveInput,
                                                                   const NuTo::ConstitutiveOutputMap &rConstitutiveOutput)
{
    Error::eError error = Error::SUCCESSFUL;
    try
    {
        for(unsigned int i=0; i<mConstitutiveLaws.size(); ++i)
        {
            error = mConstitutiveLaws[i]->Evaluate2D(rElement,
                                                     rIp,
                                                     rConstitutiveInput,
                                                     rConstitutiveOutput);
            if(error != Error::SUCCESSFUL)
            {
                throw MechanicsException(__PRETTY_FUNCTION__,"Attached constitutive law returned an error code. Can't handle this");
            }
        }
    }
    catch(Exception e)
    {
        e.AddMessage(__PRETTY_FUNCTION__,"Exception while evaluating constitutive law attached to an additive link.");
        throw e;
    }
    return error;
}










NuTo::Error::eError NuTo::ConstitutiveLawsAdditiveOutput::Evaluate3D(NuTo::ElementBase *rElement,
                                                                   int rIp,
                                                                   const NuTo::ConstitutiveInputMap &rConstitutiveInput,
                                                                   const NuTo::ConstitutiveOutputMap &rConstitutiveOutput)
{
    Error::eError error = Error::SUCCESSFUL;
    try
    {
        for(unsigned int i=0; i<mConstitutiveLaws.size(); ++i)
        {
            error = mConstitutiveLaws[i]->Evaluate3D(rElement,
                                                     rIp,
                                                     rConstitutiveInput,
                                                     rConstitutiveOutput);
            if(error != Error::SUCCESSFUL)
            {
                throw MechanicsException(__PRETTY_FUNCTION__,"Attached constitutive law returned an error code. Can't handle this");
            }
        }
    }
    catch(Exception e)
    {
        e.AddMessage(__PRETTY_FUNCTION__,"Exception while evaluating constitutive law attached to an additive link.");
        throw e;
    }
    return error;
}




NuTo::ConstitutiveInputMap NuTo::ConstitutiveLawsAdditiveOutput::GetConstitutiveInputs(const NuTo::ConstitutiveOutputMap &rConstitutiveOutput,
                                                                                     const NuTo::InterpolationType &rInterpolationType) const
{
    ConstitutiveInputMap constitutiveInputMap;


    for(unsigned int i=0; i<mConstitutiveLaws.size(); ++i)
    {

        ConstitutiveInputMap singleLawInputMap = mConstitutiveLaws[i]->GetConstitutiveInputs(rConstitutiveOutput,
                                                                                             rInterpolationType);

        constitutiveInputMap.insert(singleLawInputMap.begin(),singleLawInputMap.end());
    }

    return constitutiveInputMap;
}

void NuTo::ConstitutiveLawsAdditiveOutput::AddCalculableDofCombinations(NuTo::ConstitutiveBase *rConstitutiveLaw)
{
    std::set<Node::eDof> allDofs = Node::GetDofSet();
    for (unsigned int i=0; i<mComputableDofCombinations.size(); ++i)
    for (auto itRow : allDofs)
        for (auto itCol : allDofs)
        {
            if (rConstitutiveLaw->CheckDofCombinationComputeable(itRow,itCol,i))
                    mComputableDofCombinations[i].emplace(itRow,itCol);
        }

}
