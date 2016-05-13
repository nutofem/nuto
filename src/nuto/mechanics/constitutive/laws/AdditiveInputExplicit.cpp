#include "AdditiveInputExplicit.h"
#include "nuto/mechanics/constitutive/inputoutput/EngineeringStrain.h"
#include "nuto/mechanics/constitutive/inputoutput/EngineeringStress.h"









void NuTo::AdditiveInputExplicit::AddConstitutiveLaw(NuTo::ConstitutiveBase *rConstitutiveLaw, Constitutive::Input::eInput rModiesInput)
{
    if(mStaticDataAllocated)
        throw MechanicsException(__PRETTY_FUNCTION__,"All constitutive laws have to be attached before static data is allocated!");
    if(rModiesInput == Constitutive::Input::NONE)
    {
        if (mConstitutiveLawOutput != nullptr)
            throw MechanicsException(__PRETTY_FUNCTION__,std::string("There can be only one!!! --- This additive input law only accepts one law which calculates the output. All other laws ")+
                                     " are only allowed to modify the input to this law. Specify the modifying laws by providing the enum of the modified input as second function parameter.");
        mConstitutiveLawOutput = rConstitutiveLaw;
    }
    else
    {
        mModifiedInputs.insert(rModiesInput);
        mConstitutiveLawsModInput.push_back({rConstitutiveLaw,rModiesInput});
    }
    AddCalculableDofCombinations(rConstitutiveLaw);
}



bool NuTo::AdditiveInputExplicit::CheckDofCombinationComputable(NuTo::Node::eDof rDofRow, NuTo::Node::eDof rDofCol, int rTimeDerivative) const
{
    if(mComputableDofCombinations[rTimeDerivative].find(std::pair<Node::eDof,Node::eDof>(rDofRow,rDofCol)) != mComputableDofCombinations[rTimeDerivative].end())
        return true;
    return false;
}

template <int TDim>
NuTo::Error::eError NuTo::AdditiveInputExplicit::EvaluateAdditiveInputExplicit(NuTo::ElementBase *rElement,
                                                                               int rIp,
                                                                               const NuTo::ConstitutiveInputMap &rConstitutiveInput,
                                                                               const NuTo::ConstitutiveOutputMap &rConstitutiveOutput)
{
    static constexpr int VoigtDim = ConstitutiveIOBase::GetVoigtDim(TDim);
    Error::eError error = Error::SUCCESSFUL;

     NuTo::ConstitutiveInputMap copiedInputMap;

    // Need deep copy of all modified inputs, otherwise the modifications will possibly effect other constitutive laws that are coupled in a additive input or output law
    // So far only the engineering strain input is modified. If there are more modified Inputs in the future, a more general allocation method should be chosen
    EngineeringStrain<TDim> engineeringStrain;
    ConstitutiveVector<ConstitutiveIOBase::GetVoigtDim(TDim)> d_EngineeringStrain_D_RH;
    ConstitutiveVector<ConstitutiveIOBase::GetVoigtDim(TDim)> d_EngineeringStrain_D_WV;

    engineeringStrain.setZero();
    d_EngineeringStrain_D_RH.setZero();
    d_EngineeringStrain_D_WV.setZero();

    for (auto itInput : rConstitutiveInput)
    {
        if(mModifiedInputs.find(itInput.first) != mModifiedInputs.end())
        {
            //deep copy
            switch(itInput.first)
            {
            case Constitutive::Input::ENGINEERING_STRAIN:
                switch(TDim)
                {
                case 1:
                    engineeringStrain.AsEngineeringStrain1D() -= itInput.second->AsEngineeringStrain1D();
                    break;
                case 2:
                    engineeringStrain.AsEngineeringStrain2D() -= itInput.second->AsEngineeringStrain2D();
                    break;
                case 3:
                    engineeringStrain.AsEngineeringStrain3D() -= itInput.second->AsEngineeringStrain3D();
                    break;

                default:
                    throw MechanicsException(__PRETTY_FUNCTION__,"invalid dimension");
                }
                copiedInputMap[Constitutive::Input::ENGINEERING_STRAIN] = &engineeringStrain;
                break;

            default:
                throw MechanicsException(__PRETTY_FUNCTION__,std::string("No method provided to create a deep copy of input with enum: ")+Constitutive::InputToString(itInput.first));
            }
        }
        else
        {
            //shallow copy
            copiedInputMap[itInput.first] = itInput.second;
        }
    }





    for(unsigned int i=0; i<mConstitutiveLawsModInput.size(); ++i)
    {

        // Generate modified output map for constitutive law
        NuTo::ConstitutiveOutputMap modifiedOutputMap = rConstitutiveOutput;

        if(copiedInputMap.find(Constitutive::Input::ENGINEERING_STRAIN) != copiedInputMap.end())
        {
                modifiedOutputMap[Constitutive::Output::ENGINEERING_STRAIN] = &engineeringStrain;
        }
        if(modifiedOutputMap.find(Constitutive::Output::D_ENGINEERING_STRESS_D_RELATIVE_HUMIDITY) != modifiedOutputMap.end() &&
                                  mConstitutiveLawsModInput[i].second == Constitutive::Input::ENGINEERING_STRAIN)
        {
            modifiedOutputMap[Constitutive::Output::D_ENGINEERING_STRAIN_D_RELATIVE_HUMIDITY] = &d_EngineeringStrain_D_RH;
        }
        if(modifiedOutputMap.find(Constitutive::Output::D_ENGINEERING_STRESS_D_WATER_VOLUME_FRACTION) != modifiedOutputMap.end() &&
                                  mConstitutiveLawsModInput[i].second == Constitutive::Input::ENGINEERING_STRAIN)
        {
            modifiedOutputMap[Constitutive::Output::D_ENGINEERING_STRAIN_D_WATER_VOLUME_FRACTION] = &d_EngineeringStrain_D_WV;
        }


        // evaluate constitutive law
        switch(TDim)
        {
        case 1:
            error = mConstitutiveLawsModInput[i].first->Evaluate1D(rElement,
                                                                   rIp,
                                                                   rConstitutiveInput,
                                                                   modifiedOutputMap);
            break;

        case 2:
            error = mConstitutiveLawsModInput[i].first->Evaluate2D(rElement,
                                                                   rIp,
                                                                   rConstitutiveInput,
                                                                   modifiedOutputMap);
            break;

        case 3:
            error = mConstitutiveLawsModInput[i].first->Evaluate3D(rElement,
                                                                   rIp,
                                                                   rConstitutiveInput,
                                                                   modifiedOutputMap);
            break;

        default:
            throw MechanicsException(__PRETTY_FUNCTION__,"invalid dimension");
        }
        if(error != Error::SUCCESSFUL)
        {
            throw MechanicsException(__PRETTY_FUNCTION__,"Attached constitutive law returned an error code. Can't handle this");
        }
//        for(auto itOutput : modifiedOutputMap)
//        {
//            if(itOutput.second->GetIsCalculated() == false) // The output list has to be created differently to handle this without exceptions (for example: Thermal has no D_ENGINEERING_STRAIN_D_RELATIVE_HUMIDITY)
//                throw MechanicsException(__PRETTY_FUNCTION__,std::string("Attached constitutive law has not calculated the requested output ") + Constitutive::OutputToString(itOutput.first));

//        }
    }
    engineeringStrain.AsVector() = engineeringStrain.AsVector() * -1;

    switch(TDim)
    {
    case 1:
        error = mConstitutiveLawOutput->Evaluate1D(rElement,
                                                   rIp,
                                                   copiedInputMap,
                                                   rConstitutiveOutput);
        break;

    case 2:
        error = mConstitutiveLawOutput->Evaluate2D(rElement,
                                                   rIp,
                                                   copiedInputMap,
                                                   rConstitutiveOutput);
        break;

    case 3:
        error = mConstitutiveLawOutput->Evaluate3D(rElement,
                                                   rIp,
                                                   copiedInputMap,
                                                   rConstitutiveOutput);
        break;



    default:
        throw MechanicsException(__PRETTY_FUNCTION__,"invalid dimension");
    }

    for(auto itOutput : rConstitutiveOutput)
    {
        switch(itOutput.first)
        {
        case Constitutive::Output::D_ENGINEERING_STRESS_D_RELATIVE_HUMIDITY:
        {
            assert(rConstitutiveOutput.find(Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN)!=rConstitutiveOutput.end());
            ConstitutiveMatrix<VoigtDim, VoigtDim>& tangentStressStrain = *(static_cast<ConstitutiveMatrix<VoigtDim, VoigtDim>*>(rConstitutiveOutput.find(Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN)->second));
            if(d_EngineeringStrain_D_RH.GetIsCalculated() == false)
                throw MechanicsException(__PRETTY_FUNCTION__,std::string("Necessary value to determine ")+Constitutive::OutputToString(itOutput.first)+" was not calculated!");
            (static_cast<EngineeringStress<TDim>*>(itOutput.second))->AsVector() = tangentStressStrain.AsMatrix() * d_EngineeringStrain_D_RH.AsVector();
            break;
        }
        case Constitutive::Output::D_ENGINEERING_STRESS_D_WATER_VOLUME_FRACTION:
        {
            assert(rConstitutiveOutput.find(Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN)!=rConstitutiveOutput.end());
            ConstitutiveMatrix<VoigtDim, VoigtDim>& tangentStressStrain = *(static_cast<ConstitutiveMatrix<VoigtDim, VoigtDim>*>(rConstitutiveOutput.find(Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN)->second));
            if(d_EngineeringStrain_D_WV.GetIsCalculated() == false)
                throw MechanicsException(__PRETTY_FUNCTION__,std::string("Necessary value to determine ")+Constitutive::OutputToString(itOutput.first)+" was not calculated!");
            (static_cast<EngineeringStress<TDim>*>(itOutput.second))->AsVector() = tangentStressStrain.AsMatrix() * d_EngineeringStrain_D_WV.AsVector();
            break;
        }
        default:
            continue;
        }
        itOutput.second->SetIsCalculated(true);
    }
    return error;
}




NuTo::ConstitutiveInputMap NuTo::AdditiveInputExplicit::GetConstitutiveInputs(const NuTo::ConstitutiveOutputMap &rConstitutiveOutput,
                                                                                     const NuTo::InterpolationType &rInterpolationType) const
{
    ConstitutiveInputMap constitutiveInputMap;

    for(unsigned int i=0; i<mConstitutiveLawsModInput.size(); ++i)
    {

        ConstitutiveInputMap singleLawInputMap = mConstitutiveLawsModInput[i].first->GetConstitutiveInputs(rConstitutiveOutput,
                                                                                                           rInterpolationType);

        constitutiveInputMap.insert(singleLawInputMap.begin(),singleLawInputMap.end());
    }

    ConstitutiveInputMap singleLawInputMap = mConstitutiveLawOutput->GetConstitutiveInputs(rConstitutiveOutput,
                                                                                           rInterpolationType);

    constitutiveInputMap.insert(singleLawInputMap.begin(),singleLawInputMap.end());

    return constitutiveInputMap;
}



NuTo::Constitutive::Output::eOutput NuTo::AdditiveInputExplicit::GetOutputEnumFromInputEnum(NuTo::Constitutive::Input::eInput rInputEnum)
{
    switch(rInputEnum)
    {
    case Constitutive::Input::ENGINEERING_STRAIN:
        return Constitutive::Output::ENGINEERING_STRAIN;
    default:
        throw MechanicsException(__PRETTY_FUNCTION__,std::string("There is no output enum specified which is related to the input enum ")+Constitutive::InputToString(rInputEnum));
    }
}

void NuTo::AdditiveInputExplicit::AddCalculableDofCombinations(NuTo::ConstitutiveBase *rConstitutiveLaw)
{
    std::set<Node::eDof> allDofs = Node::GetDofSet();
    for (unsigned int i=0; i<mComputableDofCombinations.size(); ++i)
    for (auto itRow : allDofs)
        for (auto itCol : allDofs)
        {
            if (rConstitutiveLaw->CheckDofCombinationComputable(itRow,itCol,i))
                    mComputableDofCombinations[i].emplace(itRow,itCol);
        }
}
