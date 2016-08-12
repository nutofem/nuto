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
NuTo::Error::eError NuTo::AdditiveInputExplicit::EvaluateAdditiveInputExplicit(
        NuTo::ElementBase *rElement, int rIp,
        const NuTo::ConstitutiveInputMap &rConstitutiveInput,
        const NuTo::ConstitutiveOutputMap &rConstitutiveOutput)
{
    using namespace Constitutive;
    constexpr const int VoigtDim = ConstitutiveIOBase::GetVoigtDim(TDim);
    Error::eError error = Error::SUCCESSFUL;

    NuTo::ConstitutiveInputMap copiedInputMap = rConstitutiveInput;

    // TODO: Should be a vector of Outputs --> If more than one is attached, there are several possibilities at the moment,
    //       where values could be added multiple times, maybe overwritten but needed later or are just undefined!
    NuTo::ConstitutiveOutputMap modifiedOutputMap = rConstitutiveOutput;


    // TODO: Generalize and optimize!!!
    if(rConstitutiveInput.find(Input::ENGINEERING_STRAIN) != rConstitutiveInput.end())
    {
        modifiedOutputMap[Output::ENGINEERING_STRAIN] = ConstitutiveIOBase::makeConstitutiveIO<TDim>(Output::ENGINEERING_STRAIN);
    }


    // evaluate sublaws
    for (unsigned int i = 0; i < mConstitutiveLawsModInput.size(); ++i)
    {
        if(modifiedOutputMap.find(Output::D_ENGINEERING_STRESS_D_RELATIVE_HUMIDITY) != modifiedOutputMap.end() &&
                                  mConstitutiveLawsModInput[i].second == Input::ENGINEERING_STRAIN)
        {
            modifiedOutputMap[Output::D_ENGINEERING_STRAIN_D_RELATIVE_HUMIDITY] = ConstitutiveIOBase::makeConstitutiveIO<TDim>(Output::D_ENGINEERING_STRAIN_D_RELATIVE_HUMIDITY);
        }
        if(modifiedOutputMap.find(Output::D_ENGINEERING_STRESS_D_WATER_VOLUME_FRACTION) != modifiedOutputMap.end() &&
                                  mConstitutiveLawsModInput[i].second == Input::ENGINEERING_STRAIN)
        {
            modifiedOutputMap[Output::D_ENGINEERING_STRAIN_D_WATER_VOLUME_FRACTION] = ConstitutiveIOBase::makeConstitutiveIO<TDim>(Output::D_ENGINEERING_STRAIN_D_WATER_VOLUME_FRACTION);
        }
        if(modifiedOutputMap.find(Output::D_ENGINEERING_STRESS_D_TEMPERATURE) != modifiedOutputMap.end() &&
                                  mConstitutiveLawsModInput[i].second == Input::ENGINEERING_STRAIN)
        {
            modifiedOutputMap[Output::D_STRAIN_D_TEMPERATURE] = ConstitutiveIOBase::makeConstitutiveIO<TDim>(Output::D_STRAIN_D_TEMPERATURE);
        }



        // evaluate constitutive law
        mConstitutiveLawsModInput[i].first->Evaluate<TDim>(rElement, rIp, rConstitutiveInput, modifiedOutputMap);
        if(error != Error::SUCCESSFUL)
            throw MechanicsException(__PRETTY_FUNCTION__,"Attached constitutive law returned an error code. Can't handle this");



        for (const auto& it : modifiedOutputMap)
        {
            if (it.first == Output::ENGINEERING_STRAIN)
            {
                // Modify input strain for main constitutive law
                assert(modifiedOutputMap.at(Output::ENGINEERING_STRAIN)!=nullptr);
                assert(copiedInputMap.at(Input::ENGINEERING_STRAIN)!=nullptr);
                *static_cast<EngineeringStrain<TDim>*>(copiedInputMap.at(Input::ENGINEERING_STRAIN).get()) -= *static_cast<EngineeringStrain<TDim>*>(modifiedOutputMap.at(Output::ENGINEERING_STRAIN).get());
            }
            else if (rConstitutiveOutput.count(it.first) and rConstitutiveOutput.at(it.first) != nullptr)
            {
                *(rConstitutiveOutput.at(it.first)) = *(it.second);
            }
        }

    }

    // evaluate output law
    error = mConstitutiveLawOutput->Evaluate<TDim>(rElement, rIp, copiedInputMap, rConstitutiveOutput);


    // evaluate derivatives of outputs depending on sublaws
    for (auto& itOutput : rConstitutiveOutput)
    {
        switch(itOutput.first)
        {
        case Output::D_ENGINEERING_STRESS_D_RELATIVE_HUMIDITY:
        {
            assert(rConstitutiveOutput.count(Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN));
            const auto& tangentStressStrain = *static_cast<ConstitutiveMatrix<VoigtDim, VoigtDim>*>(rConstitutiveOutput.at(Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN).get());
            const auto& d_EngineeringStrain_D_RH = *static_cast<ConstitutiveVector<VoigtDim>*>(modifiedOutputMap.at(Output::D_ENGINEERING_STRAIN_D_RELATIVE_HUMIDITY).get());
            if(d_EngineeringStrain_D_RH.GetIsCalculated() == false)
                throw MechanicsException(__PRETTY_FUNCTION__,std::string("Necessary value to determine ")+OutputToString(itOutput.first)+" was not calculated!");
            (static_cast<EngineeringStress<TDim>*>(itOutput.second.get()))->AsVector() = tangentStressStrain * d_EngineeringStrain_D_RH;
            break;
        }
        case Output::D_ENGINEERING_STRESS_D_WATER_VOLUME_FRACTION:
        {
            assert(rConstitutiveOutput.count(Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN));
            const auto& tangentStressStrain = *static_cast<ConstitutiveMatrix<VoigtDim, VoigtDim>*>(rConstitutiveOutput.at(Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN).get());
            const auto& d_EngineeringStrain_D_WV = *static_cast<ConstitutiveVector<VoigtDim>*>(modifiedOutputMap.at(Output::D_ENGINEERING_STRAIN_D_WATER_VOLUME_FRACTION).get());
            if(d_EngineeringStrain_D_WV.GetIsCalculated() == false)
                throw MechanicsException(__PRETTY_FUNCTION__,std::string("Necessary value to determine ")+OutputToString(itOutput.first)+" was not calculated!");
            (static_cast<EngineeringStress<TDim>*>(itOutput.second.get()))->AsVector() = tangentStressStrain * d_EngineeringStrain_D_WV;
            break;
        }
        case Output::D_ENGINEERING_STRESS_D_TEMPERATURE:
        {
            assert(rConstitutiveOutput.count(Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN));
            const auto& tangentStressStrain = *static_cast<ConstitutiveMatrix<VoigtDim, VoigtDim>*>(rConstitutiveOutput.at(Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN).get());
            const auto& d_Strain_D_Temperature = *static_cast<ConstitutiveVector<VoigtDim>*>(modifiedOutputMap.at(Output::D_STRAIN_D_TEMPERATURE).get());
            if(d_Strain_D_Temperature.GetIsCalculated() == false)
                throw MechanicsException(__PRETTY_FUNCTION__, "Necessary value to determine " + OutputToString(itOutput.first) + " was not calculated!");
            (static_cast<EngineeringStress<TDim>*>(itOutput.second.get()))->AsVector() = tangentStressStrain * d_Strain_D_Temperature;
            break;
        }
        default:
            continue;
        }
        itOutput.second->SetIsCalculated(true);
    }
    return error;
}



NuTo::ConstitutiveInputMap NuTo::AdditiveInputExplicit::GetConstitutiveInputs(
        const NuTo::ConstitutiveOutputMap &rConstitutiveOutput,
        const NuTo::InterpolationType &rInterpolationType) const
{
    ConstitutiveInputMap constitutiveInputMap;



    // ------------------------------------------------
    // Get Inputs for output returning constitutive law
    // ------------------------------------------------

    ConstitutiveInputMap outputLawInputMap = mConstitutiveLawOutput->GetConstitutiveInputs(rConstitutiveOutput,
                                                                                           rInterpolationType);
    constitutiveInputMap.Merge(outputLawInputMap);



    // -------------------------------------------------
    // Get Inputs for input modifying constitutive laws
    // -------------------------------------------------

    // TODO: Generalize!!!

    // The input modifying laws return their modifications as output. Therefore, if the main constitutive law
    // needs a modifyable input, the other laws have to specify which inputs they need to calculate their modifcations.
    ConstitutiveOutputMap modifiedOutputMap = rConstitutiveOutput;
    if(constitutiveInputMap.find(Constitutive::Input::ENGINEERING_STRAIN) != constitutiveInputMap.end())
        modifiedOutputMap.emplace(Constitutive::Output::ENGINEERING_STRAIN,nullptr);

      for (unsigned int i = 0; i < mConstitutiveLawsModInput.size(); ++i)
    {

        ConstitutiveInputMap singleLawInputMap = mConstitutiveLawsModInput[i].first->GetConstitutiveInputs(modifiedOutputMap,
                                                                                                           rInterpolationType);
        constitutiveInputMap.Merge(singleLawInputMap);
    }

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
