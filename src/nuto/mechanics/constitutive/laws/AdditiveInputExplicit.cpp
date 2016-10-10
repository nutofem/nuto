#include "AdditiveInputExplicit.h"
#include "nuto/base/ErrorEnum.h"
#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "nuto/mechanics/constitutive/inputoutput/EngineeringStrain.h"
#include "nuto/mechanics/constitutive/inputoutput/EngineeringStress.h"
#include "nuto/mechanics/constitutive/staticData/ConstitutiveStaticDataMultipleConstitutiveLaws.h"
#include "nuto/mechanics/nodes/NodeEnum.h"






void NuTo::AdditiveInputExplicit::AddConstitutiveLaw(NuTo::ConstitutiveBase *rConstitutiveLaw, Constitutive::eInput rModiesInput)
{
    if(rConstitutiveLaw->HaveTmpStaticData())
        throw MechanicsException(__PRETTY_FUNCTION__,
                                 std::string("Constitutive law has tmp static data! The HaveTmpStaticData is only called on construction of the AdditiveInputExplicit law, but at this time, no constitutive law is attached")
                                 + "Therefore it does not know if it will have tmpstatic data or not and returns false by deafult. Find a way to update this information at the necessary code sections if a law with tmpstatic data is attached.");

    if(mStaticDataAllocated)
        throw MechanicsException(__PRETTY_FUNCTION__,
                                 "All constitutive laws have to be attached before static data is allocated!");

    if(rModiesInput == Constitutive::eInput::NONE)
    {
        if (mMainLaw != nullptr)
            throw MechanicsException(__PRETTY_FUNCTION__,
                                     std::string("There can be only one!!! --- This additive input law only accepts one law which calculates the output. All other laws ")+
                                     " are only allowed to modify the input to this law. Specify the modifying laws by providing the enum of the modified input as second function parameter.");
        mMainLaw = rConstitutiveLaw;
    }
    else
    {
        mSublaws.push_back({rConstitutiveLaw,rModiesInput});
    }
    AddCalculableDofCombinations(rConstitutiveLaw);
}

void NuTo::AdditiveInputExplicit::AddConstitutiveLaw(NuTo::ConstitutiveBase *rConstitutiveLaw)
{
    AddConstitutiveLaw(rConstitutiveLaw,Constitutive::eInput::NONE);
}








bool NuTo::AdditiveInputExplicit::CheckDofCombinationComputable(NuTo::Node::eDof rDofRow, NuTo::Node::eDof rDofCol, int rTimeDerivative) const
{
    if(mComputableDofCombinations[rTimeDerivative].find(std::pair<Node::eDof,Node::eDof>(rDofRow,rDofCol)) != mComputableDofCombinations[rTimeDerivative].end())
        return true;
    return false;
}

NuTo::ConstitutiveInputMap NuTo::AdditiveInputExplicit::GetConstitutiveInputs(const NuTo::ConstitutiveOutputMap &rConstitutiveOutput,
                                                                              const NuTo::InterpolationType &rInterpolationType) const
{
    // ------------------------------------------------
    // Get Inputs for output returning constitutive law
    // ------------------------------------------------

    ConstitutiveInputMap mainLawConstitutiveInputMap(mMainLaw->GetConstitutiveInputs(rConstitutiveOutput,
                                                                                     rInterpolationType));



    // -------------------------------------------------
    // Get Inputs for input modifying constitutive laws
    // -------------------------------------------------

    ConstitutiveInputMap sublawsConstitutiveInputMap;
    for (unsigned int i = 0; i < mSublaws.size(); ++i)
    {

        // Get the necessary inputs for the sublaws
        // ---> the modifications to the main laws inputs are returned as output from the sublaws.
        //      Therefore one first needs the sublaws output list, depending on the main laws inputs and the global outputs
        // INFO regarding the template function call:
        // The dimension can be chosen freely, because the created object is never used. Should be 1 because of the lowest construction costs.
        // Alternative would be to copy the called function and replace all objects with nullptr. But then you have to maintain 2 nearly identical functions.
        ConstitutiveOutputMap sublawOutputMap = GetSublawOutputMap<1>(mainLawConstitutiveInputMap,
                                                                      rConstitutiveOutput,
                                                                      i);

        // Don't merge the sublaw inputs directly into the main laws input map!
        // ---> When more than one sublaw is attached the inputs of the first might effect the following laws outputs - have a look at the line above!
        sublawsConstitutiveInputMap.Merge(mSublaws[i].first->GetConstitutiveInputs(sublawOutputMap,
                                                                                   rInterpolationType));
    }

    return mainLawConstitutiveInputMap.Merge(sublawsConstitutiveInputMap);
}

NuTo::Constitutive::eConstitutiveType NuTo::AdditiveInputExplicit::GetType() const
{
    return NuTo::Constitutive::eConstitutiveType::ADDITIVE_INPUT_EXPLICIT;
}








void NuTo::AdditiveInputExplicit::AddCalculableDofCombinations(NuTo::ConstitutiveBase *rConstitutiveLaw)
{
    // Not very smart, but works, feel free to do it better, but this function is only called during setup
    std::set<Node::eDof> allDofs = Node::GetDofSet();
    for (unsigned int i=0; i<mComputableDofCombinations.size(); ++i)    // vector -> time derivatives
        for (auto itRow : allDofs)
            for (auto itCol : allDofs)
                if (rConstitutiveLaw->CheckDofCombinationComputable(itRow,itCol,i))
                    mComputableDofCombinations[i].emplace(itRow,itCol);
}








template <int TDim>
NuTo::ConstitutiveStaticDataBase *NuTo::AdditiveInputExplicit::AllocateStaticDataAdditiveInputExplicit(const NuTo::ElementBase *rElement) const
{
    mStaticDataAllocated = true;    // <--- muteable member, so don't care about constness of this function

    std::vector<NuTo::ConstitutiveBase*> tempVec;
    for(unsigned int i=0; i<mSublaws.size(); ++i)
    {
        tempVec.push_back(mSublaws[i].first);
    }
    tempVec.push_back(mMainLaw);
    return new ConstitutiveStaticDataMultipleConstitutiveLaws(tempVec,rElement,TDim);
}







template <int TDim>
void NuTo::AdditiveInputExplicit::ApplySublawOutputs(const ConstitutiveInputMap& rMainLawInput,
                                                     const ConstitutiveOutputMap& rConstitutiveOutput,
                                                     const ConstitutiveOutputMap& rSublawOutput)
{
    for (const auto& it : rSublawOutput)
    {
        if (it.first == Constitutive::eOutput::ENGINEERING_STRAIN)
        {
            // Modify input strain for main constitutive law
            assert(rSublawOutput.at(Constitutive::eOutput::ENGINEERING_STRAIN)!=nullptr);
            assert(rMainLawInput.at(Constitutive::eInput::ENGINEERING_STRAIN)!=nullptr);
            *static_cast<EngineeringStrain<TDim>*>(rMainLawInput.at(Constitutive::eInput::ENGINEERING_STRAIN).get()) -= *static_cast<EngineeringStrain<TDim>*>(rSublawOutput.at(Constitutive::eOutput::ENGINEERING_STRAIN).get());
        }
        else if (rConstitutiveOutput.count(it.first) and rConstitutiveOutput.at(it.first) != nullptr)
        {
            //Copy sublaw output data direct into global output
            *(rConstitutiveOutput.at(it.first)) = *(it.second);
        }
    }
}






template <int TDim>
void NuTo::AdditiveInputExplicit::CalculateDerivatives(const ConstitutiveOutputMap& rConstitutiveOutput,
                                                       std::vector<ConstitutiveOutputMap>& rSublawOutputVec)
{
    constexpr const int VoigtDim = ConstitutiveIOBase::GetVoigtDim(TDim);
    for (auto& itOutput : rConstitutiveOutput)
    {
        switch(itOutput.first)
        {
        case Constitutive::eOutput::D_ENGINEERING_STRESS_D_RELATIVE_HUMIDITY:
        case Constitutive::eOutput::D_ENGINEERING_STRESS_D_WATER_VOLUME_FRACTION:
        case Constitutive::eOutput::D_ENGINEERING_STRESS_D_TEMPERATURE:
        {
            for(unsigned int i=0; i<mSublaws.size(); ++i)
            {
                if(mSublaws[i].second == Constitutive::eInput::ENGINEERING_STRAIN)
                {
                    // Get the corresponding sublaw output
                    Constitutive::eOutput derivativeOutputEnum = GetDerivativeEnumSublaw(Constitutive::eOutput::ENGINEERING_STRAIN,
                                                                                                 itOutput.first);

                    const ConstitutiveOutputMap::iterator& sublawOutput = rSublawOutputVec[i].find(derivativeOutputEnum);

                    if(sublawOutput == rSublawOutputVec[i].end())
                        // if current sublaw, does not provide the needed output, continue with next law
                        continue;
                    else
                    {
                        assert(rConstitutiveOutput.count(Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN));
                        assert(itOutput.second->GetIsCalculated() == false && "Currently it is not supported that multiple sublaws write to the same derivative");
                        if(sublawOutput->second->GetIsCalculated() == false)
                            throw MechanicsException(__PRETTY_FUNCTION__,
                                                     std::string("The value ") + Constitutive::OutputToString(sublawOutput->first) +", which is necessary to determine " +
                                                     Constitutive::OutputToString(itOutput.first) + " was requested  from a sublaw but has not been calculated!" );
                        const auto& tangentStressStrain = *static_cast<ConstitutiveMatrix<VoigtDim, VoigtDim>*>(rConstitutiveOutput.at(Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN).get());
                        const auto& sublawDerivative = *static_cast<ConstitutiveVector<VoigtDim>*>(sublawOutput->second.get());


                        (static_cast<EngineeringStress<TDim>*>(itOutput.second.get()))->AsVector() = tangentStressStrain * sublawDerivative;
                    }

                }
                else
                    // if law does not modify strains, continue with next sublaw
                    continue;

            }
            break;
        }

        default:
            continue;
        }
        itOutput.second->SetIsCalculated(true);
    }
}







template <int TDim>
NuTo::eError NuTo::AdditiveInputExplicit::EvaluateAdditiveInputExplicit(NuTo::ElementBase *rElement, int rIp,
                                                                               const NuTo::ConstitutiveInputMap &rConstitutiveInput,
                                                                               const NuTo::ConstitutiveOutputMap &rConstitutiveOutput)
{
    eError error = eError::SUCCESSFUL;

    // Copy inputs for main law, because they might be modified by the sublaws and these modifications will be passed above the borders of this law.
    NuTo::ConstitutiveInputMap mainLawInputMap = rConstitutiveInput;



    // ----------------
    // evaluate sublaws
    // ----------------

    std::vector<NuTo::ConstitutiveOutputMap> sublawOutputMapVec;

    for (unsigned int i = 0; i < mSublaws.size(); ++i)
    {
        // Get the sublaw specific output map depending on the main laws inputs and the global outputs
        sublawOutputMapVec.emplace_back(GetSublawOutputMap<TDim>(rConstitutiveInput,
                                                                 rConstitutiveOutput,
                                                                 i));

        // evaluate sublaw
        mSublaws[i].first->Evaluate<TDim>(rElement,
                                          rIp,
                                          rConstitutiveInput,
                                          sublawOutputMapVec[i]);
        if(error != eError::SUCCESSFUL)
            throw MechanicsException(__PRETTY_FUNCTION__,"Attached constitutive law returned an error code. Can't handle this");



        // Apply outputs to the main laws input and the global outputs
        ApplySublawOutputs<TDim>(mainLawInputMap,
                                 rConstitutiveOutput,
                                 sublawOutputMapVec[i]);

    }


    // -----------------
    // evaluate main law
    // -----------------

    error = mMainLaw->Evaluate<TDim>(rElement,
                                     rIp,
                                     mainLawInputMap,
                                     rConstitutiveOutput);



    // calculate derivatives that depend on outputs from the main law and the sublaws
    CalculateDerivatives<TDim>(rConstitutiveOutput,
                               sublawOutputMapVec);


    return error;
}







NuTo::Constitutive::eOutput NuTo::AdditiveInputExplicit::GetDerivativeEnumSublaw(NuTo::Constitutive::eOutput rParameter,
                                                                                         NuTo::Constitutive::eOutput rMainDerivative) const
{
    switch (rParameter)
    {
    case Constitutive::eOutput::ENGINEERING_STRAIN:
        switch(rMainDerivative)
        {
        case Constitutive::eOutput::D_ENGINEERING_STRESS_D_RELATIVE_HUMIDITY:
            return Constitutive::eOutput::D_ENGINEERING_STRAIN_D_RELATIVE_HUMIDITY;

        case Constitutive::eOutput::D_ENGINEERING_STRESS_D_WATER_VOLUME_FRACTION:
            return Constitutive::eOutput::D_ENGINEERING_STRAIN_D_WATER_VOLUME_FRACTION;

        case Constitutive::eOutput::D_ENGINEERING_STRESS_D_TEMPERATURE:
            return Constitutive::eOutput::D_STRAIN_D_TEMPERATURE;

        default:
            throw Exception(__PRETTY_FUNCTION__,
                            std::string("No partial derivative defined for parameter " + Constitutive::OutputToString(rParameter))
                            + " and global derivative " + Constitutive::OutputToString(rMainDerivative));
        }
        break;
    default:
        throw Exception(__PRETTY_FUNCTION__,
                        std::string("No partial derivatives defined for parameter " + Constitutive::OutputToString(rParameter)));
    }
}









template <int TDim>
NuTo::ConstitutiveOutputMap NuTo::AdditiveInputExplicit::GetSublawOutputMap(const NuTo::ConstitutiveInputMap& rMainLawInputMap,
                                                                                        const NuTo::ConstitutiveOutputMap& rMainLawOutputMap,
                                                                                        unsigned int rSublawIndex) const
{
    ConstitutiveOutputMap modifiedOutputMap;

    // ------------------------------------------------
    // Add sublaw outputs depending on main law outputs
    // ------------------------------------------------
    for(const auto& itMainLawOutput : rMainLawOutputMap)
    {
        switch(itMainLawOutput.first)
        {
        case Constitutive::eOutput::SHRINKAGE_STRAIN_VISUALIZE:
        case Constitutive::eOutput::THERMAL_STRAIN:
            modifiedOutputMap.emplace(itMainLawOutput.first,itMainLawOutput.second->clone());
            break;

        case Constitutive::eOutput::D_ENGINEERING_STRESS_D_RELATIVE_HUMIDITY:
        case Constitutive::eOutput::D_ENGINEERING_STRESS_D_WATER_VOLUME_FRACTION:
        case Constitutive::eOutput::D_ENGINEERING_STRESS_D_TEMPERATURE:
            if(mSublaws[rSublawIndex].second == Constitutive::eInput::ENGINEERING_STRAIN)
            {
                Constitutive::eOutput derivativeOutputEnum = GetDerivativeEnumSublaw(Constitutive::eOutput::ENGINEERING_STRAIN,
                                                                                             itMainLawOutput.first);
                 modifiedOutputMap.emplace(derivativeOutputEnum,
                                           ConstitutiveIOBase::makeConstitutiveIO<TDim>(derivativeOutputEnum));
            }
            break;

        default:
            break;
        }
    }


    // -----------------------------------------------
    // Add sublaw outputs depending on main law inputs
    // -----------------------------------------------
    switch(mSublaws[rSublawIndex].second)
    {
    case Constitutive::eInput::ENGINEERING_STRAIN:
        for(const auto& itMainLawInput : rMainLawInputMap)
        {
            switch(itMainLawInput.first)
            {
            case Constitutive::eInput::ENGINEERING_STRAIN:
                    modifiedOutputMap.emplace(Constitutive::eOutput::ENGINEERING_STRAIN,
                                              ConstitutiveIOBase::makeConstitutiveIO<TDim>(Constitutive::eOutput::ENGINEERING_STRAIN));
                break;

            default:
                break;
            }
        }

        break;

    default:
        break;
    }
    return modifiedOutputMap;
}







