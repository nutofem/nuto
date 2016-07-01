#include "ShrinkageCapillaryStrainBased.h"



#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveVector.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveScalar.h"
#include "nuto/physics/PhysicalConstantsSI.h"
#include "nuto/physics/PhysicalEquationsSI.h"





bool NuTo::ShrinkageCapillaryStrainBased::CheckDofCombinationComputable(NuTo::Node::eDof rDofRow, NuTo::Node::eDof rDofCol, int rTimeDerivative) const
{
    if(rTimeDerivative == 0 &&
       rDofRow == Node::DISPLACEMENTS &&
       (rDofCol == Node::RELATIVEHUMIDITY || rDofCol == Node::WATERVOLUMEFRACTION))
    {
        return true;
    }
    return false;
}

void NuTo::ShrinkageCapillaryStrainBased::CheckParameters() const
{

}





NuTo::ConstitutiveInputMap NuTo::ShrinkageCapillaryStrainBased::GetConstitutiveInputs(const NuTo::ConstitutiveOutputMap &rConstitutiveOutput,
                                                                                     const NuTo::InterpolationType &rInterpolationType) const
{
    ConstitutiveInputMap constitutiveInputMap;
    for (auto itOutput : rConstitutiveOutput)
    {
        switch (itOutput.first)
        {
        case Constitutive::Output::ENGINEERING_STRESS:      //VHIRTHAMTODO Temporary, because additiveInputExplicit isn't correct so far.
        case Constitutive::Output::ENGINEERING_STRAIN:
        case Constitutive::Output::ENGINEERING_STRAIN_VISUALIZE:
        case Constitutive::Output::ENGINEERING_STRESS_VISUALIZE:
        case Constitutive::Output::D_ENGINEERING_STRESS_D_RELATIVE_HUMIDITY:        //VHIRTHAMTODO must be the derivative of the strain instead of stress
        case Constitutive::Output::D_ENGINEERING_STRESS_D_WATER_VOLUME_FRACTION:    //VHIRTHAMTODO must be the derivative of the strain instead of stress
            constitutiveInputMap[Constitutive::Input::RELATIVE_HUMIDITY];
            constitutiveInputMap[Constitutive::Input::WATER_VOLUME_FRACTION];
            return constitutiveInputMap;

        default:
            continue;
        }
    }
    return constitutiveInputMap;
}

double NuTo::ShrinkageCapillaryStrainBased::GetParameterDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier) const
{
    switch(rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::MACROSCOPIC_BULK_MODULUS:
        return mMacroscopicBulkModulus;

    case Constitutive::eConstitutiveParameter::SOLID_PHASE_BULK_MODULUS:
        return mSolidPhaseModulus;

    case Constitutive::eConstitutiveParameter::TEMPERATURE:
        return mTemperature;

    default:
        throw MechanicsException(__PRETTY_FUNCTION__,std::string("Constitutive law does not have the parameter ")+Constitutive::ConstitutiveParameterToString(rIdentifier));
    }
}

void NuTo::ShrinkageCapillaryStrainBased::SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier, double rValue)
{
    switch(rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::MACROSCOPIC_BULK_MODULUS:
        mMacroscopicBulkModulus = rValue;
        return;

    case Constitutive::eConstitutiveParameter::SOLID_PHASE_BULK_MODULUS:
        mSolidPhaseModulus = rValue;
        return;

    case Constitutive::eConstitutiveParameter::TEMPERATURE:
        mTemperature = rValue;
        return;

    default:
        throw MechanicsException(__PRETTY_FUNCTION__,std::string("Constitutive law does not have the parameter ")+Constitutive::ConstitutiveParameterToString(rIdentifier));
    }
}

template <int TDim>
NuTo::Error::eError NuTo::ShrinkageCapillaryStrainBased::EvaluateShrinkageCapillary(NuTo::ElementBase *rElement,
                                                                                    int rIp,
                                                                                    const NuTo::ConstitutiveInputMap &rConstitutiveInput,
                                                                                    const NuTo::ConstitutiveOutputMap &rConstitutiveOutput)
{
    double relativeHumidity     = std::numeric_limits<double>::min();
    double waterVolumeFraction  = std::numeric_limits<double>::min();
    for (auto& itInput : rConstitutiveInput)
    {
        switch(itInput.first)
        {
        case Constitutive::Input::RELATIVE_HUMIDITY:
            relativeHumidity     = (*itInput.second)[0];
            break;

        case Constitutive::Input::WATER_VOLUME_FRACTION:
            waterVolumeFraction  = (*itInput.second)[0];
            break;

        default:
            continue;
        }
    }


    double bulkFactor = 1.0 / 3.0 * (1.0/mMacroscopicBulkModulus - 1.0/mSolidPhaseModulus);

    for (auto itOutput : rConstitutiveOutput)
    {

        switch(itOutput.first)
        {
        case NuTo::Constitutive::Output::ENGINEERING_STRAIN:
        {
            //Asserts
            itOutput.second->AssertIsVector<ConstitutiveIOBase::GetVoigtDim(TDim)>(itOutput.first,__PRETTY_FUNCTION__);
            assert(relativeHumidity    > std::numeric_limits<double>::min());
            assert(waterVolumeFraction > std::numeric_limits<double>::min());

            Eigen::Matrix<double, ConstitutiveIOBase::GetVoigtDim(TDim), 1>&  engineeringStrain = (*static_cast<ConstitutiveVector<ConstitutiveIOBase::GetVoigtDim(TDim)>*>(itOutput.second)).AsVector();
            //VHIRTHAMTODO --- how to handle atmospheric pressure?
            double capillaryStrain    = (//mAtmosphericPressure
                                         - waterVolumeFraction
                                         * NuTo::SI::DensityLiquidWater(mTemperature) * NuTo::SI::IdealGasConstant * mTemperature / NuTo::SI::MolarMassWater
                                         * log(relativeHumidity)) * bulkFactor;

            // The following loop does the same as multiplying the capillary stress with the kronecker delta in tensor form and adding the result to the engeneering stress
            for(unsigned int i=0; i<TDim; ++i)
            {
                engineeringStrain[i]+= capillaryStrain;
            }
            break;
        }


        case NuTo::Constitutive::Output::SHRINKAGE_STRAIN_VISUALIZE:
        {
            //Asserts
            itOutput.second->AssertIsVector<ConstitutiveIOBase::GetVoigtDim(3)>(itOutput.first,__PRETTY_FUNCTION__);
            assert(relativeHumidity    > std::numeric_limits<double>::min());
            assert(waterVolumeFraction > std::numeric_limits<double>::min());

            Eigen::Matrix<double, ConstitutiveIOBase::GetVoigtDim(3), 1>&  engineeringStrain = (*static_cast<ConstitutiveVector<ConstitutiveIOBase::GetVoigtDim(3)>*>(itOutput.second)).AsVector();
            //VHIRTHAMTODO --- how to handle atmospheric pressure?
            engineeringStrain.setZero();    //VHIRTHAMTODO Must be done on element level!
            double capillaryStrain    = (//mAtmosphericPressure
                                         - waterVolumeFraction
                                         * NuTo::SI::DensityLiquidWater(mTemperature) * NuTo::SI::IdealGasConstant * mTemperature / NuTo::SI::MolarMassWater
                                         * log(relativeHumidity)) * bulkFactor;

            // The following loop does the same as multiplying the capillary stress with the kronecker delta in tensor form and adding the result to the engeneering stress
            for(unsigned int i=0; i<TDim; ++i) //VHIRTHAMTODO Check if not 3
            {
                engineeringStrain[i]= capillaryStrain;
            }
            break;
        }

        case NuTo::Constitutive::Output::D_ENGINEERING_STRAIN_D_RELATIVE_HUMIDITY:
        {
            assert(relativeHumidity    > std::numeric_limits<double>::min());
            assert(waterVolumeFraction > std::numeric_limits<double>::min());

            Eigen::Matrix<double, ConstitutiveIOBase::GetVoigtDim(TDim), 1>&  engineeringStrain_dRH = (*static_cast<ConstitutiveVector<ConstitutiveIOBase::GetVoigtDim(TDim)>*>(itOutput.second)).AsVector();
            double capillaryStrain_dRH = - waterVolumeFraction
                                         * NuTo::SI::DensityLiquidWater(mTemperature) * NuTo::SI::IdealGasConstant * mTemperature
                                         / (NuTo::SI::MolarMassWater * relativeHumidity) * bulkFactor;
            // The following loop does the same as multiplying the capillary stress with the kronecker delta in tensor form and adding the result to the engeneering stress
            for(unsigned int i=0; i<TDim; ++i)
            {
                engineeringStrain_dRH[i]+= capillaryStrain_dRH;
            }
            break;
        }

        case NuTo::Constitutive::Output::D_ENGINEERING_STRAIN_D_WATER_VOLUME_FRACTION:
        {
            assert(relativeHumidity    > std::numeric_limits<double>::min());

            Eigen::Matrix<double, ConstitutiveIOBase::GetVoigtDim(TDim), 1>&  engineeringStrain_dWV = (*static_cast<ConstitutiveVector<ConstitutiveIOBase::GetVoigtDim(TDim)>*>(itOutput.second)).AsVector();
            double capillaryStrain_dWV = - NuTo::SI::DensityLiquidWater(mTemperature) * NuTo::SI::IdealGasConstant * mTemperature
                                         / (NuTo::SI::MolarMassWater)
                                         * log(relativeHumidity) * bulkFactor;
            // The following loop does the same as multiplying the capillary stress with the kronecker delta in tensor form and adding the result to the engeneering stress
            for(unsigned int i=0; i<TDim; ++i)
            {
                engineeringStrain_dWV[i]+= capillaryStrain_dWV;
            }
            break;
        }
        default:
            continue;
        }
        itOutput.second->SetIsCalculated(true);
    }
    return Error::SUCCESSFUL;
}
