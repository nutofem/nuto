#include "ShrinkageCapillaryStressBased.h"



#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveVector.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveScalar.h"






bool NuTo::ShrinkageCapillaryStressBased::CheckDofCombinationComputable(NuTo::Node::eDof rDofRow, NuTo::Node::eDof rDofCol, int rTimeDerivative) const
{
    if(rTimeDerivative == 0 &&
       rDofRow == Node::DISPLACEMENTS &&
       (rDofCol == Node::RELATIVEHUMIDITY || rDofCol == Node::WATERVOLUMEFRACTION))
    {
        return true;
    }
    return false;
}

void NuTo::ShrinkageCapillaryStressBased::CheckParameters() const
{

}





NuTo::ConstitutiveInputMap NuTo::ShrinkageCapillaryStressBased::GetConstitutiveInputs(const NuTo::ConstitutiveOutputMap &rConstitutiveOutput,
                                                                                     const NuTo::InterpolationType &rInterpolationType) const
{
    ConstitutiveInputMap constitutiveInputMap;
    for (auto itOutput : rConstitutiveOutput)
    {
        switch (itOutput.first)
        {
        case Constitutive::Output::ENGINEERING_STRESS:
        case Constitutive::Output::D_ENGINEERING_STRESS_D_RELATIVE_HUMIDITY:
        case Constitutive::Output::D_ENGINEERING_STRESS_D_WATER_VOLUME_FRACTION:
            constitutiveInputMap[Constitutive::Input::RELATIVE_HUMIDITY];
            constitutiveInputMap[Constitutive::Input::WATER_VOLUME_FRACTION];
            return constitutiveInputMap;

        default:
            continue;
        }
    }
    return constitutiveInputMap;
}

template <int TDim>
NuTo::Error::eError NuTo::ShrinkageCapillaryStressBased::EvaluateShrinkageCapillary(NuTo::ElementBase *rElement,
                                                                                    int rIp,
                                                                                    const NuTo::ConstitutiveInputMap &rConstitutiveInput,
                                                                                    const NuTo::ConstitutiveOutputMap &rConstitutiveOutput)
{
    double relativeHumidity     = std::numeric_limits<double>::min();
    double waterVolumeFraction  = std::numeric_limits<double>::min();
    for (auto itInput : rConstitutiveInput)
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


    for (auto itOutput : rConstitutiveOutput)
    {

        switch(itOutput.first)
        {
        case NuTo::Constitutive::Output::ENGINEERING_STRESS:
        {
            //Asserts
            itOutput.second->AssertIsVector<ConstitutiveIOBase::GetVoigtDim(TDim)>(itOutput.first,__PRETTY_FUNCTION__);
            assert(relativeHumidity    > std::numeric_limits<double>::min());
            assert(waterVolumeFraction > std::numeric_limits<double>::min());

            Eigen::Matrix<double, ConstitutiveIOBase::GetVoigtDim(TDim), 1>&  engineeringStress = (*static_cast<ConstitutiveVector<ConstitutiveIOBase::GetVoigtDim(TDim)>*>(itOutput.second)).AsVector();
            //VHIRTHAMTODO --- how to handle atmospheric pressure?
            double capillaryStress    = (//mAtmosphericPressure
                                         - waterVolumeFraction / mPorosity
                                         * mDensityWater * mIdealGasConstant * mTemperature / mMolarMassWater
                                         * log(relativeHumidity));

            // The following loop does the same as multiplying the capillary stress with the kronecker delta in tensor form and adding the result to the engeneering stress
            for(unsigned int i=0; i<TDim; ++i)
            {
                engineeringStress[i]+= capillaryStress;
            }
            break;
        }

        case NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_RELATIVE_HUMIDITY:
        {
            assert(relativeHumidity    > std::numeric_limits<double>::min());
            assert(waterVolumeFraction > std::numeric_limits<double>::min());

            Eigen::Matrix<double, ConstitutiveIOBase::GetVoigtDim(TDim), 1>&  engineeringStress_dRH = (*static_cast<ConstitutiveVector<ConstitutiveIOBase::GetVoigtDim(TDim)>*>(itOutput.second)).AsVector();
            double capillaryStress_dRH = - waterVolumeFraction / mPorosity
                                         * mDensityWater * mIdealGasConstant * mTemperature
                                         / (mMolarMassWater * relativeHumidity);
            // The following loop does the same as multiplying the capillary stress with the kronecker delta in tensor form and adding the result to the engeneering stress
            for(unsigned int i=0; i<TDim; ++i)
            {
                engineeringStress_dRH[i]+= capillaryStress_dRH;
            }
            break;
        }

        case NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_WATER_VOLUME_FRACTION:
        {
            assert(relativeHumidity    > std::numeric_limits<double>::min());

            Eigen::Matrix<double, ConstitutiveIOBase::GetVoigtDim(TDim), 1>&  engineeringStress_dWV = (*static_cast<ConstitutiveVector<ConstitutiveIOBase::GetVoigtDim(TDim)>*>(itOutput.second)).AsVector();
            double capillaryStress_dWV = - mDensityWater * mIdealGasConstant * mTemperature
                                         / (mMolarMassWater * mPorosity)
                                         * log(relativeHumidity);
            // The following loop does the same as multiplying the capillary stress with the kronecker delta in tensor form and adding the result to the engeneering stress
            for(unsigned int i=0; i<TDim; ++i)
            {
                engineeringStress_dWV[i]+= capillaryStress_dWV;
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
