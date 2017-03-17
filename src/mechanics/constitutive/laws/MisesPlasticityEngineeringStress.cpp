#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/utility.hpp>
#endif // ENABLE_SERIALIZATION

#include <iostream>
#include "base/Logger.h"

#include "mechanics/structures/StructureBase.h"
#include "mechanics/MechanicsException.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/constitutive/laws/MisesPlasticityEngineeringStress.h"
#include "mechanics/constitutive/staticData/DataMisesPlasticity.h"
#include "mechanics/constitutive/laws/EngineeringStressHelper.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "mechanics/constitutive/inputoutput/EngineeringStrain.h"
#include "mechanics/constitutive/inputoutput/EngineeringStress.h"
#include "mechanics/interpolationtypes/InterpolationType.h"
#include "mechanics/MechanicsException.h"
#include "mechanics/elements/ElementBase.h"
#include "mechanics/elements/ElementEnum.h"
#include "mechanics/nodes/NodeEnum.h"

using namespace NuTo::Constitutive;

NuTo::MisesPlasticityEngineeringStress::MisesPlasticityEngineeringStress() : ConstitutiveBase()
{
    mE = 0.;
    mNu = 0.;
    mRho = 0.;
    mSigma.resize(1);
    mH.resize(1);
    SetParametersValid();
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::MisesPlasticityEngineeringStress::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::MisesPlasticityEngineeringStress::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::MisesPlasticityEngineeringStress::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::MisesPlasticityEngineeringStress::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::MisesPlasticityEngineeringStress::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::MisesPlasticityEngineeringStress::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::MisesPlasticityEngineeringStress::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize MisesPlasticityEngineeringStress" << "\n";
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveBase)
       & BOOST_SERIALIZATION_NVP(mE)
       & BOOST_SERIALIZATION_NVP(mNu)
       & BOOST_SERIALIZATION_NVP(mSigma)
       & BOOST_SERIALIZATION_NVP(mH)
       & BOOST_SERIALIZATION_NVP(mRho)
       & BOOST_SERIALIZATION_NVP(mThermalExpansionCoefficient);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize MisesPlasticityEngineeringStress" << "\n";
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::MisesPlasticityEngineeringStress)
#endif // ENABLE_SERIALIZATION

NuTo::ConstitutiveInputMap NuTo::MisesPlasticityEngineeringStress::GetConstitutiveInputs(
        const ConstitutiveOutputMap& rConstitutiveOutput, const InterpolationType& rInterpolationType) const
{
    ConstitutiveInputMap constitutiveInputMap;

    constitutiveInputMap[Constitutive::eInput::ENGINEERING_STRAIN];
    if (rInterpolationType.IsConstitutiveInput(Node::eDof::TEMPERATURE))
        constitutiveInputMap[Constitutive::eInput::TEMPERATURE];


    for (auto& itOutput : rConstitutiveOutput)
    {
        switch (itOutput.first)
        {
        case NuTo::Constitutive::eOutput::ENGINEERING_STRESS:
        case NuTo::Constitutive::eOutput::ENGINEERING_STRESS_VISUALIZE:
        case NuTo::Constitutive::eOutput::ENGINEERING_STRAIN_VISUALIZE:
        case NuTo::Constitutive::eOutput::ENGINEERING_PLASTIC_STRAIN_VISUALIZE:
        case NuTo::Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN:
        case NuTo::Constitutive::eOutput::UPDATE_TMP_STATIC_DATA:
        case NuTo::Constitutive::eOutput::UPDATE_STATIC_DATA:
            break;
        default:
            throw MechanicsException(__PRETTY_FUNCTION__,
                    "Output object " + OutputToString(itOutput.first) + " cannot be calculated by "
                    "this constitutive law.");
        }
    }

    return constitutiveInputMap;
}

namespace NuTo
{

template<>
void NuTo::MisesPlasticityEngineeringStress::Evaluate<1>(
    const ConstitutiveInputMap& rConstitutiveInput,
    const ConstitutiveOutputMap& rConstitutiveOutput,
    Data& rStaticData)
{
    throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "not implemented for 1D.");
}

template<>
void NuTo::MisesPlasticityEngineeringStress::Evaluate<2>(
    const ConstitutiveInputMap& rConstitutiveInput,
    const ConstitutiveOutputMap& rConstitutiveOutput,
    Data& rStaticData)
{

    const auto& planeState =
        *dynamic_cast<ConstitutivePlaneState*>(rConstitutiveInput.at(Constitutive::eInput::PLANE_STATE).get());
    if (planeState.GetPlaneState() != ePlaneState::PLANE_STRAIN)
        throw MechanicsException(__PRETTY_FUNCTION__, "Only implemented for PLANE_STRAIN");

    const auto& engineeringStrain =
        rConstitutiveInput.at(Constitutive::eInput::ENGINEERING_STRAIN)->AsEngineeringStrain2D();

    StaticData::DataMisesPlasticity<3> newStaticData;
    StaticData::DataMisesPlasticity<3>* newStaticDataPtr = nullptr;

    EngineeringStress<2> engineeringStress;
    ConstitutiveIOBase* engineeringStressPtr = nullptr;
    ConstitutiveIOBase* tangent = nullptr;

    bool strainRequested = false;

    for (auto& itOutput : rConstitutiveOutput)
    {
        switch (itOutput.first)
        {
            case Constitutive::eOutput::ENGINEERING_STRESS:
            case Constitutive::eOutput::ENGINEERING_STRESS_VISUALIZE:
            {
                engineeringStressPtr = &engineeringStress;
                break;
            }

            case Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN:
            {
                tangent = itOutput.second.get();
                tangent->AssertIsMatrix<3, 3>(Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN,
                                              __FUNCTION__);
                break;
            }

            case Constitutive::eOutput::ENGINEERING_STRAIN_VISUALIZE:
                strainRequested = true;
                break;

            case Constitutive::eOutput::ENGINEERING_PLASTIC_STRAIN_VISUALIZE:
            case Constitutive::eOutput::UPDATE_STATIC_DATA:
                newStaticDataPtr = &newStaticData;
                continue;

            case Constitutive::eOutput::UPDATE_TMP_STATIC_DATA:

                continue;
            default:
                continue;
        }
        itOutput.second->SetIsCalculated(true);
    }


    if (rConstitutiveOutput.size() == 1 && strainRequested)
    {
        // return mapping can skipped, if ENGINEERING_STRAIN_VISUALIZE is the only requested output.
    }
    else
    {
        ReturnMapping2D(rStaticData.GetData(), engineeringStrain, engineeringStressPtr, tangent, newStaticDataPtr);
    }

    for (auto& itOutput : rConstitutiveOutput)
    {
        switch (itOutput.first)
        {
            case Constitutive::eOutput::ENGINEERING_STRESS:
            {
                ConstitutiveIOBase& engineeringStress2D = *itOutput.second;
                engineeringStress2D.AssertIsVector<3>(itOutput.first, __PRETTY_FUNCTION__);
                engineeringStress2D = engineeringStress;
                break;
            }
            case Constitutive::eOutput::ENGINEERING_STRESS_VISUALIZE:
            {
                ConstitutiveIOBase& engineeringStress3D = *itOutput.second;
                engineeringStress3D.AssertIsVector<6>(itOutput.first, __PRETTY_FUNCTION__);
                engineeringStress3D.SetZero();
                engineeringStress3D[0] = engineeringStress[0];
                engineeringStress3D[1] = engineeringStress[1];
                engineeringStress3D[5] = engineeringStress[2];
                break;
            }
            case Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN: // calculated via ptr in return mapping
                break;
            case Constitutive::eOutput::ENGINEERING_STRAIN_VISUALIZE:
            {
                ConstitutiveIOBase& engineeringStrain3D = *itOutput.second;
                engineeringStrain3D.AssertIsVector<6>(itOutput.first, __PRETTY_FUNCTION__);
                engineeringStrain3D.SetZero();
                engineeringStrain3D[0] = (engineeringStrain)[0];
                engineeringStrain3D[1] = (engineeringStrain)[1];
                engineeringStrain3D[5] = (engineeringStrain)[2];
                break;
            }
            case Constitutive::eOutput::ENGINEERING_PLASTIC_STRAIN_VISUALIZE:
            {
                auto& engineeringPlasticStrain = *dynamic_cast<EngineeringStrain<3>*>(itOutput.second.get());
                engineeringPlasticStrain = newStaticData.GetPlasticStrain();
                break;
            }
            case Constitutive::eOutput::UPDATE_TMP_STATIC_DATA:
                continue;
            case Constitutive::eOutput::UPDATE_STATIC_DATA:
            {
                rStaticData.SetData(newStaticData);
            }
                continue;
            default:
                continue;
        }
        itOutput.second->SetIsCalculated(true);
    }
}

template<>
void NuTo::MisesPlasticityEngineeringStress::Evaluate<3>(
    const ConstitutiveInputMap& rConstitutiveInput,
    const ConstitutiveOutputMap& rConstitutiveOutput,
    Data& rStaticData)
{

    const auto
        & engineeringStrain = rConstitutiveInput.at(Constitutive::eInput::ENGINEERING_STRAIN)->AsEngineeringStrain3D();

    ConstitutiveIOBase* engineeringStress = nullptr;
    ConstitutiveIOBase* tangent = nullptr;

    StaticData::DataMisesPlasticity<3> newStaticData;
    StaticData::DataMisesPlasticity<3>* newStaticDataPtr = nullptr;

    bool strainRequested = false;

    for (auto& itOutput : rConstitutiveOutput)
    {
        switch (itOutput.first)
        {
            case Constitutive::eOutput::ENGINEERING_STRESS:
            {
                engineeringStress = itOutput.second.get();
                engineeringStress->AssertIsVector<6>(Constitutive::eOutput::ENGINEERING_STRESS, __FUNCTION__);
                break;
            }

            case Constitutive::eOutput::ENGINEERING_STRESS_VISUALIZE:
            {
                // use the ENGINEERING_STRESS for visualization. if not, set the pointer
                if (rConstitutiveOutput.find(Constitutive::eOutput::ENGINEERING_STRESS) == rConstitutiveOutput.end())
                {
                    engineeringStress = itOutput.second.get();
                    engineeringStress
                        ->AssertIsVector<6>(Constitutive::eOutput::ENGINEERING_STRESS_VISUALIZE, __FUNCTION__);
                }
                break;
            }

            case Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN:
            {
                tangent = itOutput.second.get();
                tangent->AssertIsMatrix<6, 6>(Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN,
                                              __FUNCTION__);
                break;
            }
            case Constitutive::eOutput::ENGINEERING_STRAIN_VISUALIZE:
                strainRequested = true;
                break;
            case Constitutive::eOutput::ENGINEERING_PLASTIC_STRAIN_VISUALIZE:
            case Constitutive::eOutput::UPDATE_STATIC_DATA:
                newStaticDataPtr = &newStaticData;
                continue;
            case Constitutive::eOutput::UPDATE_TMP_STATIC_DATA:
                continue;
            default:
                continue;
        }
        itOutput.second->SetIsCalculated(true);
    }


    if (rConstitutiveOutput.size() == 1 && strainRequested)
    {
        // return mapping can skipped, if ENGINEERING_STRAIN_VISUALIZE is the only requested output.
    }
    else
    {
        ReturnMapping3D(rStaticData.GetData(), engineeringStrain, engineeringStress, tangent, newStaticDataPtr);
    }

    for (auto& itOutput : rConstitutiveOutput)
    {
        switch (itOutput.first)
        {
            case Constitutive::eOutput::ENGINEERING_STRESS:
            case Constitutive::eOutput::ENGINEERING_STRESS_VISUALIZE:
            {
                ConstitutiveIOBase& engineeringStress3D = *itOutput.second;
                engineeringStress3D = *engineeringStress;
                break;
            }
            case Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN:     // calculated via ptr in return mapping
                break;
            case Constitutive::eOutput::ENGINEERING_STRAIN_VISUALIZE:
            {
                ConstitutiveIOBase& engineeringStrain3D = *itOutput.second;
                engineeringStrain3D.AssertIsVector<6>(itOutput.first, __PRETTY_FUNCTION__);
                engineeringStrain3D = engineeringStrain;
                break;
            }
            case Constitutive::eOutput::ENGINEERING_PLASTIC_STRAIN_VISUALIZE:
            {
                auto& engineeringPlasticStrain = *dynamic_cast<EngineeringStrain<3>*>(itOutput.second.get());
                engineeringPlasticStrain = newStaticData.GetPlasticStrain();
                break;
            }
            case Constitutive::eOutput::UPDATE_TMP_STATIC_DATA:
                continue;
            case Constitutive::eOutput::UPDATE_STATIC_DATA:
            {
                rStaticData.SetData(newStaticData);
            }
                continue;
            default:
                continue;
        }
        itOutput.second->SetIsCalculated(true);
    }
}

} // namespace NuTo


bool NuTo::MisesPlasticityEngineeringStress::CheckDofCombinationComputable(Node::eDof rDofRow,
                                                                            Node::eDof rDofCol,
                                                                            int rTimeDerivative) const
{
    assert(rTimeDerivative>-1);
    if(rTimeDerivative<1 &&
            rDofRow == Node::eDof::DISPLACEMENTS &&
            rDofCol == Node::eDof::DISPLACEMENTS)
    {
        return true;
    }
    return false;
}


void NuTo::MisesPlasticityEngineeringStress::ReturnMapping2D(
        StaticData::DataMisesPlasticity<3>& oldStaticData,
        const EngineeringStrain<2>& rEngineeringStrain,
        ConstitutiveIOBase* rNewStress,
        ConstitutiveIOBase* rNewTangent,
        StaticData::DataMisesPlasticity<3>* rNewStaticData) const
{

    EngineeringStrain<3> engineeringStrain3D;
    engineeringStrain3D.SetZero();

    engineeringStrain3D[0] = rEngineeringStrain[0];
    engineeringStrain3D[1] = rEngineeringStrain[1];
    engineeringStrain3D[3] = rEngineeringStrain[2];




    EngineeringStress<3> newStress3D;
    ConstitutiveMatrix<6,6> newTangent3D;

    ReturnMapping3D(oldStaticData, engineeringStrain3D, &newStress3D, &newTangent3D, rNewStaticData);


    if (rNewStress)
    {
        (*rNewStress)[0] = newStress3D[0];
        (*rNewStress)[1] = newStress3D[1];
        (*rNewStress)[2] = newStress3D[3];
    }


    if (rNewTangent)
    {
        (*rNewTangent)(0,0) = newTangent3D(0,0);
        (*rNewTangent)(0,1) = newTangent3D(0,1);
        (*rNewTangent)(0,2) = newTangent3D(0,3);

        (*rNewTangent)(1,0) = newTangent3D(1,0);
        (*rNewTangent)(1,1) = newTangent3D(1,1);
        (*rNewTangent)(1,2) = newTangent3D(1,3);

        (*rNewTangent)(2,0) = newTangent3D(3,0);
        (*rNewTangent)(2,1) = newTangent3D(3,1);
        (*rNewTangent)(2,2) = newTangent3D(3,3);
    }
}

void NuTo::MisesPlasticityEngineeringStress::ReturnMapping3D(
        StaticData::DataMisesPlasticity<3>& oldStaticData,
        const EngineeringStrain<3>& rEngineeringStrain,
        ConstitutiveIOBase* rNewStress,
        ConstitutiveIOBase* rNewTangent,
        StaticData::DataMisesPlasticity<3>* rNewStaticData) const
{
    const double sqrt_2div3 = 0.81649658;
    const double tolerance = 1e-8;

    double sigma_trial[6],
    xi_trial[6],
    norm_dev,
    sigma_y,
    factor,
    factor2,
    yield_condition,
    epsilon_p_eq2,
    d_sigma,
    d_H,
    H,
    H2,
    mu,
    bulk_modulus,
    delta_gamma=0.,
    g,
    dg,
    df_dsigma[6],
    trace_epsilon,
    trace_epsilon_div_3;

    mu    = mE/(2.*(1.+mNu));
    bulk_modulus = mE/(3.-6.*mNu);

    // set strain data ptr
    const EngineeringStrain<3>& total_strain = rEngineeringStrain;

    const EngineeringStrain<3> plastic_strain = oldStaticData.GetPlasticStrain();
    const EngineeringStress<3> back_stress = oldStaticData.GetBackStress();

    trace_epsilon = total_strain[0] + total_strain[1] + total_strain[2];

    trace_epsilon_div_3 = trace_epsilon/3.;

    //trial stress
    sigma_trial[0] = (total_strain[0]-trace_epsilon_div_3-plastic_strain[0])*2.*mu;
    sigma_trial[1] = (total_strain[1]-trace_epsilon_div_3-plastic_strain[1])*2.*mu;
    sigma_trial[2] = (total_strain[2]-trace_epsilon_div_3-plastic_strain[2])*2.*mu;
    sigma_trial[3] = (total_strain[3]                    -plastic_strain[3])*mu; // in total strain, gamma is stored
    sigma_trial[4] = (total_strain[4]                    -plastic_strain[4])*mu; // in total strain, gamma is stored
    sigma_trial[5] = (total_strain[5]                    -plastic_strain[5])*mu; // in total strain, gamma is stored

    //subtract backstress
    xi_trial[0] = sigma_trial[0]-back_stress[0];
    xi_trial[1] = sigma_trial[1]-back_stress[1];
    xi_trial[2] = sigma_trial[2]-back_stress[2];
    xi_trial[3] = sigma_trial[3]-back_stress[3];
    xi_trial[4] = sigma_trial[4]-back_stress[4];
    xi_trial[5] = sigma_trial[5]-back_stress[5];

//    printf("total_strain %g %g %g %g %g %g\n",total_strain[0],total_strain[1],total_strain[2],total_strain[3],total_strain[4],total_strain[5]);
//    printf("Xi trial %g %g %g %g %g %g\n",xi_trial[0],xi_trial[1],xi_trial[2],xi_trial[3],xi_trial[4],xi_trial[5]);

    // norm of deviator
    norm_dev = std::sqrt(xi_trial[0]*xi_trial[0]+xi_trial[1]*xi_trial[1]+xi_trial[2]*xi_trial[2]+
                    2.*(xi_trial[3]*xi_trial[3]+xi_trial[4]*xi_trial[4]+xi_trial[5]*xi_trial[5]));

    //determine radius of yield function
    sigma_y = GetYieldStrength(oldStaticData.GetEquivalentPlasticStrain(), d_sigma);

    yield_condition = norm_dev - sqrt_2div3 * sigma_y;

    if (yield_condition<-tolerance*sigma_y)
    {
        // elastic regime
    	factor = bulk_modulus*trace_epsilon;
        if (rNewStress!=0)
        {
        	(*rNewStress)[0] = factor+sigma_trial[0];
        	(*rNewStress)[1] = factor+sigma_trial[1];
        	(*rNewStress)[2] = factor+sigma_trial[2];
        	(*rNewStress)[3] = 	   sigma_trial[3];
        	(*rNewStress)[4] = 	   sigma_trial[4];
        	(*rNewStress)[5] = 	   sigma_trial[5];
        }
        if (rNewTangent!=0)
        {
            factor = mE/(1.+mNu)/(1.-2.*mNu);

            (*rNewTangent)(0,0) = (1.-mNu)*factor;
            (*rNewTangent)(1,0) = mNu*factor;
            (*rNewTangent)(2,0) = (*rNewTangent)(1,0);
            (*rNewTangent)(3,0) = 0.;
            (*rNewTangent)(4,0) = 0.;
            (*rNewTangent)(5,0) = 0.;

            (*rNewTangent)(0,1) = (*rNewTangent)(1,0);
            (*rNewTangent)(1,1) = (*rNewTangent)(0,0);
            (*rNewTangent)(2,1) = (*rNewTangent)(1,0);
            (*rNewTangent)(3,1) = 0.;
            (*rNewTangent)(4,1) = 0.;
            (*rNewTangent)(5,1) = 0.;

            (*rNewTangent)(0,2) = (*rNewTangent)(1,0);
            (*rNewTangent)(1,2) = (*rNewTangent)(1,0);
            (*rNewTangent)(2,2) = (*rNewTangent)(0,0);
            (*rNewTangent)(3,2) = 0.;
            (*rNewTangent)(4,2) = 0.;
            (*rNewTangent)(5,2) = 0.;

            (*rNewTangent)(0,3) = 0.;
            (*rNewTangent)(1,3) = 0.;
            (*rNewTangent)(2,3) = 0.;
            (*rNewTangent)(3,3) = (0.5-mNu)*factor;
            (*rNewTangent)(4,3) = 0.;
            (*rNewTangent)(5,3) = 0.;

            (*rNewTangent)(0,4) = 0.;
            (*rNewTangent)(1,4) = 0.;
            (*rNewTangent)(2,4) = 0.;
            (*rNewTangent)(3,4) = 0.;
            (*rNewTangent)(4,4) = (*rNewTangent)(3,3);
            (*rNewTangent)(5,4) = 0.;

            (*rNewTangent)(0,5) = 0.;
            (*rNewTangent)(1,5) = 0.;
            (*rNewTangent)(2,5) = 0.;
            (*rNewTangent)(3,5) = 0.;
            (*rNewTangent)(4,5) = 0.;
            (*rNewTangent)(5,5) = (*rNewTangent)(3,3);
        }

        // static data is unchanged
        return;
    }

    //plastic loading
    H  = GetHardeningModulus(oldStaticData.GetEquivalentPlasticStrain(), d_H);
    H2 = H;

    int i=0;
    for (;i<100;i++)
    {
        g  = yield_condition - (2.*mu*delta_gamma + sqrt_2div3 * (H2-H));
        if (std::abs(g)<tolerance*sigma_y)
        {
            break;
        }
        dg = -2.* mu * (1.+(d_H+d_sigma)/(3.*mu));
        delta_gamma -=g/dg;
        epsilon_p_eq2 = oldStaticData.GetEquivalentPlasticStrain() + sqrt_2div3 * delta_gamma;

        H2  = GetHardeningModulus(epsilon_p_eq2,d_H);

        sigma_y = GetYieldStrength(epsilon_p_eq2,d_sigma);

        yield_condition = norm_dev - sqrt_2div3 * sigma_y;
    }

    if (i==100)
    {
        std::cout << "yield condition " << yield_condition << " delta_gamma " << delta_gamma << "\n";
        std::cout << "epsilon_p_eq " << oldStaticData.GetEquivalentPlasticStrain();
        std::cout << "total strain " << total_strain[0] << " " << total_strain[1] << " " << total_strain[2] << " " << total_strain[3] << " " << total_strain[4] << " " << total_strain[5] << " " <<"\n";
        std::cout << "plastic strain " << plastic_strain[0] << " " << plastic_strain[1] << " " << plastic_strain[2] << " " << plastic_strain[3] << " " << plastic_strain[4] << " " << plastic_strain[5] << " " <<"\n";
        std::cout << "back stress" << back_stress[0] << " " << back_stress[1] << " " << back_stress[2] << " " << back_stress[3] << " " << back_stress[4] << " " << back_stress[5] << " " <<"\n";
        std::cout << "[NuTo::MisesPlasticityEngineeringStress::ReturnMapping3D] No convergence after 100 steps, check the source code." << "\n";

        throw NuTo::Constitutive::DidNotConverge();

    }

    /* derivative of yield surface */
    df_dsigma[0]  = xi_trial[0]/norm_dev;
    df_dsigma[1]  = xi_trial[1]/norm_dev;
    df_dsigma[2]  = xi_trial[2]/norm_dev;
    df_dsigma[3]  = xi_trial[3]/norm_dev;
    df_dsigma[4]  = xi_trial[4]/norm_dev;
    df_dsigma[5]  = xi_trial[5]/norm_dev;

    //update static data
    if (rNewStaticData!=0)
    {
    	//update equivalent plastic strain
        double newEquivalentPlasticStrain = oldStaticData.GetEquivalentPlasticStrain() + sqrt_2div3 * delta_gamma;
    	rNewStaticData->SetEquivalentPlasticStrain(newEquivalentPlasticStrain);

        //update backstress
        factor = sqrt_2div3 * (H2-H);
        EngineeringStress<3> newBackStress;
        newBackStress[0] = back_stress[1] + factor*df_dsigma[0];
        newBackStress[1] = back_stress[1] + factor*df_dsigma[1];
        newBackStress[2] = back_stress[2] + factor*df_dsigma[2];
        newBackStress[3] = back_stress[3] + factor*df_dsigma[3];
        newBackStress[4] = back_stress[4] + factor*df_dsigma[4];
        newBackStress[5] = back_stress[5] + factor*df_dsigma[5];
        rNewStaticData->SetBackStress(newBackStress);

        //update plastic_strain
        EngineeringStrain<3> newPlasticStrain;
        newPlasticStrain[0] = plastic_strain[0] + delta_gamma*df_dsigma[0];
        newPlasticStrain[1] = plastic_strain[1] + delta_gamma*df_dsigma[1];
        newPlasticStrain[2] = plastic_strain[2] + delta_gamma*df_dsigma[2];
        newPlasticStrain[3] = plastic_strain[3] + 2.*delta_gamma*df_dsigma[3];  /* gamma */
        newPlasticStrain[4] = plastic_strain[4] + 2.*delta_gamma*df_dsigma[4];  /* gamma */
        newPlasticStrain[5] = plastic_strain[5] + 2.*delta_gamma*df_dsigma[5];  /* gamma */
        rNewStaticData->SetPlasticStrain(newPlasticStrain);
    }

    //update stress
    if (rNewStress!=0)
    {
		factor  = 2.*mu*delta_gamma;
		factor2 = bulk_modulus*trace_epsilon;
		(*rNewStress)[0] = factor2+sigma_trial[0]-factor*df_dsigma[0];
		(*rNewStress)[1] = factor2+sigma_trial[1]-factor*df_dsigma[1];
		(*rNewStress)[2] = factor2+sigma_trial[2]-factor*df_dsigma[2];
		(*rNewStress)[3] =         sigma_trial[3]-factor*df_dsigma[3];
		(*rNewStress)[4] =         sigma_trial[4]-factor*df_dsigma[4];
		(*rNewStress)[5] =         sigma_trial[5]-factor*df_dsigma[5];
    }

    //update stiffness
    if (rNewTangent!=0)
    {
        double theta     = 1.-2.*mu*delta_gamma/norm_dev;
        double theta_bar = 1./(1.+(d_sigma+d_H)/(3.*mu))-(1.-theta);
        factor    = 2.*mu*theta;
        double factor_div3  = factor/3.;
        double factor_mul2div3  = 2.*factor_div3;
        double factor2   = -2.*mu*theta_bar;
        double factor3   = bulk_modulus;

        (*rNewTangent)(0,0) =   (factor3 + factor_mul2div3 + factor2*df_dsigma[0]*df_dsigma[0]);
        (*rNewTangent)(1,0) =   (factor3 - factor_div3     + factor2*df_dsigma[0]*df_dsigma[1]);
        (*rNewTangent)(2,0) =   (factor3 - factor_div3     + factor2*df_dsigma[0]*df_dsigma[2]);
        (*rNewTangent)(3,0) =   (							 factor2*df_dsigma[0]*df_dsigma[3]);
        (*rNewTangent)(4,0) =   ( 						 factor2*df_dsigma[0]*df_dsigma[4]);
        (*rNewTangent)(5,0) =   ( 						 factor2*df_dsigma[0]*df_dsigma[5]);
        (*rNewTangent)(0,1) =   (*rNewTangent)(1,0);
        (*rNewTangent)(1,1) =   (factor3 + factor_mul2div3 + factor2*df_dsigma[1]*df_dsigma[1]);
        (*rNewTangent)(2,1) =   (factor3 - factor_div3     + factor2*df_dsigma[1]*df_dsigma[2]);
        (*rNewTangent)(3,1) =   (						 factor2*df_dsigma[1]*df_dsigma[3]);
        (*rNewTangent)(4,1) =   (						 factor2*df_dsigma[1]*df_dsigma[4]);
        (*rNewTangent)(5,1) =   (						 factor2*df_dsigma[1]*df_dsigma[5]);
        (*rNewTangent)(0,2) =   (*rNewTangent)(2,0);
        (*rNewTangent)(1,2) =   (*rNewTangent)(2,1);
        (*rNewTangent)(2,2) =   (factor3 + factor_mul2div3 + factor2*df_dsigma[2]*df_dsigma[2]);
        (*rNewTangent)(3,2) =   (						 factor2*df_dsigma[2]*df_dsigma[3]);
        (*rNewTangent)(4,2) =   (						 factor2*df_dsigma[2]*df_dsigma[4]);
        (*rNewTangent)(5,2) =   (						 factor2*df_dsigma[2]*df_dsigma[5]);
        (*rNewTangent)(0,3) =   (*rNewTangent)(3,0);
        (*rNewTangent)(1,3) =   (*rNewTangent)(3,1);
        (*rNewTangent)(2,3) =   (*rNewTangent)(3,2);
        (*rNewTangent)(3,3) =   (	 0.5*factor 		 +factor2*df_dsigma[3]*df_dsigma[3]);
        (*rNewTangent)(4,3) =   (					      factor2*df_dsigma[3]*df_dsigma[4]);
        (*rNewTangent)(5,3) =   (					      factor2*df_dsigma[3]*df_dsigma[5]);
        (*rNewTangent)(0,4) =   (*rNewTangent)(4,0);
        (*rNewTangent)(1,4) =   (*rNewTangent)(4,1);
        (*rNewTangent)(2,4) =   (*rNewTangent)(4,2);
        (*rNewTangent)(3,4) =   (*rNewTangent)(4,3);
        (*rNewTangent)(4,4) =   (	 0.5*factor 		 +factor2*df_dsigma[4]*df_dsigma[4]);
        (*rNewTangent)(5,4) =   (					      factor2*df_dsigma[4]*df_dsigma[5]);
        (*rNewTangent)(0,5) =   (*rNewTangent)(5,0);
        (*rNewTangent)(1,5) =   (*rNewTangent)(5,1);
        (*rNewTangent)(2,5) =   (*rNewTangent)(5,2);
        (*rNewTangent)(3,5) =   (*rNewTangent)(5,3);
        (*rNewTangent)(4,5) =   (*rNewTangent)(5,4);
        (*rNewTangent)(5,5) =   (	 0.5*factor 		 +factor2*df_dsigma[5]*df_dsigma[5]);
    }
}


double NuTo::MisesPlasticityEngineeringStress::GetYieldStrength(double rEpsilonPEq, double& rDSigmaDEpsilonP) const
{
    assert(mSigma.size() > 0);
    std::vector<std::pair<double, double>>::const_iterator it(mSigma.begin());
    while (rEpsilonPEq >= it->first)
    {
        it++;
        if (it == mSigma.end())
        {
            // the maximum is reached, afterwards the yield strength remains constant
        	rDSigmaDEpsilonP = 0.;
            return (it-1)->second;
        }
    }

  	rDSigmaDEpsilonP = (it->second - (it-1)->second) / (it->first - (it-1)->first);
   	return rDSigmaDEpsilonP * (rEpsilonPEq - (it-1)->first) + (it-1)->second;
}


double NuTo::MisesPlasticityEngineeringStress::GetHardeningModulus(double rEpsilonPEq, double& rDHDEpsilonP) const
{
    assert(mH.size() > 0);
	std::vector<std::pair<double, double>>::const_iterator it(mH.begin());
    double H(0);
    if (mH.size() == 1)
    {
        rDHDEpsilonP = it->second;
    	H = rEpsilonPEq * rDHDEpsilonP;
        return H;
    }
    else
    {
		do
		{
			it++;
			if (it == mH.end())
			{
				// the maximum is reached, afterwards use a constant slope
				it--;
				rDHDEpsilonP = it->second;
				H += (rEpsilonPEq - it->first) * (it->second);
				return H;
			}
			if (rEpsilonPEq > it->first)
			{
				H += (it->first - (it-1)->first) * (it-1)->second;
			}
			else
			{
				it--;
				rDHDEpsilonP = it->second;
				H += (rEpsilonPEq - it->first) * (it->second);
				return H;
			}
		}
		while(true);
    }
}


double NuTo::MisesPlasticityEngineeringStress::GetParameterDouble(eConstitutiveParameter rIdentifier) const
{
    switch(rIdentifier)
    {
    case eConstitutiveParameter::INITIAL_HARDENING_MODULUS:
        if(mH.size()==0)
            throw MechanicsException(__PRETTY_FUNCTION__, "Size of hardening modulus vector is zero.");
        return mH[0].second;
    case eConstitutiveParameter::INITIAL_YIELD_STRENGTH:
        if(mSigma.size()==0)
            throw MechanicsException(__PRETTY_FUNCTION__, "Size of yield strength vector is zero.");
        return mSigma[0].second;
    case eConstitutiveParameter::POISSONS_RATIO:
        return this->mNu;
    case eConstitutiveParameter::YOUNGS_MODULUS:
        return this->mE;
    case eConstitutiveParameter::DENSITY:
    	return this->mRho;
    default:
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Constitutive law does not have the requested variable");
    }
    }
}


void NuTo::MisesPlasticityEngineeringStress::SetParameterDouble(eConstitutiveParameter rIdentifier, double rValue)
{

    switch(rIdentifier)
    {
    case eConstitutiveParameter::INITIAL_HARDENING_MODULUS:
    {
        if(mH.size()==0)
            throw MechanicsException(__PRETTY_FUNCTION__, "Size of hardening modulus vector is zero.");
        if (rValue<0)
            throw MechanicsException(__PRETTY_FUNCTION__, "Initial hardening modulus must not be negative.");
        mH[0].second = rValue;
        break;
    }
    case eConstitutiveParameter::INITIAL_YIELD_STRENGTH:
    {
        if(mSigma.size()==0)
            throw MechanicsException(__PRETTY_FUNCTION__, "Size of yield strength vector is zero.");
        if (rValue<=0)
            throw MechanicsException(__PRETTY_FUNCTION__, "Initial yield strength has to be positive.");
        mSigma[0].second = rValue;
        break;
    }
    case eConstitutiveParameter::POISSONS_RATIO:
    {
        this->mNu = rValue;
        break;
    }
    case eConstitutiveParameter::YOUNGS_MODULUS:
    {
        this->mE = rValue;
        break;
    }
    case eConstitutiveParameter::DENSITY:
    {
        this->mRho = rValue;
        break;
    }
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Constitutive law does not have the requested variable.");
    }

    this->SetParametersValid();
}


Eigen::MatrixXd NuTo::MisesPlasticityEngineeringStress::GetYieldStrength() const
{
	Eigen::MatrixXd returnMatrix(mSigma.size(),2);
	for (unsigned int count=0; count<mSigma.size(); count++)
	{
		returnMatrix(count,0) = mSigma[count].first;
		returnMatrix(count,1) = mSigma[count].second;
	}
	return returnMatrix;
}


void NuTo::MisesPlasticityEngineeringStress::AddYieldStrength(double rEpsilon, double rSigma)
{
	if (rEpsilon<=0)
		throw MechanicsException(__PRETTY_FUNCTION__, "Equivalente strain has to be positive.");
	if (rSigma<=0)
		throw MechanicsException(__PRETTY_FUNCTION__, "Yield strength has to be positive.");
	std::vector<std::pair<double, double> >::iterator it;
	for (it = mSigma.begin(); it!=mSigma.end(); it++)
	{
		if (it->first > rEpsilon)
		{
			if (it != mSigma.begin())
			{
				if ((it-1)->second > (it->second))
					throw MechanicsException(__PRETTY_FUNCTION__, "The yield strength can only increase for "
                            "increasing epsilon equivalent.");
			}
			break;
		}
	}
	mSigma.insert(it,1, std::pair<double,double>(rEpsilon, rSigma));
    this->SetParametersValid();
}


Eigen::MatrixXd NuTo::MisesPlasticityEngineeringStress::GetHardeningModulus() const
{
	Eigen::MatrixXd returnMatrix(mH.size(),2);
	for (unsigned int count=0; count<mH.size(); count++)
	{
		returnMatrix(count,0) = mH[count].first;
		returnMatrix(count,1) = mH[count].second;
	}
	return returnMatrix;
}


void NuTo::MisesPlasticityEngineeringStress::AddHardeningModulus(double rEpsilon, double rH)
{
	if (rEpsilon<=0)
		throw MechanicsException(__PRETTY_FUNCTION__, "Equivalente strain has to be positive.");
	if (rH<0)
		throw MechanicsException(__PRETTY_FUNCTION__, "Hardening modul must not be negative.");
	std::vector<std::pair<double, double> >::iterator it;
	for (it = mH.begin(); it!=mH.end(); it++)
	{
		if (it->first>rEpsilon)
		{
			break;
		}
	}
	mH.insert(it, 1, std::pair<double,double>(rEpsilon, rH));
    this->SetParametersValid();
}


NuTo::Constitutive::eConstitutiveType NuTo::MisesPlasticityEngineeringStress::GetType() const
{
    return Constitutive::eConstitutiveType::MISES_PLASTICITY_ENGINEERING_STRESS;
}


void NuTo::MisesPlasticityEngineeringStress::CheckYieldStrength(std::vector<std::pair<double, double> > rSigma) const
{
	if (rSigma.size()==0)
    {
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "At least an initial yield strength is required.");
    }
	if (rSigma[0].second <= 0.0)
    {
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "The initial yield strength must be a positive value.");
    }
	rSigma[0].first = 0.;

	for (unsigned int count=1; count<rSigma.size(); count++)
	{
		if (rSigma[count-1].first>=rSigma[count].first)
	        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "For multilinear plasticity, the epsilon should always increase.");
		if (rSigma[count-1].second>rSigma[count].second)
	        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "For multilinear plasticity, the yield strength should always increase.");
	}
}


void NuTo::MisesPlasticityEngineeringStress::CheckHardeningModulus(std::vector<std::pair<double, double> > rH) const
{
	if (rH.size()==0)
    {
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "At least an initial hardening modulus is required.");
    }
	if (rH[0].second < 0.0)
    {
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "The initial hardening modulus must not be a negative value.");
    }
	rH[0].first = 0.;

	for (unsigned int count=1; count<rH.size(); count++)
	{
		if (rH[count-1].first>=rH[count].first)
	        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "For multilinear plasticity, the epsilon should always increase.");
		if (rH[count].second<0)
	        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "For multilinear plasticity, the hardening modulus should always be positive.");
	}
}


void NuTo::MisesPlasticityEngineeringStress::Info(unsigned short rVerboseLevel, Logger& rLogger) const
{
    this->ConstitutiveBase::Info(rVerboseLevel, rLogger);
    rLogger << "    Young's modulus: " << this->mE << "\n";
    rLogger << "    Poisson's ratio: " << this->mNu << "\n";
    rLogger << "    multilinear yield strength: (interval,epsilon,sigma)" << "\n";
    for (unsigned int count=0; count<mSigma.size(); count++)
        rLogger << "       " << count<< " : " << this->mSigma[count].first << "    " << this->mSigma[count].second << "\n";

    rLogger << "    multilinear hardening modulus: (interval,epsilon,H')" << "\n";
    for (unsigned int count=0; count<mH.size(); count++)
    	rLogger << "       " << count<< " : " << this->mH[count].first << "    " << this->mH[count].second << "\n";
}

// check parameters
void NuTo::MisesPlasticityEngineeringStress::CheckParameters()const
{
    ConstitutiveBase::CheckParameterDouble(eConstitutiveParameter::YOUNGS_MODULUS, mE);
    ConstitutiveBase::CheckParameterDouble(eConstitutiveParameter::POISSONS_RATIO, mNu);
    ConstitutiveBase::CheckParameterDouble(eConstitutiveParameter::DENSITY, mRho);

    this->CheckYieldStrength(this->mSigma);
    this->CheckHardeningModulus(this->mH);
}
