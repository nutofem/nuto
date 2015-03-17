// $Id: GradientDamageEngineeringStress.cpp 612 2012-08-13 07:31:23Z unger3 $
// GradientDamageEngineeringStress.cpp
// created Apr 26, 2010 by Joerg F. Unger
#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include <eigen3/Eigen/LU>

#include "nuto/base/Logger.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveStaticDataBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveTangentNonlocal.h"
#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataGradientDamage1D.h"
#include "nuto/mechanics/constitutive/mechanics/Damage.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient1D.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient2D.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient3D.h"
#include "nuto/mechanics/constitutive/mechanics/GradientDamageEngineeringStress.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress1D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress2D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress3D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain1D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain2D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain3D.h"
#include "nuto/mechanics/constitutive/mechanics/LocalEqStrain.h"
#include "nuto/mechanics/constitutive/mechanics/NonlocalEqStrain.h"
#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/sections/SectionBase.h"
#include "nuto/mechanics/sections/SectionEnum.h"

#define MAX_OMEGA 0.99
//#define ENABLE_DEBUG

NuTo::GradientDamageEngineeringStress::GradientDamageEngineeringStress() : ConstitutiveBase()
{
    mRho = 0.;
    mE = 0.;
    mNu = 0.;
    mNonlocalRadius = 0.;
    mNonlocalRadiusParameter = 0.;
    mThermalExpansionCoefficient = 0.;
    mTensileStrength = 0.;
    mFractureEnergy = 0.;

    mDamageLawType = Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING;


    SetParametersValid();
#ifdef ENABLE_DEBUG
    std::cout << "NuTo::GradientDamageEngineeringStress::GradientDamageEngineeringStress debug active" << std::endl;
#else
    std::cout << "NuTo::GradientDamageEngineeringStress::GradientDamageEngineeringStress debug inactive" << std::endl;
#endif

}

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void NuTo::GradientDamageEngineeringStress::serialize(Archive & ar, const unsigned int version)
    {
#ifdef DEBUG_SERIALIZATION
       std::cout << "start serialize GradientDamageEngineeringStress" << "\n";
#endif
       ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveBase)
          & BOOST_SERIALIZATION_NVP(mRho)
          & BOOST_SERIALIZATION_NVP(mE)
          & BOOST_SERIALIZATION_NVP(mNu)
          & BOOST_SERIALIZATION_NVP(mNonlocalRadius)
          & BOOST_SERIALIZATION_NVP(mThermalExpansionCoefficient)
          & BOOST_SERIALIZATION_NVP(mDamageLawType)
          & BOOST_SERIALIZATION_NVP(mDamageLawParameters);
#ifdef DEBUG_SERIALIZATION
       std::cout << "finish serialize GradientDamageEngineeringStress" << "\n";
#endif
    }
    BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::GradientDamageEngineeringStress)
#endif // ENABLE_SERIALIZATION

//! @brief ... evaluate the constitutive relation in 1D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
//! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
NuTo::Error::eError NuTo::GradientDamageEngineeringStress::Evaluate1D(ElementBase* rElement, int rIp,
        const std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*>& rConstitutiveInput,
        std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>& rConstitutiveOutput)
{
        // get section information determining which input on the constitutive level should be used
        const SectionBase* section(rElement->GetSection());


        // check if parameters are valid
        if (this->mParametersValid == false)
        {
            //throw an exception giving information related to the wrong parameter
            CheckParameters();
        }

        if (section->GetType()!=Section::TRUSS)
            throw MechanicsException("[NuTo::GradientDamageEngineeringStress::Evaluate1D] only truss sections are implemented.");

        // calculate engineering strain
        EngineeringStrain1D strain1D;
        if(rConstitutiveInput.find(NuTo::Constitutive::Input::DEFORMATION_GRADIENT_1D)==rConstitutiveInput.end())
            throw MechanicsException("[NuTo::GradientDamageEngineeringStress::Evaluate1D] deformation gradient 1d needed to evaluate engineering strain1d.");
        const DeformationGradient1D& deformationGradient(rConstitutiveInput.find(NuTo::Constitutive::Input::DEFORMATION_GRADIENT_1D)->second->GetDeformationGradient1D());
        deformationGradient.GetEngineeringStrain(strain1D);

        //Get previous ip_data
        ConstitutiveStaticDataGradientDamage1D *oldStaticData = (rElement->GetStaticData(rIp))->AsGradientDamage1D();

        // subtract thermal strain
        if (section->GetInputConstitutiveIsTemperature())
        {
            std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*>::const_iterator itInput(rConstitutiveInput.find(NuTo::Constitutive::Input::TEMPERATURE));
            if (itInput==rConstitutiveInput.end())
                throw MechanicsException("[NuTo::GradientDamageEngineeringStress::Evaluate1D] temperature needed to evaluate thermal engineering strain2d.");
            double temperature(itInput->second->GetTemperature());
            double deltaStrain(mThermalExpansionCoefficient * temperature);
            strain1D[0] -= deltaStrain;
            throw MechanicsException("[NuTo::GradientDamageEngineeringStress::Evaluate1D] add temperature history.");
/*            double prevTemperature(oldStaticData->mPrevTemperature);
            double deltaStrain(mThermalExpansionCoefficient * prevTemperature);
            prevEngineeringStrain3D.mEngineeringStrain[0] -= deltaStrain;
            prevEngineeringStrain3D.mEngineeringStrain[1] -= deltaStrain;
            prevEngineeringStrain3D.mEngineeringStrain[2] -= deltaStrain;
*/
        }

        // nonlocal eq strain
        if(rConstitutiveInput.find(NuTo::Constitutive::Input::NONLOCAL_EQ_STRAIN)==rConstitutiveInput.end())
            throw MechanicsException("[NuTo::GradientDamageEngineeringStress::Evaluate1D] nonlocal eq strain needed to evaluate engineering stress 1d.");
        const double nonlocalEqStrain = rConstitutiveInput.find(NuTo::Constitutive::Input::NONLOCAL_EQ_STRAIN)->second->GetNonlocalEqStrain().GetValue(0);

        // last converged kappa
        const double kappaLastConverged = oldStaticData->mKappa;

        // kappa
        const double kappa = std::max(nonlocalEqStrain, kappaLastConverged);

        //calculate damage
        double omega = CalculateDamage(kappa);

        // calculate stress
        EngineeringStress1D stress1D;
        stress1D(0) = mE * strain1D(0);

        // calculate local eq strain and its derivative
        NuTo::FullVector<double, 1> localEqStrainDerivative;
        LocalEqStrain localEqStrain;

        localEqStrain[0] = strain1D[0];
        localEqStrainDerivative[0] = 1.;





        //set this to true, if update is in the map, perform the update after all other outputs have been calculated
        bool performUpdateAtEnd(false);

        for (std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>::iterator itOutput = rConstitutiveOutput.begin();
                itOutput != rConstitutiveOutput.end(); itOutput++)
        {
            switch(itOutput->first)
            {
            case NuTo::Constitutive::Output::ENGINEERING_STRESS_1D:
            {
                EngineeringStress1D& engineeringStress1D = itOutput->second->GetEngineeringStress1D();
                engineeringStress1D = (1.-omega)*stress1D;
            }
            break;
            case NuTo::Constitutive::Output::ENGINEERING_STRESS_3D:
            {
                //this is for the visualize routines
                EngineeringStress3D& engineeringStress3D(itOutput->second->GetEngineeringStress3D());

                engineeringStress3D[0] = (1.-omega)*stress1D[0];
                engineeringStress3D[1] = 0.;
                engineeringStress3D[2] = 0.;
                engineeringStress3D[3] = 0.;
                engineeringStress3D[4] = 0.;
                engineeringStress3D[5] = 0.;
            }
            break;
            case NuTo::Constitutive::Output::ENGINEERING_STRAIN_3D:
            {
                EngineeringStrain3D& engineeringStrain3D(itOutput->second->GetEngineeringStrain3D());
                engineeringStrain3D[0] = strain1D[0];
                engineeringStrain3D[1] = -mNu*strain1D[0];
                engineeringStrain3D[2] = -mNu*strain1D[0];
                engineeringStrain3D[3] = 0.;
                engineeringStrain3D[4] = 0.;
                engineeringStrain3D[5] = 0.;
            }
            break;
            case NuTo::Constitutive::Output::LOCAL_EQ_STRAIN:
            {
                LocalEqStrain& localEqStrainOut = itOutput->second->GetLocalEqStrain();
                localEqStrainOut = localEqStrain;
            }
            break;
            case NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_1D:
            {
                ConstitutiveTangentLocal<1,1>& tangent(itOutput->second->AsConstitutiveTangentLocal_1x1());
                tangent.SetSymmetry(true);
                tangent(0,0) = (1.-omega) * mE;
            }
            break;
            case NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_NONLOCAL_EQ_STRAIN_1D:
            {
                ConstitutiveTangentLocal<1,1>& tangent(itOutput->second->AsConstitutiveTangentLocal_1x1());
                tangent.SetSymmetry(true);
                if (nonlocalEqStrain > kappaLastConverged)
                {
                    // loading

                    tangent(0,0) = - CalculateDerivativeDamage(nonlocalEqStrain) * stress1D[0];

                } else {
                    // unloading
                    tangent(0,0) = 0.;
                }
            }
            break;
            case NuTo::Constitutive::Output::D_LOCAL_EQ_STRAIN_D_STRAIN_1D:
            {
                ConstitutiveTangentLocal<1,1>& tangent(itOutput->second->AsConstitutiveTangentLocal_1x1());
                tangent(0,0) = localEqStrainDerivative[0];
            }
            break;
            case NuTo::Constitutive::Output::D_LOCAL_EQ_STRAIN_XI_D_STRAIN_1D:
            {
                ConstitutiveTangentLocal<1,1>& tangent(itOutput->second->AsConstitutiveTangentLocal_1x1());
                tangent.SetSymmetry(true);

                double xi = mNonlocalRadius;
                double dXi = 0;

                if (mNonlocalRadiusParameter != 0.)
                {
                    double e_xi = mTensileStrength / mE * mNonlocalRadiusParameter;
                    if (localEqStrain[0] <= e_xi)
                    {
                        double c0 = 0.01;
                        xi = c0 + (mNonlocalRadius - c0) * (localEqStrain[0] / e_xi);
                        dXi = (mNonlocalRadius - c0)/ e_xi;
                    }
                }

                double factor = 1./xi - (localEqStrain[0] - nonlocalEqStrain)/(xi*xi)*dXi;

                tangent(0,0) = factor * localEqStrainDerivative[0];
            }
            break;
            case NuTo::Constitutive::Output::NONLOCAL_PARAMETER_XI:
            {
                ConstitutiveTangentLocal<1,1>& variableNonlocalRadius(itOutput->second->AsConstitutiveTangentLocal_1x1());
                variableNonlocalRadius[0] = mNonlocalRadius;
                if (mNonlocalRadiusParameter != 0.)
                {
                    double e_xi = mTensileStrength / mE * mNonlocalRadiusParameter;
                    if (localEqStrain[0] <= e_xi)
                    {
                        double c0 = 0.01;
                        variableNonlocalRadius[0] = c0 + (mNonlocalRadius - c0) * (localEqStrain[0] / e_xi);
                    }
                }
            }
            break;
            case NuTo::Constitutive::Output::DAMAGE:
            {
                itOutput->second->GetDamage().SetDamage(omega);
            }
            break;
            case NuTo::Constitutive::Output::UPDATE_TMP_STATIC_DATA:
            {
                   throw MechanicsException("[NuTo::GradientDamageEngineeringStress::Evaluate1D] tmp_static_data has to be updated without any other outputs, call it separately.");
            }
            break;
            case NuTo::Constitutive::Output::UPDATE_STATIC_DATA:
            {
                performUpdateAtEnd = true;
            }
            break;
            default:
                throw MechanicsException(std::string("[NuTo::GradientDamageEngineeringStress::Evaluate1D] output object ") +
                        NuTo::Constitutive::OutputToString(itOutput->first) +
                        std::string(" could not be calculated, check the allocated material law and the section behavior."));
            }
        }

        //update history variables
        if (performUpdateAtEnd)
        {
            oldStaticData->mKappa = kappa;
        }
        return Error::SUCCESSFUL;
}

//! @brief ... evaluate the constitutive relation in 2D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
//! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
NuTo::Error::eError NuTo::GradientDamageEngineeringStress::Evaluate2D(ElementBase* rElement, int rIp,
        const std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*>& rConstitutiveInput,
        std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>& rConstitutiveOutput)
{
    throw MechanicsException("[NuTo::GradientDamageEngineeringStress::Evaluate2D] not implemented for 2D.");
    return Error::SUCCESSFUL;
}

//! @brief ... evaluate the constitutive relation in 2D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
//! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
NuTo::Error::eError NuTo::GradientDamageEngineeringStress::Evaluate3D(ElementBase* rElement, int rIp,
        const std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*>& rConstitutiveInput,
        std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>& rConstitutiveOutput)
{
    throw NuTo::MechanicsException("[NuTo::GradientDamageEngineeringStress::Evaluate3D] not implemented for 3D.");
}

//! @brief ... create new static data object for an integration point
//! @return ... pointer to static data object
NuTo::ConstitutiveStaticDataBase* NuTo::GradientDamageEngineeringStress::AllocateStaticDataEngineeringStress_EngineeringStrain1D(
        const ElementBase* rElement) const
{
    return new ConstitutiveStaticDataGradientDamage1D();
}


//! @brief ... create new static data object for an integration point
//! @return ... pointer to static data object
NuTo::ConstitutiveStaticDataBase* NuTo::GradientDamageEngineeringStress::AllocateStaticDataEngineeringStress_EngineeringStrain2D(
        const ElementBase* rElement) const
{
    throw MechanicsException("[NuTo::GradientDamageEngineeringStress::AllocateStaticDataEngineeringStress_EngineeringStrain2D] Nonlocal damage plasticity model not implemented for 2D.");
}


//! @brief ... create new static data object for an integration point
//! @return ... pointer to static data object
NuTo::ConstitutiveStaticDataBase* NuTo::GradientDamageEngineeringStress::AllocateStaticDataEngineeringStress_EngineeringStrain3D(
        const ElementBase* rElement) const
{
    throw MechanicsException("[NuTo::GradientDamageEngineeringStress::AllocateStaticDataEngineeringStress_EngineeringStrain3D] Nonlocal damage plasticity model not implemented for 3D.");
}


// calculate coefficients of the material matrix
void NuTo::GradientDamageEngineeringStress::CalculateCoefficients3D(double& C11, double& C12, double& C44) const
{
    double factor = this->mE/((1.0 + this->mNu) * (1.0 - 2.0 * this->mNu));
    C11 = factor * (1.0 - this->mNu);
    C12 = factor * this->mNu;
    C44 = this->mE/(2*(1.0 + this->mNu));
}

// parameters /////////////////////////////////////////////////////////////
//! @brief ... get density
//! @return ... density
double NuTo::GradientDamageEngineeringStress::GetDensity() const
{
    return this->mRho;
}

//! @brief ... set density
//! @param rRho ... density
void NuTo::GradientDamageEngineeringStress::SetDensity(double rRho)
{
    this->CheckDensity(rRho);
    this->mRho = rRho;
    this->SetParametersValid();
}

//! @brief ... get Young's modulus
//! @return ... Young's modulus
double NuTo::GradientDamageEngineeringStress::GetYoungsModulus() const
{
    return mE;
}


//! @brief ... set Young's modulus
//! @param rE ... Young's modulus
void NuTo::GradientDamageEngineeringStress::SetYoungsModulus(double rE)
{
    this->CheckYoungsModulus(rE);
    this->mE = rE;
    this->SetParametersValid();
}


//! @brief ... get Poisson's ratio
//! @return ... Poisson's ratio
double NuTo::GradientDamageEngineeringStress::GetPoissonsRatio() const
{
    return mNu;
}

//! @brief ... set Poisson's ratio
//! @param rNu ... Poisson's ratio
void NuTo::GradientDamageEngineeringStress::SetPoissonsRatio(double rNu)
{
    this->CheckPoissonsRatio(rNu);
    this->mNu = rNu;
    this->SetParametersValid();
}

//! @brief ... get nonlocal radius
//! @return ... nonlocal radius
double NuTo::GradientDamageEngineeringStress::GetNonlocalRadius() const
{
    return mNonlocalRadius;
}

//! @brief ... set nonlocal radius
//! @param rRadius...  nonlocal radius
void NuTo::GradientDamageEngineeringStress::SetNonlocalRadius(double rNonlocalRadius)
{
    this->CheckNonlocalRadius(rNonlocalRadius);
    this->mNonlocalRadius = rNonlocalRadius;
    this->SetParametersValid();
}

//! @brief ... get nonlocal radius parameter
//! @return ... nonlocal radius parameter
double NuTo::GradientDamageEngineeringStress::GetNonlocalRadiusParameter() const
{
    return mNonlocalRadiusParameter;
}

//! @brief ... set nonlocal radius parameter
//! @param rRadius ...  nonlocal radius parameter
void NuTo::GradientDamageEngineeringStress::SetNonlocalRadiusParameter(double rRadiusParameter)
{
    this->CheckNonlocalRadiusParameter(rRadiusParameter);
    mNonlocalRadiusParameter = rRadiusParameter;
    this->SetParametersValid();
}

//! @brief ... get thermal expansion coefficient
//! @return ... thermal expansion coefficient
double NuTo::GradientDamageEngineeringStress::GetThermalExpansionCoefficient() const
{
    return mThermalExpansionCoefficient;
}

//! @brief ... set thermal expansion coefficient
//! @param rNu ... thermal expansion coefficient
void NuTo::GradientDamageEngineeringStress::SetThermalExpansionCoefficient(double rAlpha)
{
    this->CheckThermalExpansionCoefficient(rAlpha);
    this->mThermalExpansionCoefficient = rAlpha;
    this->SetParametersValid();
}

//! @brief ... get fracture energy
//! @return ... fracture energy
double NuTo::GradientDamageEngineeringStress::GetFractureEnergy() const
{
    return mFractureEnergy;
}

//! @brief ... set fracture energy
//! @param rFractureEnergy... fracture energy
void NuTo::GradientDamageEngineeringStress::SetFractureEnergy(double rFractureEnergy)
{
    this->CheckFractureEnergy(rFractureEnergy);
    this->mFractureEnergy = rFractureEnergy;
    this->SetParametersValid();
}

//! @brief ... get tensile strength
//! @return ... tensile strength
double NuTo::GradientDamageEngineeringStress::GetTensileStrength() const
{
    return mTensileStrength;
}

//! @brief ... set tensile strength
//! @param rTensileStrength...  tensile strength
void NuTo::GradientDamageEngineeringStress::SetTensileStrength(double rTensileStrength)
{
    this->CheckTensileStrength(rTensileStrength);
    this->mTensileStrength = rTensileStrength;
    this->SetParametersValid();
}

//! @brief ... get damage law
//! @return ... damage law
NuTo::FullVector<double, Eigen::Dynamic> NuTo::GradientDamageEngineeringStress::GetDamageLaw() const
{
    NuTo::FullVector<double, Eigen::Dynamic> damageLaw(mDamageLawParameters.rows() + 1);

    damageLaw[0] = static_cast<double>(mDamageLawType);
    damageLaw.SetBlock(1,0, mDamageLawParameters);
    return damageLaw;
}

//! @brief ... set damage law parameters
//! @param rDamageLaw ... damage law parameters
void NuTo::GradientDamageEngineeringStress::SetDamageLaw(
        const NuTo::FullVector<double, Eigen::Dynamic> rDamageLaw)
{
    this->CheckDamageLaw(rDamageLaw);
    mDamageLawType = static_cast<int>(rDamageLaw[0]);
    int numDamageLawParameters = rDamageLaw.rows() - 1;
    if (numDamageLawParameters > 0)
        mDamageLawParameters = rDamageLaw.GetBlock(1,0,numDamageLawParameters,1);

    this->SetParametersValid();

}

///////////////////////////////////////////////////////////////////////////


//! @brief ... get type of constitutive relationship
//! @return ... type of constitutive relationship
//! @sa eConstitutiveType
NuTo::Constitutive::eConstitutiveType NuTo::GradientDamageEngineeringStress::GetType() const
{
    return NuTo::Constitutive::GRADIENT_DAMAGE_ENGINEERING_STRESS;
}


//! @brief ... check compatibility between element type and type of constitutive relationship
//! @param rElementType ... element type
//! @return ... <B>true</B> if the element is compatible with the constitutive relationship, <B>false</B> otherwise.
bool NuTo::GradientDamageEngineeringStress::CheckElementCompatibility(NuTo::Element::eElementType rElementType) const
{
    switch (rElementType)
    {
    case NuTo::Element::BRICK8N:
        return false;
    case NuTo::Element::PLANE2D3N:
        return false;
    case NuTo::Element::PLANE2D4N:
        return false;
    case NuTo::Element::PLANE2D6N:
        return false;
    case NuTo::Element::TETRAHEDRON4N:
        return false;
    case NuTo::Element::TETRAHEDRON10N:
        return false;
    case NuTo::Element::TRUSS1D2N:
        return true;
    case NuTo::Element::TRUSS1D3N:
        return true;
    case NuTo::Element::TRUSS1D4NDISP3NX:
        return true;
    case NuTo::Element::BOUNDARYGRADIENTDAMAGE1D:
        return true;
    default:
        return false;
    }
}

//! @brief ... check if density is non negative
//! @param rE ... Young's modulus
void NuTo::GradientDamageEngineeringStress::CheckDensity(double rRho) const
{
    if (rRho < 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::GradientDamageEngineeringStress::CheckDensity] The density must not be negative.");
    }
}

//! @brief ... check if Young's modulus is positive
//! @param rE ... Young's modulus
void NuTo::GradientDamageEngineeringStress::CheckYoungsModulus(double rE) const
{
    if (rE <= 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::GradientDamageEngineeringStress::CheckYoungsModulus] The Young's modulus must be a positive value.");
    }
}

//! @brief ... check if Poisson's ratio is valid \f$ (-1.0 < \nu < 0.5) \f$
//! @param rE ... Poisson's ratio
void NuTo::GradientDamageEngineeringStress::CheckPoissonsRatio(double rNu) const
{
    if (rNu <= -1.0)
    {
        throw NuTo::MechanicsException("[NuTo::GradientDamageEngineeringStress::CheckPoissonsRatio] Poisson's ratio must be greater or equal to -1.0.");
    }
    if (rNu >= 0.5)
    {
        throw NuTo::MechanicsException("[NuTo::GradientDamageEngineeringStress::CheckPoissonsRatio] Poisson's ratio must be smaller or equal to 0.5.");
    }
}

//! @brief ... check if the nonlocal radius is positive
//! @param rRadius ... nonlocal radius
void NuTo::GradientDamageEngineeringStress::CheckNonlocalRadius(double rRadius) const
{
    if (rRadius <= 0)
    {
        throw NuTo::MechanicsException("[NuTo::GradientDamageEngineeringStress::CheckNonlocalRadius] Nonlocal radius must be positive.");
    }
}

//! @brief ... check if the nonlocal radius parameter is greater than 1
//! @param rRadius ... nonlocal radius
void NuTo::GradientDamageEngineeringStress::CheckNonlocalRadiusParameter(double rRadiusParameter) const
{
    if (rRadiusParameter <= 1)
    {
        throw NuTo::MechanicsException("[NuTo::GradientDamageEngineeringStress::CheckNonlocalRadius] Nonlocal radius parameter must be greater than 1.");
    }
}


//! @brief ... check thermal expansion coefficient
//! @param rAlpha ... thermal expansion coefficient
void NuTo::GradientDamageEngineeringStress::CheckThermalExpansionCoefficient(double rAlpha) const
{
}

//! @brief ... check if tensile strength is positive
//! @param rTensileStrength ... nonlocal radius
void NuTo::GradientDamageEngineeringStress::CheckTensileStrength(double rTensileStrength) const
{
    if (rTensileStrength <= 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::GradientDamageEngineeringStress::CheckTensileStrength] The tensile strength must be a positive value.");
    }
}

//! @brief ... check if fracture energy is positive
//! @param rFractureEnergy ... fracture energy
void NuTo::GradientDamageEngineeringStress::CheckFractureEnergy(double rFractureEnergy) const
{
    if (rFractureEnergy <= 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::NonlocalDamagePlasticityEngineeringStress::CheckFractureEnergy] The fracture energy must be a positive value.");
    }
}

//! @brief ... check damage law parameters
//! @param rDamageLawParameters ... damage law parameters
void NuTo::GradientDamageEngineeringStress::CheckDamageLaw(
        const NuTo::FullVector<double, Eigen::Dynamic>& rDamageLaw) const
{
    int damageLawType = static_cast<int>(rDamageLaw[0]);
    int numDamageLawParameters = rDamageLaw.rows() - 1;

    switch (damageLawType) {
    case Constitutive::eDamageLawType::ISOTROPIC_NO_SOFTENING:
        break;
    case Constitutive::eDamageLawType::ISOTROPIC_LINEAR_SOFTENING:
        break;
    case Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING:
        break;
    case Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING_RES_LOAD:
        break;
    case Constitutive::eDamageLawType::ISOTROPIC_CUBIC_HERMITE:
        break;
    case Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING_SMOOTH:
    {
        if (numDamageLawParameters != 3)
            throw NuTo::MechanicsException(
                    std::string("[NuTo::GradientDamageEngineeringStress::CheckDamageLaw] ") +
                    std::string("ISOTROPIC_EXPONENTIAL_SOFTENING_SMOOTH needs exactly 3 parameters."));

    }
        break;
    default:
    {
        throw NuTo::MechanicsException(
                std::string("[NuTo::GradientDamageEngineeringStress::CheckDamageLaw] ") +
                std::string("The required damage law is not implemented. "));
        break;
    }

    }

}

//! @brief ... print information about the object
//! @param rVerboseLevel ... verbosity of the information
//! @param rLogger stream for the output
void NuTo::GradientDamageEngineeringStress::Info(unsigned short rVerboseLevel, Logger& rLogger) const
{
    this->ConstitutiveBase::Info(rVerboseLevel, rLogger);
    rLogger << "    Young's modulus      : " << this->mE << "\n";
    rLogger << "    Poisson's ratio      : " << this->mNu << "\n";
    rLogger << "    nonlocal radius      : " << this->mNonlocalRadius << "\n";
    rLogger << "    damage law params    : " << this->mDamageLawParameters << "\n";
    rLogger << "    thermal expansion coeff : " << this->mThermalExpansionCoefficient << "\n";
}

// check parameters
void NuTo::GradientDamageEngineeringStress::CheckParameters()const
{
    this->CheckDensity(mRho);
    this->CheckYoungsModulus(mE);
    this->CheckPoissonsRatio(mNu);
    this->CheckNonlocalRadius(mNonlocalRadius);
    this->CheckTensileStrength(mTensileStrength);
    this->CheckFractureEnergy(mFractureEnergy);
    this->CheckDamageLaw(GetDamageLaw());
    this->CheckThermalExpansionCoefficient(mThermalExpansionCoefficient);
}


double NuTo::GradientDamageEngineeringStress::CalculateDamage(double rKappa)const
{
    double omega = 0;

    double e_0 = mTensileStrength / mE;
    double e_f = mFractureEnergy / mTensileStrength;
    double e_c = 2* e_f; // or something


    switch (mDamageLawType) {
        case Constitutive::eDamageLawType::ISOTROPIC_NO_SOFTENING:
        {
            if (rKappa > e_0)
            {
                omega = 1 - e_0 / rKappa;
            }
            break;
        }
        case Constitutive::eDamageLawType::ISOTROPIC_LINEAR_SOFTENING:
        {
            if (rKappa > e_0)
            {
                omega = e_c / rKappa * (rKappa - e_0) / (e_c - e_0);
                omega = std::min(omega, MAX_OMEGA);
            }
            break;
        }
        case Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING:
        {
            if (rKappa > e_0)
            {
                omega = 1 - e_0 / rKappa * exp( (e_0 - rKappa) / e_f);
            }
            break;
        }
        case Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING_RES_LOAD:
        {
            if (rKappa > e_0)
            {
                omega = 1 - e_0 / rKappa *(1 - MAX_OMEGA + MAX_OMEGA* exp( (e_0 - rKappa) / e_f));
            }
            break;
        }
        case Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING_SMOOTH:
        {
            double alpha = mDamageLawParameters[0];
            double beta = mDamageLawParameters[1];
            double gamma = mDamageLawParameters[2];

            double e_Dev_e0 = rKappa / e_0;
            double e0_Dev_e = e_0/rKappa;

            if (e_Dev_e0 <= alpha)
                break;

            if (e_Dev_e0 < beta)
            {
                // pre peak polynomial smoothing
                omega = 1. - e0_Dev_e*((1-alpha)*pow((e_Dev_e0-beta)/(beta-alpha),3.) +1 );
                break;
            }


            double term = 2. * e_f / (e_0*(gamma-beta));
            double A = -1. / (1+term);
            if (e_Dev_e0 < gamma)
            {
                // post peak polynomial smoothing
                omega = 1. - e0_Dev_e*( pow((e_Dev_e0-beta)/(gamma-beta),2.) * A  +1);
                break;
            }
            double e_s = e_0*gamma+e_f*log(-A*term);


            // else
            omega = 1. - e0_Dev_e * exp((e_s-rKappa)/e_f);

            break;
        }
        case Constitutive::eDamageLawType::ISOTROPIC_CUBIC_HERMITE:
        {
            if (rKappa > e_0)
            {
                if (rKappa < e_c)
                {
                    double kappa_scaled = (rKappa - e_0) / (e_c-e_0);
                    omega = 1-e_0/rKappa * ( 2*kappa_scaled*kappa_scaled*kappa_scaled - 3* kappa_scaled*kappa_scaled +1);
                } else {
                    omega = 1.;
                }
            }
            break;
        }
        default:
            throw NuTo::MechanicsException(
                    std::string("[NuTo::GradientDamageEngineeringStress::CalculateDamage] ") +
                    std::string("The required damage law is not implemented. "));

            break;
    }


    return omega;
}

double NuTo::GradientDamageEngineeringStress::CalculateDerivativeDamage(double rKappa)const
{
    double DomegaDkappa = 0;

    double e_0 = mTensileStrength / mE;
    double e_f = mFractureEnergy / mTensileStrength;
    double e_c = 2* e_f; // or something

    switch (mDamageLawType) {
        case Constitutive::eDamageLawType::ISOTROPIC_NO_SOFTENING:
        {
            if (rKappa > e_0)
            {
                DomegaDkappa = e_0 / (rKappa * rKappa);
            }
            break;
        }
        case Constitutive::eDamageLawType::ISOTROPIC_LINEAR_SOFTENING:
        {
            double termA = e_c / MAX_OMEGA / (e_c-e_0);
            double kappa_max = termA*e_0 / (termA-1);

            if (rKappa > kappa_max)
                std::cout << CalculateDamage(kappa_max) << std::endl;

            if (rKappa > e_0 && rKappa < kappa_max)
            {
                DomegaDkappa = e_c * e_0 / (rKappa * rKappa * (e_c - e_0) );
            }
            break;
        }
        case Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING:
        {
            if (rKappa > e_0)
            {
                DomegaDkappa = e_0/rKappa*(1/rKappa + 1/e_f) * exp((e_0 - rKappa)/e_f);
            }
            break;
        }
        case Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING_RES_LOAD:
        {
            if (rKappa > e_0)
            {
                DomegaDkappa = e_0/rKappa*((1/rKappa + 1/e_f) * MAX_OMEGA * exp((e_0 - rKappa)/e_f) + (1-MAX_OMEGA)/rKappa) ;
            }
            break;
        }
        case Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING_SMOOTH:
        {
            double alpha = mDamageLawParameters[0];
            double beta = mDamageLawParameters[1];
            double gamma = mDamageLawParameters[2];

            double e_Dev_e0 = rKappa / e_0;
            double e0_dev_e2 = e_0/pow(rKappa,2.);

            if (e_Dev_e0 <= alpha)
                break;

            if (e_Dev_e0 < beta)
            {
                // pre peak polynomial smoothing
                DomegaDkappa = e0_dev_e2 * (1. + (1.-alpha)/pow(beta-alpha,3.) * pow(e_Dev_e0-beta,2.) * (-2.*e_Dev_e0-beta) );
                break;
            }


            double term = 2. * e_f / (e_0*(gamma-beta));
            double A = -1. / (1+term);
            if (e_Dev_e0 < gamma)
            {
                // post peak polynomial smoothing
                DomegaDkappa = e0_dev_e2 * (1. + (beta*beta-e_Dev_e0*e_Dev_e0)/((gamma-beta)*(gamma-beta-2.*e_f/e_0)));
                break;
            }
            double e_s = e_0*gamma+e_f*log(-A*term);


            // else
            DomegaDkappa = e0_dev_e2 * (rKappa/e_f+1.)*exp((e_s-rKappa)/e_f);

            break;
        }
        case Constitutive::eDamageLawType::ISOTROPIC_CUBIC_HERMITE:
        {
            if (rKappa > e_0 && rKappa < e_c)
            {
                double kappa_scaled = (rKappa - e_0) / (e_c-e_0);
                DomegaDkappa = - 6*e_0/rKappa/(e_c-e_0)*(kappa_scaled*kappa_scaled -kappa_scaled)
                        + e_0/(rKappa*rKappa)*(2*kappa_scaled*kappa_scaled*kappa_scaled-3*kappa_scaled*kappa_scaled+1);

            }
            break;
        }
        default:
            throw NuTo::MechanicsException(
                    std::string("[NuTo::GradientDamageEngineeringStress::CalculateDerivativeDamage] ") +
                    std::string("The required damage law is not implemented. "));

            break;
    }


    return DomegaDkappa;
}

//! @brief ... returns true, if a material model has tmp static data (which has to be updated before stress or stiffness are calculated)
//! @return ... see brief explanation
bool NuTo::GradientDamageEngineeringStress::HaveTmpStaticData() const
{
    return false;
}
