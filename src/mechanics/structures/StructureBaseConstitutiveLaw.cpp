#include "mechanics/structures/StructureBase.h"
#include "mechanics/MechanicsException.h"

#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/constitutive/laws/AdditiveInputExplicit.h"
#include "mechanics/constitutive/laws/AdditiveInputImplicit.h"
#include "mechanics/constitutive/laws/AdditiveOutput.h"
#include "mechanics/constitutive/laws/FibreMatrixBondStressSlip.h"
#include "mechanics/constitutive/laws/GradientDamageEngineeringStress.h"
#include "mechanics/constitutive/laws/GradientDamageFatigueEngineeringStress.h"
#include "mechanics/constitutive/laws/HeatConduction.h"
#include "mechanics/constitutive/laws/LinearDampingEngineeringStress.h"
#include "mechanics/constitutive/laws/LinearElasticEngineeringStress.h"
#include "mechanics/constitutive/laws/LocalDamageModel.h"
#include "mechanics/constitutive/laws/MisesPlasticityEngineeringStress.h"
#include "mechanics/constitutive/laws/MoistureTransport.h"
#include "mechanics/constitutive/laws/PhaseField.h"
#include "mechanics/constitutive/laws/ShrinkageCapillaryStrainBased.h"
#include "mechanics/constitutive/laws/ShrinkageCapillaryStressBased.h"
#include "mechanics/constitutive/laws/ThermalStrains.h"

// create a new constitutive law
int NuTo::StructureBase::ConstitutiveLawCreate(const std::string& rType)
{
    int constitutiveNumber = GetUnusedId(mConstitutiveLawMap);

    this->ConstitutiveLawCreate(constitutiveNumber, Constitutive::ConstitutiveTypeToEnum(rType));
    return constitutiveNumber;
}

// create a new constitutive law
int NuTo::StructureBase::ConstitutiveLawCreate(Constitutive::eConstitutiveType rType)
{
    int constitutiveNumber = GetUnusedId(mConstitutiveLawMap);
    this->ConstitutiveLawCreate(constitutiveNumber, rType);
    return constitutiveNumber;
}

// create a new constitutive law
void NuTo::StructureBase::ConstitutiveLawCreate(int rIdent, Constitutive::eConstitutiveType rType)
{
    using NuTo::Constitutive::eConstitutiveType;
    // check if constitutive law identifier exists
    boost::ptr_map<int, ConstitutiveBase>::iterator it = this->mConstitutiveLawMap.find(rIdent);
    if (it == this->mConstitutiveLawMap.end())
    {
        // create new constitutive law
        ConstitutiveBase* ConstitutiveLawPtr;
        switch (rType)
        {
        case eConstitutiveType::ADDITIVE_INPUT_EXPLICIT:
            ConstitutiveLawPtr = new NuTo::AdditiveInputExplicit(mNumTimeDerivatives);
            break;

        case eConstitutiveType::ADDITIVE_INPUT_IMPLICIT:
            ConstitutiveLawPtr = new NuTo::AdditiveInputImplicit(mNumTimeDerivatives);
            break;

        case eConstitutiveType::ADDITIVE_OUTPUT:
            ConstitutiveLawPtr = new NuTo::AdditiveOutput(mNumTimeDerivatives);
            break;

        case eConstitutiveType::GRADIENT_DAMAGE_ENGINEERING_STRESS:
            ConstitutiveLawPtr = new NuTo::GradientDamageEngineeringStress();
            break;

        case eConstitutiveType::GRADIENT_DAMAGE_FATIGUE_ENGINEERING_STRESS:
            ConstitutiveLawPtr = new NuTo::GradientDamageFatigueEngineeringStress();
            break;

        case eConstitutiveType::HEAT_CONDUCTION:
            ConstitutiveLawPtr = new NuTo::HeatConduction();
            break;

        case eConstitutiveType::LINEAR_DAMPING_ENGINEERING_STRESS:
            ConstitutiveLawPtr = new NuTo::LinearDampingEngineeringStress();
            break;

        case eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS:
            ConstitutiveLawPtr = new NuTo::LinearElasticEngineeringStress();
            break;

        case eConstitutiveType::LOCAL_DAMAGE_MODEL:
            ConstitutiveLawPtr = new NuTo::LocalDamageModel();
            break;

        case eConstitutiveType::MISES_PLASTICITY_ENGINEERING_STRESS:
            ConstitutiveLawPtr = new NuTo::MisesPlasticityEngineeringStress();
            break;

        case eConstitutiveType::MOISTURE_TRANSPORT:
            ConstitutiveLawPtr = new NuTo::MoistureTransport();
            break;

        case eConstitutiveType::FIBRE_MATRIX_BOND_STRESS_SLIP:
            // parameter is the global dimension
            ConstitutiveLawPtr = new NuTo::FibreMatrixBondStressSlip(GetDimension());
            break;

        case eConstitutiveType::SHRINKAGE_CAPILLARY_STRAIN_BASED:
            ConstitutiveLawPtr = new NuTo::ShrinkageCapillaryStrainBased();
            break;

        case eConstitutiveType::SHRINKAGE_CAPILLARY_STRESS_BASED:
            ConstitutiveLawPtr = new NuTo::ShrinkageCapillaryStressBased();
            break;

        case eConstitutiveType::THERMAL_STRAINS:
            ConstitutiveLawPtr = new NuTo::ThermalStrains();
            break;

         default:
            throw NuTo::MechanicsException(__PRETTY_FUNCTION__,
                    "Constitutive law " + Constitutive::ConstitutiveTypeToString(rType) + " currently not supported.");
        }

        // add section to map (insert does not allow const keys!!!!)
        this->mConstitutiveLawMap.insert(rIdent, ConstitutiveLawPtr);
        if (ConstitutiveLawPtr->HaveTmpStaticData())
            mHaveTmpStaticData = true;
    } else
    {
        throw NuTo::MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] Constitutive law already exists.");
    }
}

int NuTo::StructureBase::AddConstitutiveLaw(ConstitutiveBase* rConstitutiveLawPtr)
{
    int constitutiveID = GetUnusedId(mConstitutiveLawMap);
    // add constitutive law to map
    this->mConstitutiveLawMap.insert(constitutiveID, rConstitutiveLawPtr);
    if (rConstitutiveLawPtr->HaveTmpStaticData())
        mHaveTmpStaticData = true;

    return constitutiveID;
}


void NuTo::StructureBase::ConstitutiveLawDelete(int rIdent)
{
    // find constitutive law identifier in map
    boost::ptr_map<int, ConstitutiveBase>::iterator it = this->mConstitutiveLawMap.find(rIdent);
    if (it == this->mConstitutiveLawMap.end())
    {
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Constitutive law does not exist.");
    } else
    {
        this->mConstitutiveLawMap.erase(it);
    }
}

NuTo::ConstitutiveBase* NuTo::StructureBase::ConstitutiveLawGetConstitutiveLawPtr(int rIdent)
{
    boost::ptr_map<int, ConstitutiveBase>::iterator it = this->mConstitutiveLawMap.find(rIdent);
    if (it == this->mConstitutiveLawMap.end())
    {
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Constitutive law does not exist.");
    }
    return it->second;
}

const NuTo::ConstitutiveBase* NuTo::StructureBase::ConstitutiveLawGetConstitutiveLawPtr(int rIdent) const
{
    boost::ptr_map<int, ConstitutiveBase>::const_iterator it = this->mConstitutiveLawMap.find(rIdent);
    if (it == this->mConstitutiveLawMap.end())
    {
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Constitutive law does not exist.");
    }
    return it->second;
}

int NuTo::StructureBase::ConstitutiveLawGetId(const NuTo::ConstitutiveBase* rConstitutiveLawPtr) const
{
    for (boost::ptr_map<int, ConstitutiveBase>::const_iterator it = mConstitutiveLawMap.begin(); it != mConstitutiveLawMap.end(); it++)
    {
        if (it->second == rConstitutiveLawPtr)
        {
            return it->first;
        }
    }
    throw MechanicsException(__PRETTY_FUNCTION__, "Constitutive law does not exist.");
}

void NuTo::StructureBase::ConstitutiveLawInfo(unsigned short rVerboseLevel) const
{
    std::cout << "Number of constitutive laws: " << this->GetNumConstitutiveLaws() << std::endl;
    for (boost::ptr_map<int, ConstitutiveBase>::const_iterator it = mConstitutiveLawMap.begin(); it != mConstitutiveLawMap.end(); it++)
    {
        std::cout << "  Constitutive law: " << it->first << std::endl;
        it->second->Info(rVerboseLevel, mLogger);
    }
}

void NuTo::StructureBase::ConstitutiveLawInfo(int rIdent, unsigned short rVerboseLevel) const
{
    const ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
    std::cout << "  Constitutive law: " << rIdent << std::endl;
    ConstitutiveLawPtr->Info(rVerboseLevel, mLogger);
}

bool NuTo::StructureBase::ConstitutiveLawGetParameterBool(int rIdent, const std::string &rIdentifier) const
{
    return ConstitutiveLawGetParameterBool(rIdent, Constitutive::ConstitutiveParameterToEnum(rIdentifier));
}

void NuTo::StructureBase::ConstitutiveLawSetParameterBool(int rIdent, const std::string &rIdentifier, bool rValue)
{
   ConstitutiveLawSetParameterBool(rIdent, Constitutive::ConstitutiveParameterToEnum(rIdentifier), rValue);
}

double NuTo::StructureBase::ConstitutiveLawGetParameterDouble(int rIdent, const std::string &rIdentifier) const
{
    return ConstitutiveLawGetParameterDouble(rIdent, Constitutive::ConstitutiveParameterToEnum(rIdentifier));
}

void NuTo::StructureBase::ConstitutiveLawSetParameterDouble(int rIdent, const std::string &rIdentifier, double rValue)
{
    ConstitutiveLawSetParameterDouble(rIdent, Constitutive::ConstitutiveParameterToEnum(rIdentifier), rValue);
}

Eigen::VectorXd NuTo::StructureBase::ConstitutiveLawGetParameterFullVectorDouble(int rIdent, const std::string &rIdentifier) const
{
    return ConstitutiveLawGetParameterFullVectorDouble(rIdent, Constitutive::ConstitutiveParameterToEnum(rIdentifier));
}

void NuTo::StructureBase::ConstitutiveLawSetParameterFullVectorDouble(int rIdent, const std::string &rIdentifier, Eigen::VectorXd rValue)
{
    ConstitutiveLawSetParameterFullVectorDouble(rIdent, Constitutive::ConstitutiveParameterToEnum(rIdentifier), rValue);
}

bool NuTo::StructureBase::ConstitutiveLawGetParameterBool(int rIdent, NuTo::Constitutive::eConstitutiveParameter rIdentifier) const
{
    try
    {
        const ConstitutiveBase* constitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        return constitutiveLawPtr->GetParameterBool(rIdentifier);
    } catch (NuTo::MechanicsException& e)
    {
        e.AddMessage(__PRETTY_FUNCTION__, "error getting requested value.");
        throw;
    }
}

void NuTo::StructureBase::ConstitutiveLawSetParameterBool(int rIdent, NuTo::Constitutive::eConstitutiveParameter rIdentifier, bool rValue)
{
    try
    {
        ConstitutiveBase* constitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        constitutiveLawPtr->SetParameterBool(rIdentifier, rValue);
    } catch (NuTo::MechanicsException& e)
    {
        e.AddMessage(__PRETTY_FUNCTION__, "error setting requested value.");
        throw;
    }
}


double NuTo::StructureBase::ConstitutiveLawGetParameterDouble(int rIdent, NuTo::Constitutive::eConstitutiveParameter rIdentifier) const
{
    try
    {
        const ConstitutiveBase* constitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        return constitutiveLawPtr->GetParameterDouble(rIdentifier);
    } catch (NuTo::MechanicsException& e)
    {
        e.AddMessage(__PRETTY_FUNCTION__, "error getting requested value.");
        throw;
    }
}

void NuTo::StructureBase::ConstitutiveLawSetParameterDouble(int rIdent, NuTo::Constitutive::eConstitutiveParameter rIdentifier, double rValue)
{
    try
    {
        ConstitutiveBase* constitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        constitutiveLawPtr->SetParameterDouble(rIdentifier, rValue);
    } catch (NuTo::MechanicsException& e)
    {
        e.AddMessage(__PRETTY_FUNCTION__, "error setting requested value.");
        throw;
    }
}

Eigen::VectorXd NuTo::StructureBase::ConstitutiveLawGetParameterFullVectorDouble(int rIdent, NuTo::Constitutive::eConstitutiveParameter rIdentifier) const
{
    try
    {
        const ConstitutiveBase* constitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        return constitutiveLawPtr->GetParameterFullVectorDouble(rIdentifier);
    } catch (NuTo::MechanicsException& e)
    {
        e.AddMessage(__PRETTY_FUNCTION__, "error getting requested value.");
        throw;
    }
}
void NuTo::StructureBase::ConstitutiveLawSetParameterFullVectorDouble(int rIdent, NuTo::Constitutive::eConstitutiveParameter rIdentifier, Eigen::VectorXd rValue)
{
    try
    {
        ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        ConstitutiveLawPtr->SetParameterFullVectorDouble(rIdentifier, rValue);
    } catch (NuTo::MechanicsException& e)
    {
        e.AddMessage(__PRETTY_FUNCTION__, "error setting requested value.");
        throw;
    }
}


void NuTo::StructureBase::ConstitutiveLawSetDamageLaw(int lawId, std::shared_ptr<Constitutive::DamageLaw> damageLaw)
{
    try
    {
        ConstitutiveBase* constitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(lawId);
        constitutiveLawPtr->SetDamageLaw(damageLaw);
    } catch (NuTo::MechanicsException& e)
    {
        e.AddMessage(__PRETTY_FUNCTION__, "error setting damage law.");
        throw;
    }
}

double NuTo::StructureBase::ConstitutiveLawGetEquilibriumWaterVolumeFraction(int rIdent, double rRelativeHumidity, Eigen::VectorXd rCoeffs) const
{
    double EquilibriumWaterVolumeFraction = 0.0;
    try
    {
        const ConstitutiveBase* constitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        EquilibriumWaterVolumeFraction = constitutiveLawPtr->GetEquilibriumWaterVolumeFraction(rRelativeHumidity, rCoeffs);
    } catch (NuTo::MechanicsException& e)
    {
        e.AddMessage(__PRETTY_FUNCTION__, "error getting the equilibrium water volume fraction.");
        throw;
    }
    return EquilibriumWaterVolumeFraction;
}
