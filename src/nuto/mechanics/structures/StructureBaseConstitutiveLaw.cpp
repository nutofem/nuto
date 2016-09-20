#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/MechanicsException.h"

#include "nuto/mechanics/constitutive/laws/AdditiveInputExplicit.h"
#include "nuto/mechanics/constitutive/laws/AdditiveInputImplicit.h"
#include "nuto/mechanics/constitutive/laws/AdditiveOutput.h"
#include "nuto/mechanics/constitutive/laws/GradientDamageEngineeringStress.h"
#include "nuto/mechanics/constitutive/laws/HeatConduction.h"
#include "nuto/mechanics/constitutive/laws/FibreMatrixBondStressSlip.h"
#include "nuto/mechanics/constitutive/laws/LinearElasticEngineeringStress.h"
#include "nuto/mechanics/constitutive/laws/MoistureTransport.h"
#include "nuto/mechanics/constitutive/laws/MisesPlasticityEngineeringStress.h"
#include "nuto/mechanics/constitutive/laws/PhaseField.h"
#include "nuto/mechanics/constitutive/laws/ShrinkageCapillaryStrainBased.h"
#include "nuto/mechanics/constitutive/laws/ShrinkageCapillaryStressBased.h"
#include "nuto/mechanics/constitutive/laws/ThermalStrains.h"

// create a new constitutive law
int NuTo::StructureBase::ConstitutiveLawCreate(const std::string& rType)
{
    //find unused integer id
    int constitutiveNumber = mConstitutiveLawMap.size() == 0 ? 0 : mConstitutiveLawMap.rbegin()->first +1;

    this->ConstitutiveLawCreate(constitutiveNumber, Constitutive::ConstitutiveTypeToEnum(rType));
    return constitutiveNumber;
}

// create a new constitutive law
int NuTo::StructureBase::ConstitutiveLawCreate(Constitutive::eConstitutiveType rType)
{
    //find unused integer id
    int constitutiveNumber = mConstitutiveLawMap.size();
    boost::ptr_map<int, ConstitutiveBase>::iterator it = mConstitutiveLawMap.find(constitutiveNumber);
    while (it != mConstitutiveLawMap.end())
    {
        constitutiveNumber++;
        it = mConstitutiveLawMap.find(constitutiveNumber);
    }

    this->ConstitutiveLawCreate(constitutiveNumber, rType);
    return constitutiveNumber;
}

// create a new constitutive law
void NuTo::StructureBase::ConstitutiveLawCreate(int rIdent, Constitutive::eConstitutiveType rType)
{
    // check if constitutive law identifier exists
    boost::ptr_map<int, ConstitutiveBase>::iterator it = this->mConstitutiveLawMap.find(rIdent);
    if (it == this->mConstitutiveLawMap.end())
    {
        // create new constitutive law
        ConstitutiveBase* ConstitutiveLawPtr;
        switch (rType)
        {
        case NuTo::Constitutive::ADDITIVE_INPUT_EXPLICIT:
            ConstitutiveLawPtr = new NuTo::AdditiveInputExplicit();
            break;

        case NuTo::Constitutive::ADDITIVE_INPUT_IMPLICIT:
            ConstitutiveLawPtr = new NuTo::AdditiveInputImplicit();
            break;
        case NuTo::Constitutive::CONSTITUTIVE_LAWS_ADDITIVE_OUTPUT:
            ConstitutiveLawPtr = new NuTo::ConstitutiveLawsAdditiveOutput();
            break;

        case NuTo::Constitutive::DAMAGE_VISCO_PLASTICITY_ENGINEERING_STRESS:
            throw NuTo::MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] constitutive law " + Constitutive::ConstitutiveTypeToString(rType) + " currently not supported.");
//            ConstitutiveLawPtr = new NuTo::DamageViscoPlasticityEngineeringStress();
            break;

        case NuTo::Constitutive::DAMAGE_VISCO_PLASTICITY_HARDENING_ENGINEERING_STRESS:
            throw NuTo::MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] constitutive law " + Constitutive::ConstitutiveTypeToString(rType) + " currently not supported.");
//            ConstitutiveLawPtr = new NuTo::DamageViscoPlasticityHardeningEngineeringStress();
            break;

        case NuTo::Constitutive::GRADIENT_DAMAGE_ENGINEERING_STRESS:
            ConstitutiveLawPtr = new NuTo::GradientDamageEngineeringStress();
            break;

        case NuTo::Constitutive::GRADIENT_DAMAGE_ENGINEERING_STRESS_FATIGUE:
            throw NuTo::MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] constitutive law " + Constitutive::ConstitutiveTypeToString(rType) + " currently not supported.");
//            ConstitutiveLawPtr = new NuTo::GradientDamageEngineeringStressFatigue();
            break;

        case NuTo::Constitutive::GRADIENT_DAMAGE_PLASTICITY_ENGINEERING_STRESS:
            throw NuTo::MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] constitutive law " + Constitutive::ConstitutiveTypeToString(rType) + " currently not supported.");
//            ConstitutiveLawPtr = new NuTo::GradientDamagePlasticityEngineeringStress();
            break;

        case NuTo::Constitutive::HEAT_CONDUCTION:
            ConstitutiveLawPtr = new NuTo::HeatConduction();
            break;

        case NuTo::Constitutive::LINEAR_ELASTIC_ENGINEERING_STRESS:
            ConstitutiveLawPtr = new NuTo::LinearElasticEngineeringStress();
            break;

        case NuTo::Constitutive::LINEAR_SPRING:
            throw NuTo::MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] constitutive law " + Constitutive::ConstitutiveTypeToString(rType) + " currently not supported.");
//            ConstitutiveLawPtr = new NuTo::LinearSpring();
            break;

        case NuTo::Constitutive::MISES_PLASTICITY_ENGINEERING_STRESS:
            ConstitutiveLawPtr = new NuTo::MisesPlasticityEngineeringStress();
            break;

        case NuTo::Constitutive::MOISTURE_TRANSPORT:
            ConstitutiveLawPtr = new NuTo::MoistureTransport();
            break;

        case NuTo::Constitutive::NONLOCAL_DAMAGE_PLASTICITY_ENGINEERING_STRESS:
            throw NuTo::MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] constitutive law " + Constitutive::ConstitutiveTypeToString(rType) + " currently not supported.");
//            ConstitutiveLawPtr = new NuTo::NonlocalDamagePlasticityEngineeringStress();
            break;

        case NuTo::Constitutive::FIBRE_MATRIX_BOND_STRESS_SLIP:
            ConstitutiveLawPtr = new NuTo::FibreMatrixBondStressSlip(GetDimension());
            break;

        case NuTo::Constitutive::PHASE_FIELD:
            throw NuTo::MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] constitutive law " + Constitutive::ConstitutiveTypeToString(rType) + " currently not supported.");
            //ConstitutiveLawPtr = new NuTo::PhaseField();
            break;

        case NuTo::Constitutive::SHRINKAGE_CAPILLARY_STRAIN_BASED:
            ConstitutiveLawPtr = new NuTo::ShrinkageCapillaryStrainBased();
            break;

        case NuTo::Constitutive::SHRINKAGE_CAPILLARY_STRESS_BASED:
            ConstitutiveLawPtr = new NuTo::ShrinkageCapillaryStressBased();
            break;

        case NuTo::Constitutive::STRAIN_GRADIENT_DAMAGE_PLASTICITY_ENGINEERING_STRESS:
            throw NuTo::MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] constitutive law " + Constitutive::ConstitutiveTypeToString(rType) + " currently not supported.");
//            ConstitutiveLawPtr = new NuTo::StrainGradientDamagePlasticityEngineeringStress();
            break;

        case NuTo::Constitutive::THERMAL_STRAINS:
            ConstitutiveLawPtr = new NuTo::ThermalStrains();
            break;



         default:
            throw NuTo::MechanicsException("[NuTo::StructureBase::ConstitutiveLawCreate] invalid type of constitutive law.");
        }

        // add section to map (insert does not allow const keys!!!!)
        this->mConstitutiveLawMap.insert(rIdent, ConstitutiveLawPtr);
        if (ConstitutiveLawPtr->HaveTmpStaticData())
            mHaveTmpStaticData = true;
    } else
    {
        throw NuTo::MechanicsException("[NuTo::StructureBase::ConstitutiveLawCreate] Constitutive law already exists.");
    }
}

int NuTo::StructureBase::AddConstitutiveLaw(ConstitutiveBase* rConstitutiveLawPtr)
{
    //find unused integer id
    int constitutiveID = mConstitutiveLawMap.size();
    boost::ptr_map<int, ConstitutiveBase>::iterator it = mConstitutiveLawMap.find(constitutiveID);
    while (it != mConstitutiveLawMap.end())
    {
        constitutiveID++;
        it = mConstitutiveLawMap.find(constitutiveID);
    }

    // add constitutive law to map
    this->mConstitutiveLawMap.insert(constitutiveID, rConstitutiveLawPtr);
    if (rConstitutiveLawPtr->HaveTmpStaticData())
        mHaveTmpStaticData = true;

    return constitutiveID;
}


// delete an existing constitutive law
void NuTo::StructureBase::ConstitutiveLawDelete(int rIdent)
{
    // find constitutive law identifier in map
    boost::ptr_map<int, ConstitutiveBase>::iterator it = this->mConstitutiveLawMap.find(rIdent);
    if (it == this->mConstitutiveLawMap.end())
    {
        throw NuTo::MechanicsException("[NuTo::StructureBase::ConstitutiveLawDelete] Constitutive law does not exist.");
    } else
    {
        this->mConstitutiveLawMap.erase(it);
    }
}

// get constitutive law pointer from constitutive law identifier
NuTo::ConstitutiveBase* NuTo::StructureBase::ConstitutiveLawGetConstitutiveLawPtr(int rIdent)
{
    boost::ptr_map<int, ConstitutiveBase>::iterator it = this->mConstitutiveLawMap.find(rIdent);
    if (it == this->mConstitutiveLawMap.end())
    {
        throw NuTo::MechanicsException("[NuTo::StructureBase::ConstitutiveLawGetConstitutiveLawPtr] Constitutive law does not exist.");
    }
    return it->second;
}

// get constitutive law pointer from constitutive law identifier
const NuTo::ConstitutiveBase* NuTo::StructureBase::ConstitutiveLawGetConstitutiveLawPtr(int rIdent) const
{
    boost::ptr_map<int, ConstitutiveBase>::const_iterator it = this->mConstitutiveLawMap.find(rIdent);
    if (it == this->mConstitutiveLawMap.end())
    {
        throw NuTo::MechanicsException("[NuTo::StructureBase::ConstitutiveLawGetConstitutiveLawPtr] Constitutive law does not exist.");
    }
    return it->second;
}

// get constitutive law identifier from constitutive law pointer
int NuTo::StructureBase::ConstitutiveLawGetId(const NuTo::ConstitutiveBase* rConstitutiveLawPtr) const
{
    for (boost::ptr_map<int, ConstitutiveBase>::const_iterator it = mConstitutiveLawMap.begin(); it != mConstitutiveLawMap.end(); it++)
    {
        if (it->second == rConstitutiveLawPtr)
        {
            return it->first;
        }
    }
    throw MechanicsException("[NuTo::StructureBase::ConstitutiveLawGetId] Constitutive law does not exist.");
}

// info routines
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

//! @brief ... gets a variable of the constitutive law which is selected by a string
//! @param rIdent ... constitutive law identifier
//! @param rIdentifier ... String to identify the requested variable
//! @return ... value of the requested variable
bool NuTo::StructureBase::ConstitutiveLawGetParameterBool(int rIdent, const std::string &rIdentifier) const
{
    return ConstitutiveLawGetParameterBool(rIdent, Constitutive::ConstitutiveParameterToEnum(rIdentifier));
}

//! @brief ... sets a variable of the constitutive law which is selected by a string
//! @param rIdent ... constitutive law identifier
//! @param rIdentifier ... String to identify the requested variable
//! @param rValue ... new value for requested variable
void NuTo::StructureBase::ConstitutiveLawSetParameterBool(int rIdent, const std::string &rIdentifier, bool rValue)
{
   ConstitutiveLawSetParameterBool(rIdent, Constitutive::ConstitutiveParameterToEnum(rIdentifier), rValue);
}

//! @brief ... gets a variable of the constitutive law which is selected by a string
//! @param rIdent ... constitutive law identifier
//! @param rIdentifier ... String to identify the requested variable
//! @return ... value of the requested variable
double NuTo::StructureBase::ConstitutiveLawGetParameterDouble(int rIdent, const std::string &rIdentifier) const
{
    return ConstitutiveLawGetParameterDouble(rIdent, Constitutive::ConstitutiveParameterToEnum(rIdentifier));
}

//! @brief ... sets a variable of the constitutive law which is selected by a sting
//! @param rIdent ... constitutive law identifier
//! @param rIdentifier ... String to identify the requested variable
//! @param rValue ... new value for requested variable
void NuTo::StructureBase::ConstitutiveLawSetParameterDouble(int rIdent, const std::string &rIdentifier, double rValue)
{
    ConstitutiveLawSetParameterDouble(rIdent, Constitutive::ConstitutiveParameterToEnum(rIdentifier), rValue);
}

//! @brief ... gets a variable of the constitutive law which is selected by a string
//! @param rIdent ... constitutive law identifier
//! @param rIdentifier ... String to identify the requested variable
//! @return ... value of the requested variable
NuTo::FullVector<double, Eigen::Dynamic> NuTo::StructureBase::ConstitutiveLawGetParameterFullVectorDouble(int rIdent, const std::string &rIdentifier) const
{
    return ConstitutiveLawGetParameterFullVectorDouble(rIdent, Constitutive::ConstitutiveParameterToEnum(rIdentifier));
}

//! @brief ... sets a variable of the constitutive law which is selected by a string
//! @param rIdent ... constitutive law identifier
//! @param rIdentifier ... String to identify the requested variable
//! @param rValue ... new value for requested variable
void NuTo::StructureBase::ConstitutiveLawSetParameterFullVectorDouble(int rIdent, const std::string &rIdentifier, NuTo::FullVector<double, Eigen::Dynamic> rValue)
{
    ConstitutiveLawSetParameterFullVectorDouble(rIdent, Constitutive::ConstitutiveParameterToEnum(rIdentifier), rValue);
}

//! @brief ... gets a parameter of the constitutive law which is selected by an enum
//! @param rIdent ... constitutive law identifier
//! @param rIdentifier ... Enum to identify the requested parameter
//! @return ... value of the requested variable
bool NuTo::StructureBase::ConstitutiveLawGetParameterBool(int rIdent, NuTo::Constitutive::eConstitutiveParameter rIdentifier) const
{
    bool requestedValue;
    try
    {
        const ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        requestedValue = ConstitutiveLawPtr->GetParameterBool(rIdentifier);
    } catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawGetParameterBool] error getting requested value.");
        throw e;
    }
    return requestedValue;
}

//! @brief ... sets a parameter of the constitutive law which is selected by an enum
//! @param rIdent ... constitutive law identifier
//! @param rIdentifier ... Enum to identify the requested parameter
//! @param rValue ... new value for requested variable
void NuTo::StructureBase::ConstitutiveLawSetParameterBool(int rIdent, NuTo::Constitutive::eConstitutiveParameter rIdentifier, bool rValue)
{
    try
    {
        ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        ConstitutiveLawPtr->SetParameterBool(rIdentifier, rValue);
    } catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawSetParameterBool] error setting requested value.");
        throw e;
    }
}



//! @brief ... gets a parameter of the constitutive law which is selected by an enum
//! @param rIdent ... constitutive law identifier
//! @param rIdentifier ... Enum to identify the requested parameter
//! @return ... value of the requested variable
double NuTo::StructureBase::ConstitutiveLawGetParameterDouble(int rIdent, NuTo::Constitutive::eConstitutiveParameter rIdentifier) const
{
    double requestedValue;
    try
    {
        const ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        requestedValue = ConstitutiveLawPtr->GetParameterDouble(rIdentifier);
    } catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawGetParameterDouble] error getting requested value.");
        throw e;
    }
    return requestedValue;
}

//! @brief ... sets a parameter of the constitutive law which is selected by an enum
//! @param rIdent ... constitutive law identifier
//! @param rIdentifier ... Enum to identify the requested parameter
//! @param rValue ... new value for requested variable
void NuTo::StructureBase::ConstitutiveLawSetParameterDouble(int rIdent, NuTo::Constitutive::eConstitutiveParameter rIdentifier, double rValue)
{
    try
    {
        ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        ConstitutiveLawPtr->SetParameterDouble(rIdentifier, rValue);
    } catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawSetParameterDouble] error setting requested value.");
        throw e;
    }
}

//! @brief ... gets a parameter of the constitutive law which is selected by an enum
//! @param rIdent ... constitutive law identifier
//! @param rIdentifier ... Enum to identify the requested parameter
//! @return ... value of the requested variable
NuTo::FullVector<double, Eigen::Dynamic> NuTo::StructureBase::ConstitutiveLawGetParameterFullVectorDouble(int rIdent, NuTo::Constitutive::eConstitutiveParameter rIdentifier) const
{
    NuTo::FullVector<double, Eigen::Dynamic> requestedValue;
    try
    {
        const ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        requestedValue = ConstitutiveLawPtr->GetParameterFullVectorDouble(rIdentifier);
    } catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawGetParameterFullVectorDouble] error getting requested value.");
        throw e;
    }
    return requestedValue;
}

//! @brief ... sets a parameter of the constitutive law which is selected by an enum
//! @param rIdent ... constitutive law identifier
//! @param rIdentifier ... Enum to identify the requested parameter
//! @param rValue ... new value for requested variable
void NuTo::StructureBase::ConstitutiveLawSetParameterFullVectorDouble(int rIdent, NuTo::Constitutive::eConstitutiveParameter rIdentifier, NuTo::FullVector<double, Eigen::Dynamic> rValue)
{
    try
    {
        ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        ConstitutiveLawPtr->SetParameterFullVectorDouble(rIdentifier, rValue);
    } catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawSetParameterFullVectorDouble] error setting requested value.");
        throw e;
    }
}

//! @brief ... get yield strength for multilinear response
//! @return ... first column: equivalent plastic strain
//! @return ... second column: corresponding yield strength
NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> NuTo::StructureBase::ConstitutiveLawGetYieldStrength(int rIdent) const
{
    try
    {
        const ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        return ConstitutiveLawPtr->GetYieldStrength();
    } catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawGetYieldStrength] error getting yield strength.");
        throw e;
    }
}

//! @brief ... add yield strength
//! @param rEpsilon ...  equivalent plastic strain
//! @param rSigma ...  yield strength
void NuTo::StructureBase::ConstitutiveLawAddYieldStrength(int rIdent, double rEpsilon, double rSigma)
{
    try
    {
        ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        ConstitutiveLawPtr->AddYieldStrength(rEpsilon, rSigma);
    } catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawAddYieldStrength] error adding yield strength.");
        throw e;
    }
}

//! @brief ... get hardening modulus for multilinear response
//! @return ... first column: equivalent plastic strain
//! @return ... second column: corresponding hardening modulus
NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> NuTo::StructureBase::ConstitutiveLawGetHardeningModulus(int rIdent) const
{
    try
    {
        const ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        return ConstitutiveLawPtr->GetHardeningModulus();
    } catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawGetHardeningModulus's ratio] error getting hardening moduli.");
        throw e;
    }
}

//! @brief ... add hardening modulus
//! @param rEpsilon ...  equivalent plastic strain
//! @param rSigma ...  hardening modulus
void NuTo::StructureBase::ConstitutiveLawAddHardeningModulus(int rIdent, double rEpsilon, double rH)
{
    try
    {
        ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        ConstitutiveLawPtr->AddHardeningModulus(rEpsilon, rH);
    } catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawAddHardeningModulus] error adding hardening modulus.");
        throw e;
    }
}


void NuTo::StructureBase::ConstitutiveLawSetDamageLaw(int rIdent, Constitutive::eDamageLawType rDamageLaw)
{
    try
    {
        ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        ConstitutiveLawPtr->SetParameterDouble(Constitutive::eConstitutiveParameter::DAMAGE_LAW, rDamageLaw);
    } catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawSetDamageLaw] error setting damage law.");
        throw e;
    }
}

//****************************

//! @brief ... gets the equilibrium water volume fraction depend on the relative humidity
//! @param rIdent ... constitutive law identifier
//! @param rRelativeHumidity ... relative humidity
//! @param rCoeffs ... polynomial coefficients of the sorption curve
//! @return ... equilibrium water volume fraction
double NuTo::StructureBase::ConstitutiveLawGetEquilibriumWaterVolumeFraction(int rIdent, double rRelativeHumidity, NuTo::FullVector<double, Eigen::Dynamic> rCoeffs) const
{
    double EquilibriumWaterVolumeFraction = 0.0;
    try
    {
        const ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        EquilibriumWaterVolumeFraction = ConstitutiveLawPtr->GetEquilibriumWaterVolumeFraction(rRelativeHumidity, rCoeffs);
    } catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawGetEquilibriumWaterVolumeFraction] error getting the equilibrium water volume fraction.");
        throw e;
    }
    return EquilibriumWaterVolumeFraction;
}

//VHIRTHAMTODO Delete???
////! @brief ... adds a constitutive law to a multi physics model
////! @param rIdentMultiPhysics ... multi physics constitutive law to which the constitutive law should be added
////! @param rIdentConstitutiveLaw ... constitutive law which should be added to the multi physics model
//void NuTo::StructureBase::ConstitutiveLawMultiPhysicsAddConstitutiveLaw(int rIdentMultiPhysics, int rIdentConstitutiveLaw)
//{
//    try
//    {
//        ConstitutiveBase* MultiPhysicsLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdentMultiPhysics);
//        ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdentConstitutiveLaw);

//        MultiPhysicsLawPtr->MultiPhysicsAddConstitutiveLaw(ConstitutiveLawPtr);
//    }
//    catch (NuTo::MechanicsException& e)
//    {
//        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawMultiPhysicsAddConstitutiveLaw] error adding constitutive law to multi physics model.");
//        throw e;
//    }
//}
