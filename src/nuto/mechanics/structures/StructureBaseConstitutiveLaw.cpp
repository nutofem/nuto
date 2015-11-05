// $Id$


#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/constitutive/mechanics/LinearSpring.h"
#include "nuto/mechanics/constitutive/mechanics/LinearElasticEngineeringStress.h"
#include "nuto/mechanics/constitutive/mechanics/DamageViscoPlasticityEngineeringStress.h"
#include "nuto/mechanics/constitutive/mechanics/DamageViscoPlasticityHardeningEngineeringStress.h"
#include "nuto/mechanics/constitutive/mechanics/MisesPlasticityEngineeringStress.h"
#include "nuto/mechanics/constitutive/mechanics/GradientDamagePlasticityEngineeringStress.h"
#include "nuto/mechanics/constitutive/mechanics/FibreMatrixBondStressSlip.h"
#include "nuto/mechanics/constitutive/mechanics/GradientDamageEngineeringStress.h"
#include "nuto/mechanics/constitutive/mechanics/NonlocalDamagePlasticityEngineeringStress.h"
#include "nuto/mechanics/constitutive/mechanics/StrainGradientDamagePlasticityEngineeringStress.h"
#include "nuto/mechanics/constitutive/moistureTransport/MoistureTransport.h"
#include "nuto/mechanics/constitutive/multiPhysics/ConstitutiveMultiPhysics.h"
#include <nuto/mechanics/constitutive/shrinkage/DryingShrinkage.h>
#include "nuto/mechanics/constitutive/thermal/LinearHeatFlux.h"

// create a new constitutive law
int NuTo::StructureBase::ConstitutiveLawCreate(const std::string& rType)
{
    // convert section type string to upper case
    std::string ConstitutiveLawTypeString;
    std::transform(rType.begin(), rType.end(), std::back_inserter(ConstitutiveLawTypeString), (int (*)(int)) toupper);

    // get section type from string
Constitutive    ::eConstitutiveType ConstitutiveLawType;
    if (ConstitutiveLawTypeString == "LINEARELASTICENGINEERINGSTRESS")
    {
        ConstitutiveLawType = Constitutive::LINEAR_ELASTIC_ENGINEERING_STRESS;
    } else if (ConstitutiveLawTypeString == "MISESPLASTICITYENGINEERINGSTRESS")
    {
        ConstitutiveLawType = Constitutive::MISES_PLASTICITY_ENGINEERING_STRESS;
    } else if (ConstitutiveLawTypeString == "NONLOCALDAMAGEPLASTICITYENGINEERINGSTRESS")
    {
        ConstitutiveLawType = Constitutive::NONLOCAL_DAMAGE_PLASTICITY_ENGINEERING_STRESS;
    } else if (ConstitutiveLawTypeString == "GRADIENTDAMAGEPLASTICITYENGINEERINGSTRESS")
    {
        ConstitutiveLawType = Constitutive::GRADIENT_DAMAGE_PLASTICITY_ENGINEERING_STRESS;
    } else if (ConstitutiveLawTypeString == "GRADIENTDAMAGEENGINEERINGSTRESS")
    {
        ConstitutiveLawType = Constitutive::GRADIENT_DAMAGE_ENGINEERING_STRESS;
    } else if (ConstitutiveLawTypeString == "LINEARHEATFLUX")
    {
        ConstitutiveLawType = Constitutive::LINEAR_HEAT_FLUX;
    } else if (ConstitutiveLawTypeString == "STRAINGRADIENTDAMAGEPLASTICITYENGINEERINGSTRESS")
    {
        ConstitutiveLawType = Constitutive::STRAIN_GRADIENT_DAMAGE_PLASTICITY_ENGINEERING_STRESS;
    } else if (ConstitutiveLawTypeString == "DAMAGEVISCOPLASTICITYENGINEERINGSTRESS")
    {
        ConstitutiveLawType = Constitutive::DAMAGE_VISCO_PLASTICITY_ENGINEERING_STRESS;
    } else if (ConstitutiveLawTypeString == "DAMAGEVISCOPLASTICITYHARDENINGENGINEERINGSTRESS")
    {
        ConstitutiveLawType = Constitutive::DAMAGE_VISCO_PLASTICITY_HARDENING_ENGINEERING_STRESS;
    } else if (ConstitutiveLawTypeString == "MOISTURETRANSPORT")
    {
        ConstitutiveLawType = Constitutive::MOISTURE_TRANSPORT;
    } else if (ConstitutiveLawTypeString == "MULTIPHYSICS")
    {
        ConstitutiveLawType = Constitutive::MULTI_PHYSICS;
    } else if (ConstitutiveLawTypeString == "DRYINGSHRINKAGE")
    {
        ConstitutiveLawType = Constitutive::DRYING_SHRINKAGE;
    } else
    {
        throw NuTo::MechanicsException("[NuTo::StructureBase::ConstitutiveLawCreate] invalid type of constitutive law.");
    }

    //find unused integer id
    int constitutiveNumber(mConstitutiveLawMap.size());
    boost::ptr_map<int, ConstitutiveBase>::iterator it = mConstitutiveLawMap.find(constitutiveNumber);
    while (it != mConstitutiveLawMap.end())
    {
        constitutiveNumber++;
        it = mConstitutiveLawMap.find(constitutiveNumber);
    }

    this->ConstitutiveLawCreate(constitutiveNumber, ConstitutiveLawType);
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
        case NuTo::Constitutive::LINEAR_ELASTIC_ENGINEERING_STRESS:
            ConstitutiveLawPtr = new NuTo::LinearElasticEngineeringStress();
            break;
        case NuTo::Constitutive::LINEAR_SPRING:
            ConstitutiveLawPtr = new NuTo::LinearSpring();
            break;
        case NuTo::Constitutive::MISES_PLASTICITY_ENGINEERING_STRESS:
            ConstitutiveLawPtr = new NuTo::MisesPlasticityEngineeringStress();
            break;
        case NuTo::Constitutive::NONLOCAL_DAMAGE_PLASTICITY_ENGINEERING_STRESS:
            ConstitutiveLawPtr = new NuTo::NonlocalDamagePlasticityEngineeringStress();
            break;
        case NuTo::Constitutive::LINEAR_HEAT_FLUX:
            ConstitutiveLawPtr = new NuTo::LinearHeatFlux();
            break;
        case NuTo::Constitutive::GRADIENT_DAMAGE_PLASTICITY_ENGINEERING_STRESS:
            ConstitutiveLawPtr = new NuTo::GradientDamagePlasticityEngineeringStress();
            break;
        case NuTo::Constitutive::GRADIENT_DAMAGE_ENGINEERING_STRESS:
            ConstitutiveLawPtr = new NuTo::GradientDamageEngineeringStress();
            break;
        case NuTo::Constitutive::STRAIN_GRADIENT_DAMAGE_PLASTICITY_ENGINEERING_STRESS:
            ConstitutiveLawPtr = new NuTo::StrainGradientDamagePlasticityEngineeringStress();
            break;
        case NuTo::Constitutive::DAMAGE_VISCO_PLASTICITY_ENGINEERING_STRESS:
            ConstitutiveLawPtr = new NuTo::DamageViscoPlasticityEngineeringStress();
            break;
        case NuTo::Constitutive::DAMAGE_VISCO_PLASTICITY_HARDENING_ENGINEERING_STRESS:
            ConstitutiveLawPtr = new NuTo::DamageViscoPlasticityHardeningEngineeringStress();
            break;
        case NuTo::Constitutive::MOISTURE_TRANSPORT:
            ConstitutiveLawPtr = new NuTo::MoistureTransport();
            break;
        case NuTo::Constitutive::MULTI_PHYSICS:
            ConstitutiveLawPtr = new NuTo::ConstitutiveMultiPhysics();
            break;
        case NuTo::Constitutive::FIBRE_MATRIX_BOND_STRESS_SLIP:
            ConstitutiveLawPtr = new NuTo::FibreMatrixBondStressSlip();
            break;
        case NuTo::Constitutive::DRYING_SHRINKAGE:
            ConstitutiveLawPtr = new NuTo::DryingShrinkage();
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
    std::string upperCaseIdentifier;
    NuTo::Constitutive::eConstitutiveParameter IdentifierEnum;

    // transform string to uppercase
    std::transform(rIdentifier.begin(), rIdentifier.end(), std::back_inserter(upperCaseIdentifier), (int (*)(int)) toupper);

    // Get enum
IdentifierEnum    = Constitutive::GetConstitutiveVariableFromString(upperCaseIdentifier);

    return ConstitutiveLawGetParameterBool(rIdent, IdentifierEnum);
}

//! @brief ... sets a variable of the constitutive law which is selected by a string
//! @param rIdent ... constitutive law identifier
//! @param rIdentifier ... String to identify the requested variable
//! @param rValue ... new value for requested variable
void NuTo::StructureBase::ConstitutiveLawSetParameterBool(int rIdent, const std::string &rIdentifier, bool rValue)
{
    std::string upperCaseIdentifier;
    NuTo::Constitutive::eConstitutiveParameter IdentifierEnum;

    // transform string to uppercase
    std::transform(rIdentifier.begin(), rIdentifier.end(), std::back_inserter(upperCaseIdentifier), (int (*)(int)) toupper);

    // Get enum
IdentifierEnum    = Constitutive::GetConstitutiveVariableFromString(upperCaseIdentifier);

    ConstitutiveLawSetParameterBool(rIdent, IdentifierEnum, rValue);
}

//! @brief ... gets a variable of the constitutive law which is selected by a string
//! @param rIdent ... constitutive law identifier
//! @param rIdentifier ... String to identify the requested variable
//! @return ... value of the requested variable
double NuTo::StructureBase::ConstitutiveLawGetParameterDouble(int rIdent, const string &rIdentifier) const
{
    std::string upperCaseIdentifier;
    NuTo::Constitutive::eConstitutiveParameter IdentifierEnum;

    // transform string to uppercase
    std::transform(rIdentifier.begin(), rIdentifier.end(), std::back_inserter(upperCaseIdentifier), (int (*)(int)) toupper);

    // Get enum
IdentifierEnum    = Constitutive::GetConstitutiveVariableFromString(upperCaseIdentifier);

    return ConstitutiveLawGetParameterDouble(rIdent, IdentifierEnum);
}

//! @brief ... sets a variable of the constitutive law which is selected by a sting
//! @param rIdent ... constitutive law identifier
//! @param rIdentifier ... String to identify the requested variable
//! @param rValue ... new value for requested variable
void NuTo::StructureBase::ConstitutiveLawSetParameterDouble(int rIdent, const string &rIdentifier, double rValue)
{
    std::string upperCaseIdentifier;
    NuTo::Constitutive::eConstitutiveParameter IdentifierEnum;

    // transform string to uppercase
    std::transform(rIdentifier.begin(), rIdentifier.end(), std::back_inserter(upperCaseIdentifier), (int (*)(int)) toupper);

    // Get enum
IdentifierEnum    = Constitutive::GetConstitutiveVariableFromString(upperCaseIdentifier);

    ConstitutiveLawSetParameterDouble(rIdent, IdentifierEnum, rValue);
}

//! @brief ... gets a variable of the constitutive law which is selected by a string
//! @param rIdent ... constitutive law identifier
//! @param rIdentifier ... String to identify the requested variable
//! @return ... value of the requested variable
NuTo::FullVector<double, Eigen::Dynamic> NuTo::StructureBase::ConstitutiveLawGetParameterFullVectorDouble(int rIdent, const string &rIdentifier) const
{
    std::string upperCaseIdentifier;
    NuTo::Constitutive::eConstitutiveParameter IdentifierEnum;

    // transform string to uppercase
    std::transform(rIdentifier.begin(), rIdentifier.end(), std::back_inserter(upperCaseIdentifier), (int (*)(int)) toupper);

    // Get enum
IdentifierEnum    = Constitutive::GetConstitutiveVariableFromString(upperCaseIdentifier);

    return ConstitutiveLawGetParameterFullVectorDouble(rIdent, IdentifierEnum);
}

//! @brief ... sets a variable of the constitutive law which is selected by a string
//! @param rIdent ... constitutive law identifier
//! @param rIdentifier ... String to identify the requested variable
//! @param rValue ... new value for requested variable
void NuTo::StructureBase::ConstitutiveLawSetParameterFullVectorDouble(int rIdent, const string &rIdentifier, NuTo::FullVector<double, Eigen::Dynamic> rValue)
{
    std::string upperCaseIdentifier;
    NuTo::Constitutive::eConstitutiveParameter IdentifierEnum;

    // transform string to uppercase
    std::transform(rIdentifier.begin(), rIdentifier.end(), std::back_inserter(upperCaseIdentifier), (int (*)(int)) toupper);

    // Get enum
IdentifierEnum    = Constitutive::GetConstitutiveVariableFromString(upperCaseIdentifier);

    ConstitutiveLawSetParameterFullVectorDouble(rIdent, IdentifierEnum, rValue);
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

//! @brief ... get damage law
//! @return ... damage law
NuTo::FullVector<double, Eigen::Dynamic> NuTo::StructureBase::ConstitutiveLawGetDamageLaw(int rIdent) const
{
    NuTo::FullVector<double, Eigen::Dynamic> damageLaw;
    try
    {
        const ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        damageLaw = ConstitutiveLawPtr->GetParameterFullVectorDouble(Constitutive::eConstitutiveParameter::DAMAGE_LAW);
    } catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawGetDamageLaw] error getting damage law.");
        throw e;
    }
    return damageLaw;
}

//! @brief ... set damage law
//! @param rDamageLaw ... damage law <BR>
//! ============================================================================================<BR>
//! rDamageLaw[0] = Constitutive::Constitutive::eDamageLawType::ISOTROPIC_NO_SOFTENING <BR>
//! w(k) = 1 - e_0/k <BR>
//! rDamageLaw[1] = e_0 // strain at elastic limit <BR>
//! ============================================================================================<BR>
//! rDamageLaw[0] = Constitutive::Constitutive::eDamageLawType::ISOTROPIC_LINEAR_SOFTENING <BR>
//! w(k) = e_c/k * (k-e_0) / (e_c-e_0) <BR>
//! rDamageLaw[1] = e_0 // strain at elastic limit <BR>
//! rDamageLaw[2] = e_c // strain at full damage <BR>
//! ============================================================================================<BR>
//! rDamageLaw[0] = Constitutive::Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING <BR>
//! w(k) = 1 - e_0/k exp{ (e_0 - k) / e_f } <BR>
//! rDamageLaw[1] = e_0 // strain at elastic limit <BR>
//! rDamageLaw[2] = e_f // post-peak slope parameter <BR>
void NuTo::StructureBase::ConstitutiveLawSetDamageLaw(int rIdent, const NuTo::FullVector<double, Eigen::Dynamic> rDamageLaw)
{
    try
    {
        ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        ConstitutiveLawPtr->SetParameterFullVectorDouble(Constitutive::eConstitutiveParameter::DAMAGE_LAW, rDamageLaw);
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


//! @brief ... adds a constitutive law to a multi physics model
//! @param rIdentMultiPhysics ... multi physics constitutive law to which the constitutive law should be added
//! @param rIdentConstitutiveLaw ... constitutive law which should be added to the multi physics model
void NuTo::StructureBase::ConstitutiveLawMultiPhysicsAddConstitutiveLaw(int rIdentMultiPhysics, int rIdentConstitutiveLaw)
{
    try
    {
        ConstitutiveBase* MultiPhysicsLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdentMultiPhysics);
        ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdentConstitutiveLaw);

        MultiPhysicsLawPtr->MultiPhysicsAddConstitutiveLaw(ConstitutiveLawPtr);
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawMultiPhysicsAddConstitutiveLaw] error adding constitutive law to multi physics model.");
        throw e;
    }
}


