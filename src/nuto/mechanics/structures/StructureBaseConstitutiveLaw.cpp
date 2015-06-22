// $Id$

#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/constitutive/mechanics/LinearElasticEngineeringStress.h"
#include "nuto/mechanics/constitutive/mechanics/DamageViscoPlasticityEngineeringStress.h"
#include "nuto/mechanics/constitutive/mechanics/DamageViscoPlasticityHardeningEngineeringStress.h"
#include "nuto/mechanics/constitutive/mechanics/MisesPlasticityEngineeringStress.h"
#include "nuto/mechanics/constitutive/mechanics/GradientDamagePlasticityEngineeringStress.h"
#include "nuto/mechanics/constitutive/mechanics/GradientDamageEngineeringStress.h"
#include "nuto/mechanics/constitutive/mechanics/NonlocalDamagePlasticityEngineeringStress.h"
#include "nuto/mechanics/constitutive/mechanics/StrainGradientDamagePlasticityEngineeringStress.h"
#include "nuto/mechanics/constitutive/thermal/LinearHeatFlux.h"
#include "nuto/mechanics/constitutive/moistureTransport/MoistureTransport.h"

// create a new constitutive law
int NuTo::StructureBase::ConstitutiveLawCreate(const std::string& rType)
{
    // convert section type string to upper case
    std::string ConstitutiveLawTypeString;
    std::transform(rType.begin(), rType.end(), std::back_inserter(ConstitutiveLawTypeString), (int(*)(int)) toupper);

    // get section type from string
    Constitutive::eConstitutiveType ConstitutiveLawType;
    if (ConstitutiveLawTypeString == "LINEARELASTICENGINEERINGSTRESS")
    {
        ConstitutiveLawType = Constitutive::LINEAR_ELASTIC_ENGINEERING_STRESS;
    }
    else if (ConstitutiveLawTypeString == "MISESPLASTICITYENGINEERINGSTRESS")
    {
        ConstitutiveLawType = Constitutive::MISES_PLASTICITY_ENGINEERING_STRESS;
    }
    else if (ConstitutiveLawTypeString == "NONLOCALDAMAGEPLASTICITYENGINEERINGSTRESS")
    {
        ConstitutiveLawType = Constitutive::NONLOCAL_DAMAGE_PLASTICITY_ENGINEERING_STRESS;
    }
    else if (ConstitutiveLawTypeString == "GRADIENTDAMAGEPLASTICITYENGINEERINGSTRESS")
    {
        ConstitutiveLawType = Constitutive::GRADIENT_DAMAGE_PLASTICITY_ENGINEERING_STRESS;
    }
    else if (ConstitutiveLawTypeString == "GRADIENTDAMAGEENGINEERINGSTRESS")
    {
        ConstitutiveLawType = Constitutive::GRADIENT_DAMAGE_ENGINEERING_STRESS;
    }
    else if (ConstitutiveLawTypeString == "LINEARHEATFLUX")
    {
        ConstitutiveLawType = Constitutive::LINEAR_HEAT_FLUX;
    }
    else if (ConstitutiveLawTypeString == "STRAINGRADIENTDAMAGEPLASTICITYENGINEERINGSTRESS")
    {
        ConstitutiveLawType = Constitutive::STRAIN_GRADIENT_DAMAGE_PLASTICITY_ENGINEERING_STRESS;
    }
    else if (ConstitutiveLawTypeString == "DAMAGEVISCOPLASTICITYENGINEERINGSTRESS")
    {
        ConstitutiveLawType = Constitutive::DAMAGE_VISCO_PLASTICITY_ENGINEERING_STRESS;
    }
    else if (ConstitutiveLawTypeString == "DAMAGEVISCOPLASTICITYHARDENINGENGINEERINGSTRESS")
    {
        ConstitutiveLawType = Constitutive::DAMAGE_VISCO_PLASTICITY_HARDENING_ENGINEERING_STRESS;
    }
    else if (ConstitutiveLawTypeString == "MOISTURETRANSPORT")
    {
        ConstitutiveLawType = Constitutive::MOISTURE_TRANSPORT;
    }
    else
    {
        throw NuTo::MechanicsException("[NuTo::StructureBase::ConstitutiveLawCreate] invalid type of constitutive law.");
    }

	//find unused integer id
	int constitutiveNumber(mConstitutiveLawMap.size());
	boost::ptr_map<int,ConstitutiveBase>::iterator it = mConstitutiveLawMap.find(constitutiveNumber);
	while (it!=mConstitutiveLawMap.end())
	{
		constitutiveNumber++;
		it = mConstitutiveLawMap.find(constitutiveNumber);
	}

    this->ConstitutiveLawCreate(constitutiveNumber, ConstitutiveLawType);
    return constitutiveNumber;
}

// create a new constitutive law
void NuTo::StructureBase::ConstitutiveLawCreate(int rIdent, Constitutive::eConstitutiveType rType)
{
    // check if constitutive law identifier exists
    boost::ptr_map<int,ConstitutiveBase>::iterator it = this->mConstitutiveLawMap.find(rIdent);
    if (it == this->mConstitutiveLawMap.end())
    {
        // create new constitutive law
        ConstitutiveBase* ConstitutiveLawPtr;
        switch (rType)
        {
        case NuTo::Constitutive::LINEAR_ELASTIC_ENGINEERING_STRESS:
            ConstitutiveLawPtr = new NuTo::LinearElasticEngineeringStress();
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
         default:
            throw NuTo::MechanicsException("[NuTo::StructureBase::ConstitutiveLawCreate] invalid type of constitutive law.");
        }

        // add section to map (insert does not allow const keys!!!!)
        this->mConstitutiveLawMap.insert(rIdent, ConstitutiveLawPtr);
        if (ConstitutiveLawPtr->HaveTmpStaticData())
            mHaveTmpStaticData = true;
    }
    else
    {
        throw NuTo::MechanicsException("[NuTo::StructureBase::ConstitutiveLawCreate] Constitutive law already exists.");
    }
}

// delete an existing constitutive law
void NuTo::StructureBase::ConstitutiveLawDelete(int rIdent)
{
    // find constitutive law identifier in map
    boost::ptr_map<int,ConstitutiveBase>::iterator it = this->mConstitutiveLawMap.find(rIdent);
    if (it == this->mConstitutiveLawMap.end())
    {
        throw NuTo::MechanicsException("[NuTo::StructureBase::ConstitutiveLawDelete] Constitutive law does not exist.");
    }
    else
    {
        this->mConstitutiveLawMap.erase(it);
    }
}

// get constitutive law pointer from constitutive law identifier
NuTo::ConstitutiveBase* NuTo::StructureBase::ConstitutiveLawGetConstitutiveLawPtr(int rIdent)
{
    boost::ptr_map<int,ConstitutiveBase>::iterator it = this->mConstitutiveLawMap.find(rIdent);
    if (it == this->mConstitutiveLawMap.end())
    {
        throw NuTo::MechanicsException("[NuTo::StructureBase::ConstitutiveLawGetConstitutiveLawPtr] Constitutive law does not exist.");
    }
    return it->second;
}

// get constitutive law pointer from constitutive law identifier
const NuTo::ConstitutiveBase* NuTo::StructureBase::ConstitutiveLawGetConstitutiveLawPtr(int rIdent) const
{
    boost::ptr_map<int,ConstitutiveBase>::const_iterator it = this->mConstitutiveLawMap.find(rIdent);
    if (it == this->mConstitutiveLawMap.end())
    {
        throw NuTo::MechanicsException("[NuTo::StructureBase::ConstitutiveLawGetConstitutiveLawPtr] Constitutive law does not exist.");
    }
    return it->second;
}

// get constitutive law identifier from constitutive law pointer
int NuTo::StructureBase::ConstitutiveLawGetId(const NuTo::ConstitutiveBase* rConstitutiveLawPtr) const
{
    for (boost::ptr_map<int,ConstitutiveBase>::const_iterator it = mConstitutiveLawMap.begin(); it!= mConstitutiveLawMap.end(); it++)
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
    for (boost::ptr_map<int,ConstitutiveBase>::const_iterator it = mConstitutiveLawMap.begin(); it!= mConstitutiveLawMap.end(); it++)
    {
        std::cout << "  Constitutive law: " << it->first << std::endl;
        it->second->Info(rVerboseLevel,mLogger);
    }
}

void NuTo::StructureBase::ConstitutiveLawInfo(int rIdent, unsigned short rVerboseLevel) const
{
    const ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
    std::cout << "  Constitutive law: " << rIdent << std::endl;
    ConstitutiveLawPtr->Info(rVerboseLevel,mLogger);
}

//! @brief ... gets a variable of the constitutive law which is selected by a string
//! @param rIdent ... constitutive law identifier
//! @param rIdentifier ... String to identify the requested variable
//! @return ... value of the requested variable
bool NuTo::StructureBase::ConstitutiveLawGetVariableBool(int rIdent, const std::string &rIdentifier) const
{
    std::string upperCaseIdentifier;
    NuTo::Constitutive::eConstitutiveVariable IdentifierEnum;

    // transform string to uppercase
    std::transform(rIdentifier.begin(), rIdentifier.end(), std::back_inserter(upperCaseIdentifier), (int(*)(int)) toupper);

    // Get enum
    IdentifierEnum = Constitutive::GetConstitutiveVariableFromString(upperCaseIdentifier);

    return ConstitutiveLawGetVariableBool(rIdent,IdentifierEnum);
}

//! @brief ... sets a variable of the constitutive law which is selected by a string
//! @param rIdent ... constitutive law identifier
//! @param rIdentifier ... String to identify the requested variable
//! @param rValue ... new value for requested variable
void NuTo::StructureBase::ConstitutiveLawSetVariableBool(int rIdent, const std::string &rIdentifier, bool rValue)
{
    std::string upperCaseIdentifier;
    NuTo::Constitutive::eConstitutiveVariable IdentifierEnum;

    // transform string to uppercase
    std::transform(rIdentifier.begin(), rIdentifier.end(), std::back_inserter(upperCaseIdentifier), (int(*)(int)) toupper);

    // Get enum
    IdentifierEnum = Constitutive::GetConstitutiveVariableFromString(upperCaseIdentifier);

    ConstitutiveLawSetVariableBool(rIdent,IdentifierEnum,rValue);
}

//! @brief ... gets a variable of the constitutive law which is selected by a string
//! @param rIdent ... constitutive law identifier
//! @param rIdentifier ... String to identify the requested variable
//! @return ... value of the requested variable
double NuTo::StructureBase::ConstitutiveLawGetVariableDouble(int rIdent, const string &rIdentifier) const
{
    std::string upperCaseIdentifier;
    NuTo::Constitutive::eConstitutiveVariable IdentifierEnum;

    // transform string to uppercase
    std::transform(rIdentifier.begin(), rIdentifier.end(), std::back_inserter(upperCaseIdentifier), (int(*)(int)) toupper);

    // Get enum
    IdentifierEnum = Constitutive::GetConstitutiveVariableFromString(upperCaseIdentifier);

    return ConstitutiveLawGetVariableDouble(rIdent,IdentifierEnum);
}


//! @brief ... sets a variable of the constitutive law which is selected by a sting
//! @param rIdent ... constitutive law identifier
//! @param rIdentifier ... String to identify the requested variable
//! @param rValue ... new value for requested variable
void NuTo::StructureBase::ConstitutiveLawSetVariableDouble(int rIdent, const string &rIdentifier, double rValue)
{
    std::string upperCaseIdentifier;
    NuTo::Constitutive::eConstitutiveVariable IdentifierEnum;

    // transform string to uppercase
    std::transform(rIdentifier.begin(), rIdentifier.end(), std::back_inserter(upperCaseIdentifier), (int(*)(int)) toupper);

    // Get enum
    IdentifierEnum = Constitutive::GetConstitutiveVariableFromString(upperCaseIdentifier);

    ConstitutiveLawSetVariableDouble(rIdent,IdentifierEnum,rValue);
}

//! @brief ... gets a variable of the constitutive law which is selected by a string
//! @param rIdent ... constitutive law identifier
//! @param rIdentifier ... String to identify the requested variable
//! @return ... value of the requested variable
NuTo::FullVector<double, Eigen::Dynamic> NuTo::StructureBase::ConstitutiveLawGetVariableFullVectorDouble(int rIdent, const string &rIdentifier) const
{
    std::string upperCaseIdentifier;
    NuTo::Constitutive::eConstitutiveVariable IdentifierEnum;

    // transform string to uppercase
    std::transform(rIdentifier.begin(), rIdentifier.end(), std::back_inserter(upperCaseIdentifier), (int(*)(int)) toupper);

    // Get enum
    IdentifierEnum = Constitutive::GetConstitutiveVariableFromString(upperCaseIdentifier);

    return ConstitutiveLawGetVariableFullVectorDouble(rIdent,IdentifierEnum);
}

//! @brief ... sets a variable of the constitutive law which is selected by a string
//! @param rIdent ... constitutive law identifier
//! @param rIdentifier ... String to identify the requested variable
//! @param rValue ... new value for requested variable
void NuTo::StructureBase::ConstitutiveLawSetVariableFullVectorDouble(int rIdent, const string &rIdentifier, NuTo::FullVector<double, Eigen::Dynamic> rValue)
{
    std::string upperCaseIdentifier;
    NuTo::Constitutive::eConstitutiveVariable IdentifierEnum;

    // transform string to uppercase
    std::transform(rIdentifier.begin(), rIdentifier.end(), std::back_inserter(upperCaseIdentifier), (int(*)(int)) toupper);

    // Get enum
    IdentifierEnum = Constitutive::GetConstitutiveVariableFromString(upperCaseIdentifier);

    ConstitutiveLawSetVariableFullVectorDouble(rIdent,IdentifierEnum,rValue);
}

//! @brief ... gets a variable of the constitutive law which is selected by an enum
//! @param rIdent ... constitutive law identifier
//! @param rIdentifier ... Enum to identify the requested variable
//! @return ... value of the requested variable
bool NuTo::StructureBase::ConstitutiveLawGetVariableBool(int rIdent, NuTo::Constitutive::eConstitutiveVariable rIdentifier) const
{
    bool requestedValue;
    try
    {
        const ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        requestedValue = ConstitutiveLawPtr->GetVariableBool(rIdentifier);
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawGetVariableBool] error getting requested value.");
        throw e;
    }
    return requestedValue;
}


//! @brief ... sets a variable of the constitutive law which is selected by an enum
//! @param rIdent ... constitutive law identifier
//! @param rIdentifier ... Enum to identify the requested variable
//! @param rValue ... new value for requested variable
void NuTo::StructureBase::ConstitutiveLawSetVariableBool(int rIdent, NuTo::Constitutive::eConstitutiveVariable rIdentifier, bool rValue)
{
    try
    {
        ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        ConstitutiveLawPtr->SetVariableBool(rIdentifier,rValue);
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawSetVariableBool] error setting requested value.");
        throw e;
    }
}


//! @brief ... gets a variable of the constitutive law which is selected by an enum
//! @param rIdent ... constitutive law identifier
//! @param rIdentifier ... Enum to identify the requested variable
//! @return ... value of the requested variable
double NuTo::StructureBase::ConstitutiveLawGetVariableDouble(int rIdent, NuTo::Constitutive::eConstitutiveVariable rIdentifier) const
{
    double requestedValue;
    try
    {
        const ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        requestedValue = ConstitutiveLawPtr->GetVariableDouble(rIdentifier);
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawGetVariableDouble] error getting requested value.");
        throw e;
    }
    return requestedValue;
}

//! @brief ... sets a variable of the constitutive law which is selected by an enum
//! @param rIdent ... constitutive law identifier
//! @param rIdentifier ... Enum to identify the requested variable
//! @param rValue ... new value for requested variable
void NuTo::StructureBase::ConstitutiveLawSetVariableDouble(int rIdent, NuTo::Constitutive::eConstitutiveVariable rIdentifier, double rValue)
{
    try
    {
        ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        ConstitutiveLawPtr->SetVariableDouble(rIdentifier,rValue);
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawSetVariableDouble] error setting requested value.");
        throw e;
    }
}

//! @brief ... gets a variable of the constitutive law which is selected by an enum
//! @param rIdent ... constitutive law identifier
//! @param rIdentifier ... Enum to identify the requested variable
//! @return ... value of the requested variable
NuTo::FullVector<double, Eigen::Dynamic> NuTo::StructureBase::ConstitutiveLawGetVariableFullVectorDouble(int rIdent, NuTo::Constitutive::eConstitutiveVariable rIdentifier) const
{
    NuTo::FullVector<double,Eigen::Dynamic> requestedValue;
    try
    {
        const ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        requestedValue = ConstitutiveLawPtr->GetVariableFullVectorDouble(rIdentifier);
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawGetVariableFullVectorDouble] error getting requested value.");
        throw e;
    }
    return requestedValue;
}

//! @brief ... sets a variable of the constitutive law which is selected by an enum
//! @param rIdent ... constitutive law identifier
//! @param rIdentifier ... Enum to identify the requested variable
//! @param rValue ... new value for requested variable
void NuTo::StructureBase::ConstitutiveLawSetVariableFullVectorDouble(int rIdent, NuTo::Constitutive::eConstitutiveVariable rIdentifier, NuTo::FullVector<double, Eigen::Dynamic> rValue)
{
    try
    {
        ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        ConstitutiveLawPtr->SetVariableFullVectorDouble(rIdentifier,rValue);
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawSetVariableFullVectorDouble] error setting requested value.");
        throw e;
    }
}



// set density
void NuTo::StructureBase::ConstitutiveLawSetDensity(int rIdent, double rRho)
{
    try
    {
        ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        ConstitutiveLawPtr->SetDensity(rRho);
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawSetDensity] error setting density.");
        throw e;
    }
}

// get Young's modulus
double NuTo::StructureBase::ConstitutiveLawGetDensity(int rIdent) const
{
    double density;
    try
    {
        const ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        density = ConstitutiveLawPtr->GetDensity();
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawGetDensity] error getting density.");
        throw e;
    }
    return density;
}

// set Young's modulus
void NuTo::StructureBase::ConstitutiveLawSetYoungsModulus(int rIdent, double rE)
{
    try
    {
        ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        ConstitutiveLawPtr->SetYoungsModulus(rE);
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawSetYoungsModulus] error setting Young's modulus.");
        throw e;
    }
}

// get Young's modulus
double NuTo::StructureBase::ConstitutiveLawGetYoungsModulus(int rIdent) const
{
    double youngs_modulus;
    try
    {
        const ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        youngs_modulus = ConstitutiveLawPtr->GetYoungsModulus();
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawGetYoungsModulus] error getting Young's modulus.");
        throw e;
    }
    return youngs_modulus;
}

// set Poisson's ratio
void NuTo::StructureBase::ConstitutiveLawSetPoissonsRatio(int rIdent, double rNu)
{
    try
    {
        ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        ConstitutiveLawPtr->SetPoissonsRatio(rNu);
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawSetPoissonsRatio] error setting Poisson's ratio.");
        throw e;
    }
}

// get Poisson's ratio
double NuTo::StructureBase::ConstitutiveLawGetPoissonsRatio(int rIdent) const
{
    double nu;
    try
    {
        const ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        nu = ConstitutiveLawPtr->GetPoissonsRatio();
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawGetPoissonsRatio] error getting Poisson's ratio.");
        throw e;
    }
    return nu;
}

//! @brief ... get initial yield strength
//! @return ... yield strength
double NuTo::StructureBase::ConstitutiveLawGetInitialYieldStrength(int rIdent) const
{
	try
	{
		const ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
		return ConstitutiveLawPtr->GetInitialYieldStrength();
	}
	catch (NuTo::MechanicsException& e)
	{
		e.AddMessage("[NuTo::StructureBase::ConstitutiveLawGetInitialYieldStrength] error getting initial yield strength.");
		throw e;
	}
}

//! @brief ... set initial yield strength
//! @param rSigma ...  yield strength
void NuTo::StructureBase::ConstitutiveLawSetInitialYieldStrength(int rIdent, double rSigma)
{
    try
    {
        ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        ConstitutiveLawPtr->SetInitialYieldStrength(rSigma);
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawSetInitialYieldStrength] error setting initial yield strength.");
        throw e;
    }
}

//! @brief ... get yield strength for multilinear response
//! @return ... first column: equivalent plastic strain
//! @return ... second column: corresponding yield strength
NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> NuTo::StructureBase::ConstitutiveLawGetYieldStrength(int rIdent) const
{
	try
	{
		const ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
		return ConstitutiveLawPtr->GetYieldStrength();
	}
	catch (NuTo::MechanicsException& e)
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
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawAddYieldStrength] error adding yield strength.");
        throw e;
    }
}

//! @brief ... get initial hardening modulus
//! @return ... hardening modulus
double NuTo::StructureBase::ConstitutiveLawGetInitialHardeningModulus(int rIdent) const
{
	try
	{
		const ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
		return ConstitutiveLawPtr->GetInitialHardeningModulus();
	}
	catch (NuTo::MechanicsException& e)
	{
		e.AddMessage("[NuTo::StructureBase::ConstitutiveLawGetInitialHardeningModulus] error getting initial hardening modulus.");
		throw e;
	}
}

//! @brief ... set hardening modulus
//! @param rH ...  hardening modulus
void NuTo::StructureBase::ConstitutiveLawSetInitialHardeningModulus(int rIdent, double rH)
{
    try
    {
        ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        ConstitutiveLawPtr->SetInitialHardeningModulus(rH);
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawSetInitialHardeningModulus] error setting initial hardening modulus.");
        throw e;
    }
}

//! @brief ... get hardening value
//! @return ... hardening value
double NuTo::StructureBase::ConstitutiveLawGetHardeningValue(int rIdent) const
{
	try
	{
		const ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
		return ConstitutiveLawPtr->GetHardeningValue();
	}
	catch (NuTo::MechanicsException& e)
	{
		e.AddMessage("[NuTo::StructureBase::ConstitutiveLawGetHardeningValue] error getting hardening value.");
		throw e;
	}
}

//! @brief ... set hardening value
//! @param rHardening ...  hardening value
void NuTo::StructureBase::ConstitutiveLawSetHardeningValue(int rIdent, double rHardening)
{
    try
    {
        ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        ConstitutiveLawPtr->SetHardeningValue(rHardening);
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawSetHardeningValue] error setting hardening value.");
        throw e;
    }
}

//! @brief ... get hardening exponent
//! @return ... hardening exponent
double NuTo::StructureBase::ConstitutiveLawGetHardeningExponent(int rIdent) const
{
	try
	{
		const ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
		return ConstitutiveLawPtr->GetHardeningExponent();
	}
	catch (NuTo::MechanicsException& e)
	{
		e.AddMessage("[NuTo::StructureBase::ConstitutiveLawGetHardeningExponent] error getting hardening exponent.");
		throw e;
	}
}

//! @brief ... set hardening exponent
//! @param rHardeningExponent ...  hardening exponent
void NuTo::StructureBase::ConstitutiveLawSetHardeningExponent(int rIdent, double rHardeningExponent)
{
    try
    {
        ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        ConstitutiveLawPtr->SetHardeningExponent(rHardeningExponent);
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawSetHardeningExponent] error setting hardening exponent.");
        throw e;
    }
}

//! @brief ... get hardening modulus for multilinear response
//! @return ... first column: equivalent plastic strain
//! @return ... second column: corresponding hardening modulus
NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> NuTo::StructureBase::ConstitutiveLawGetHardeningModulus(int rIdent) const
{
	try
	{
		const ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
		return ConstitutiveLawPtr->GetHardeningModulus();
	}
	catch (NuTo::MechanicsException& e)
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
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawAddHardeningModulus] error adding hardening modulus.");
        throw e;
    }
}

//! @brief ... get nonlocal radius
//! @return ... nonlocal radius
double NuTo::StructureBase::ConstitutiveLawGetNonlocalRadius(int rIdent) const
{
	try
	{
		const ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
		return ConstitutiveLawPtr->GetNonlocalRadius();
	}
	catch (NuTo::MechanicsException& e)
	{
		e.AddMessage("[NuTo::StructureBase::ConstitutiveLawGetNonlocalRadius] error getting nonlocal radius.");
		throw e;
	}
}

//! @brief ... set nonlocal radius
//! @param rRadius ...  nonlocal radius
void NuTo::StructureBase::ConstitutiveLawSetNonlocalRadius(int rIdent, double rRadius)
{
    try
    {
        ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        ConstitutiveLawPtr->SetNonlocalRadius(rRadius);
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawSetNonlocalRadius] error setting nonlocal radius.");
        throw e;
    }
}

//! @brief ... get nonlocal radius parameter
//! @return ... nonlocal radius parameter
double NuTo::StructureBase::ConstitutiveLawGetNonlocalRadiusParameter(int rIdent) const
{
    try
    {
        const ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        return ConstitutiveLawPtr->GetNonlocalRadiusParameter();
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawGetNonlocalRadiusParameter] error getting nonlocal radius.");
        throw e;
    }
}

//! @brief ... set nonlocal radius parameter
//! @param rRadius ...  nonlocal radius parameter
void NuTo::StructureBase::ConstitutiveLawSetNonlocalRadiusParameter(int rIdent, double rRadiusParameter)
{
    try
    {
        ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        ConstitutiveLawPtr->SetNonlocalRadiusParameter(rRadiusParameter);
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawSetNonlocalRadiusParameter] error setting nonlocal radius.");
        throw e;
    }
}

//! @brief ... get tensile strength
//! @param rTensileStrength ...  tensile strength
double NuTo::StructureBase::ConstitutiveLawGetTensileStrength(int rIdent)
{
	try
	{
		const ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
		return ConstitutiveLawPtr->GetTensileStrength();
	}
	catch (NuTo::MechanicsException& e)
	{
		e.AddMessage("[NuTo::StructureBase::ConstitutiveLawGetTensileStrength] error getting tensile strength.");
		throw e;
	}
}

//! @brief ... set tensile strength
//! @param rTensileStrength ...  tensile strength
void NuTo::StructureBase::ConstitutiveLawSetTensileStrength(int rIdent, double rTensileStrength)
{
    try
    {
        ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        ConstitutiveLawPtr->SetTensileStrength(rTensileStrength);
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawSetTensileStrength] error setting tensile strength.");
        throw e;
    }
}

//! @brief ... get shear strength
//! @param rIdent ...  constitutive model
double NuTo::StructureBase::ConstitutiveLawGetShearStrength(int rIdent)
{
	try
	{
		const ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
		return ConstitutiveLawPtr->GetShearStrength();
	}
	catch (NuTo::MechanicsException& e)
	{
		e.AddMessage("[NuTo::StructureBase::ConstitutiveLawGetshearStrength] error getting shear strength.");
		throw e;
	}
}

//! @brief ... set shear strength
//! @param rShearStrength ...  shear strength
void NuTo::StructureBase::ConstitutiveLawSetShearStrength(int rIdent, double rShearStrength)
{
    try
    {
        ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        ConstitutiveLawPtr->SetShearStrength(rShearStrength);
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawSetShearStrength] error setting shear strength.");
        throw e;
    }
}

//! @brief ... get compressive strength
//! @param rCompressiveStrength ...  compressive strength
double NuTo::StructureBase::ConstitutiveLawGetCompressiveStrength(int rIdent)
{
	try
	{
		const ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
		return ConstitutiveLawPtr->GetCompressiveStrength();
	}
	catch (NuTo::MechanicsException& e)
	{
		e.AddMessage("[NuTo::StructureBase::ConstitutiveLawGetCompressiveStrength] error getting compressive strength.");
		throw e;
	}
}

//! @brief ... set compressive strength
//! @param rCompressiveStrength ...  compressive strength
void NuTo::StructureBase::ConstitutiveLawSetCompressiveStrength(int rIdent, double rCompressiveStrength)
{
    try
    {
        ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        ConstitutiveLawPtr->SetCompressiveStrength(rCompressiveStrength);
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLaSetCompressiveStrength] error setting compressive strength.");
        throw e;
    }
}

//! @brief ... get biaxial compressive strength
//! @param rBiaxialCompressiveStrength ...  biaxial compressive strength
double NuTo::StructureBase::ConstitutiveLawGetBiaxialCompressiveStrength(int rIdent)
{
	try
	{
		const ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
		return ConstitutiveLawPtr->GetBiaxialCompressiveStrength();
	}
	catch (NuTo::MechanicsException& e)
	{
		e.AddMessage("[NuTo::StructureBase::ConstitutiveLawGetBiaxialCompressiveStrength] error getting biaxial compressive strength.");
		throw e;
	}
}

//! @brief ... set biaxial compressive strength
//! @param rBiaxialCompressiveStrength ...  biaxial compressive strength
void NuTo::StructureBase::ConstitutiveLawSetBiaxialCompressiveStrength(int rIdent, double rBiaxialCompressiveStrength)
{
    try
    {
        ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        ConstitutiveLawPtr->SetBiaxialCompressiveStrength(rBiaxialCompressiveStrength);
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawSetBiaxialCompressiveStrength] error setting biaxial compressive strength.");
        throw e;
    }
}

//! @brief ... get fracture energy
//! @param rFractureEnergy ...  fracture energy
double NuTo::StructureBase::ConstitutiveLawGetFractureEnergy(int rIdent)
{
	try
	{
		const ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
		return ConstitutiveLawPtr->GetFractureEnergy();
	}
	catch (NuTo::MechanicsException& e)
	{
		e.AddMessage("[NuTo::StructureBase::ConstitutiveLawGetFractureEnergy] error getting fracture energy.");
		throw e;
	}
}

//! @brief ... set fracture energy
//! @param rFractureEnergy ...  fracture energy
void NuTo::StructureBase::ConstitutiveLawSetFractureEnergy(int rIdent, double rFractureEnergy)
{
    try
    {
        ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        ConstitutiveLawPtr->SetFractureEnergy(rFractureEnergy);
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawSetFractureEnergy] error setting fracture energy.");
        throw e;
    }
}

//! @brief ... get friction coefficient
//! @param rIdent ...  constitutive model
double NuTo::StructureBase::ConstitutiveLawGetFrictionCoefficient(int rIdent)
{
	try
	{
		const ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
		return ConstitutiveLawPtr->GetFrictionCoefficient();
	}
	catch (NuTo::MechanicsException& e)
	{
		e.AddMessage("[NuTo::StructureBase::ConstitutiveLawGetFrictionCoefficient] error getting friction coefficient.");
		throw e;
	}
}

//! @brief ... set fracture energy
//! @param rFractureEnergy ...  friction coefficient
void NuTo::StructureBase::ConstitutiveLawSetFrictionCoefficient(int rIdent, double rFrictionCoefficient)
{
    try
    {
        ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        ConstitutiveLawPtr->SetFrictionCoefficient(rFrictionCoefficient);
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawSetFrictionCoefficient] error setting friction coefficient.");
        throw e;
    }
}

// set heat cpacity
void NuTo::StructureBase::ConstitutiveLawSetHeatCapacity(int rIdent, double rHeatCapacity)
{
    try
    {
        ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        ConstitutiveLawPtr->SetHeatCapacity(rHeatCapacity);
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawSetHeatCapacity] error setting density.");
        throw e;
    }
}

// get heat cpacity
double NuTo::StructureBase::ConstitutiveLawGetheatCapacity(int rIdent) const
{
    double heatCapacity;
    try
    {
        const ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        heatCapacity = ConstitutiveLawPtr->GetHeatCapacity();
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawGetDensity] error getting density.");
        throw e;
    }
    return heatCapacity;
}

// set thermal conductivity
void NuTo::StructureBase::ConstitutiveLawSetThermalConductivity(int rIdent, double rThermalConductivity)
{
    try
    {
        ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        ConstitutiveLawPtr->SetThermalConductivity(rThermalConductivity);
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawSetThermalConductivity] error setting thermal conductivity.");
        throw e;
    }
}

// get thermal conductivity
double NuTo::StructureBase::ConstitutiveLawGetThermalConductivity(int rIdent) const
{
    double thermalConductivity;
    try
    {
        const ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        thermalConductivity = ConstitutiveLawPtr->GetThermalConductivity();
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawGetThermalConductivity] error getting thermal conductivity.");
        throw e;
    }
    return thermalConductivity;
}

// set viscosity
void NuTo::StructureBase::ConstitutiveLawSetViscosity(int rIdent, double rViscosity)
{
    try
    {
        ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        ConstitutiveLawPtr->SetViscosity(rViscosity);
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawSetViscosity] error setting viscosity.");
        throw e;
    }
}

// get viscosity
double NuTo::StructureBase::ConstitutiveLawGetViscosity(int rIdent) const
{
    double Viscosity;
    try
    {
        const ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        Viscosity = ConstitutiveLawPtr->GetViscosity();
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawGetViscosity] error getting viscosity.");
        throw e;
    }
    return Viscosity;
}

// set viscosity exponent
void NuTo::StructureBase::ConstitutiveLawSetViscosityExponent(int rIdent, double rViscosityExponent)
{
    try
    {
        ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        ConstitutiveLawPtr->SetViscosityExponent(rViscosityExponent);
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawSetViscosityExponent] error setting viscosity exponent.");
        throw e;
    }
}

// get viscosity exponent
double NuTo::StructureBase::ConstitutiveLawGetViscosityExponent(int rIdent) const
{
    double ViscosityExponent;
    try
    {
        const ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        ViscosityExponent = ConstitutiveLawPtr->GetViscosityExponent();
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawGetViscosityExponent] error getting viscosity exponent.");
        throw e;
    }
    return ViscosityExponent;
}

// set damage distribution
void NuTo::StructureBase::ConstitutiveLawSetDamageDistribution(int rIdent, double rDamageDistribution)
{
    try
    {
        ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        ConstitutiveLawPtr->SetDamageDistribution(rDamageDistribution);
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawSetDamageDistribution] error setting damage distribution.");
        throw e;
    }
}

// get damage distribution
double NuTo::StructureBase::ConstitutiveLawGetDamageDistribution(int rIdent) const
{
    double DamageDistribution;
    try
    {
        const ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        DamageDistribution = ConstitutiveLawPtr->GetDamageDistribution();
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawGetDamageDistribution] error getting damage distribution.");
        throw e;
    }
    return DamageDistribution;
}

//! @brief ... get damage law
//! @return ... damage law
NuTo::FullVector<double, Eigen::Dynamic> NuTo::StructureBase::ConstitutiveLawGetDamageLaw(int rIdent) const
{
    NuTo::FullVector<double, Eigen::Dynamic> damageLaw;
    try
    {
        const ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        damageLaw = ConstitutiveLawPtr->GetDamageLaw();
    }
    catch (NuTo::MechanicsException& e)
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
        ConstitutiveLawPtr->SetDamageLaw(rDamageLaw);
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawSetDamageLaw] error setting damage law.");
        throw e;
    }
}

//****************************
// set viscosity
void NuTo::StructureBase::ConstitutiveLawSetViscoplasticYieldSurfaceOffset(int rIdent, double rViscoplasticYieldSurfaceOffset)
{
    try
    {
        ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        ConstitutiveLawPtr->SetViscoplasticYieldSurfaceOffset(rViscoplasticYieldSurfaceOffset);
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawSetViscoplasticYieldSurfaceOffset] error setting viscoplastic yield surface offset.");
        throw e;
    }
}




// get viscosity
double NuTo::StructureBase::ConstitutiveLawGetViscoplasticYieldSurfaceOffset(int rIdent) const
{
    double ViscoplasticYieldSurfaceOffset;
    try
    {
        const ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        ViscoplasticYieldSurfaceOffset = ConstitutiveLawPtr->GetViscoplasticYieldSurfaceOffset();
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawGetViscoplasticYieldSurfaceOffset] error getting viscoplastic yield surface offset.");
        throw e;
    }
    return ViscoplasticYieldSurfaceOffset;
}

//! @brief ... gets the equilibrium water volume fraction depend on the relative humidity
//! @param rIdent ... constitutive law identifier
//! @param rRelativeHumidity ... relative humidity
//! @param rCoeffs ... polynomial coefficients of the sorption curve
//! @return ... equilibrium water volume fraction
double NuTo::StructureBase::ConstitutiveLawGetEquilibriumWaterVolumeFraction(int rIdent, double rRelativeHumidity, NuTo::FullVector<double,Eigen::Dynamic> rCoeffs) const
{
    double EquilibriumWaterVolumeFraction = 0.0;
    try
    {
        const ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        EquilibriumWaterVolumeFraction = ConstitutiveLawPtr->GetEquilibriumWaterVolumeFraction(rRelativeHumidity,rCoeffs);
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawGetEquilibriumWaterVolumeFraction] error getting the equilibrium water volume fraction.");
        throw e;
    }
    return EquilibriumWaterVolumeFraction;
}

// set fatigue flag
void NuTo::StructureBase::ConstitutiveLawSetFatigueExtrapolation(int rIdent, bool rFatigueExtrapolation)
{
	try
	{
		ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
		if (ConstitutiveLawPtr->GetType() == NuTo::Constitutive::DAMAGE_VISCO_PLASTICITY_HARDENING_ENGINEERING_STRESS) {
			ConstitutiveLawPtr->SetFatigueExtrapolation(rFatigueExtrapolation);
		} else {
			std::cout << "[NuTo::StructureBase::ConstitutiveLawSetFatigueExtrapolation] the constitutive law " << ConstitutiveLawPtr->GetType() << " is currently not implemented for fatigue extrapolation." << std::endl;
			throw MechanicsException(" ");
		}
	}
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawSetFatigueExtrapolation] error setting fatigue extrapolation flag.");
        throw e;
    }
}
// get fatigue flag
bool NuTo::StructureBase::ConstitutiveLawGetFatigueExtrapolation(int rIdent) const
{
	bool FatigueExtrapolation;
	try
	{
		const ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
		if (ConstitutiveLawPtr->GetType() == NuTo::Constitutive::DAMAGE_VISCO_PLASTICITY_HARDENING_ENGINEERING_STRESS) {
			FatigueExtrapolation = ConstitutiveLawPtr->GetFatigueExtrapolation();
		} else {
			std::cout << "[NuTo::StructureBase::ConstitutiveLawGetFatigueExtrapolation] the constitutive law " << ConstitutiveLawPtr->GetType() << " is currently not implemented for fatigue extrapolation." << std::endl;
			throw MechanicsException(" ");
		}
	}
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawGetFatigueExtrapolation] error getting fatigue extrapolation flag.");
        throw e;
    }
    return FatigueExtrapolation;
}

