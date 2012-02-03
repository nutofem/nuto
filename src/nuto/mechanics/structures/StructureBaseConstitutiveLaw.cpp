// $Id$

#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/constitutive/mechanics/ConstitutiveLatticeConcrete.h"
#include "nuto/mechanics/constitutive/mechanics/LinearElastic.h"
#include "nuto/mechanics/constitutive/mechanics/ConstitutiveMisesPlasticity.h"
#include "nuto/mechanics/constitutive/mechanics/Multiscale.h"
#include "nuto/mechanics/constitutive/mechanics/NonlocalDamagePlasticity.h"

// create a new constitutive law
int NuTo::StructureBase::ConstitutiveLawCreate(const std::string& rType)
{
    // convert section type string to upper case
    std::string ConstitutiveLawTypeString;
    std::transform(rType.begin(), rType.end(), std::back_inserter(ConstitutiveLawTypeString), (int(*)(int)) toupper);

    // get section type from string
    Constitutive::eConstitutiveType ConstitutiveLawType;
    if (ConstitutiveLawTypeString == "LINEARELASTIC")
    {
        ConstitutiveLawType = Constitutive::LINEAR_ELASTIC;
    }
    else if (ConstitutiveLawTypeString == "MISESPLASTICITY")
    {
        ConstitutiveLawType = Constitutive::MISES_PLASTICITY;
    }
    else if (ConstitutiveLawTypeString == "NONLOCALDAMAGEPLASTICITY")
    {
        ConstitutiveLawType = Constitutive::NONLOCAL_DAMAGE_PLASTICITY;
    }
    else if (ConstitutiveLawTypeString == "MULTISCALE")
    {
        ConstitutiveLawType = Constitutive::MULTISCALE;
    }
    else if (ConstitutiveLawTypeString == "LATTICECONCRETE")
    {
        ConstitutiveLawType = Constitutive::LATTICE_CONCRETE;
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
        case NuTo::Constitutive::LINEAR_ELASTIC:
            ConstitutiveLawPtr = new NuTo::LinearElastic();
            break;
        case NuTo::Constitutive::MISES_PLASTICITY:
            ConstitutiveLawPtr = new NuTo::ConstitutiveMisesPlasticity();
            break;
        case NuTo::Constitutive::NONLOCAL_DAMAGE_PLASTICITY:
            ConstitutiveLawPtr = new NuTo::NonlocalDamagePlasticity();
            break;
        case NuTo::Constitutive::MULTISCALE:
            ConstitutiveLawPtr = new NuTo::Multiscale();
            break;
        case NuTo::Constitutive::LATTICE_CONCRETE:
            ConstitutiveLawPtr = new NuTo::ConstitutiveLatticeConcrete();
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
NuTo::FullMatrix<double> NuTo::StructureBase::ConstitutiveLawGetYieldStrength(int rIdent) const
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

//! @brief ... get hardening modulus for multilinear response
//! @return ... first column: equivalent plastic strain
//! @return ... second column: corresponding hardening modulus
NuTo::FullMatrix<double> NuTo::StructureBase::ConstitutiveLawGetHardeningModulus(int rIdent) const
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

//! @brief ... get elastic stiffness
//! @param rFractureEnergy ...  fracture energy
NuTo::FullMatrix<double> NuTo::StructureBase::ConstitutiveLawGetElasticStiffness(int rIdent)
{
	try
	{
		const ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
		return ConstitutiveLawPtr->GetElasticStiffness();
	}
	catch (NuTo::MechanicsException& e)
	{
		e.AddMessage("[NuTo::StructureBase::ConstitutiveLawGetElasticStiffness] error getting elastic stiffness.");
		throw e;
	}
}

//! @brief ... set fracture energy
//! @param rFractureEnergy ...  fracture energy
void NuTo::StructureBase::ConstitutiveLawSetElasticStiffness(int rIdent, NuTo::FullMatrix<double> rElasticStiffness)
{
    try
    {
        ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        ConstitutiveLawPtr->SetElasticStiffness(rElasticStiffness);
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawSetElasticStiffness] error setting elastic stiffness.");
        throw e;
    }
}

//! @brief ... get elastic stiffness
//! @param rFractureEnergy ...  fracture energy
std::string NuTo::StructureBase::ConstitutiveLawGetMultiscaleFile(int rIdent)
{
	try
	{
		const ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
		return ConstitutiveLawPtr->GetMultiscaleFile();
	}
	catch (NuTo::MechanicsException& e)
	{
		e.AddMessage("[NuTo::StructureBase::ConstitutiveLawGetMultiscaleFile] error getting multiscale file.");
		throw e;
	}
}

//! @brief ... set fracture energy
//! @param rFractureEnergy ...  fracture energy
void NuTo::StructureBase::ConstitutiveLawSetMultiscaleFile(int rIdent, std::string rFileName)
{
    try
    {
        ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        ConstitutiveLawPtr->SetMultiscaleFile(rFileName);
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawSetMultiscaleFile] error setting multiscale file.");
        throw e;
    }
}

//! @brief ... get crack transition radius
//! @param rIdent ...  identifier
double NuTo::StructureBase::ConstitutiveLawGetCrackTransitionRadius(int rIdent)
{
	try
	{
		const ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
		return ConstitutiveLawPtr->GetCrackTransitionRadius();
	}
	catch (NuTo::MechanicsException& e)
	{
		e.AddMessage("[NuTo::StructureBase::ConstitutiveLawGetCrackTransitionRadius] error getting crack transition radius.");
		throw e;
	}
}

//! @brief ... set crack transition radius
//! @param rCrackTransitionRadius ...  fracture energy
void NuTo::StructureBase::ConstitutiveLawSetCrackTransitionRadius(int rIdent, double rCrackTransitionRadius)
{
    try
    {
        ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        ConstitutiveLawPtr->SetCrackTransitionRadius(rCrackTransitionRadius);
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawSetCrackTransitionRadius] error setting crack transition radius.");
        throw e;
    }
}

//! @brief ... get scaling factor for the crack angle
//! @param rIdent ...  identifier
double NuTo::StructureBase::ConstitutiveLawGetScalingFactorCrackAngle(int rIdent)
{
	try
	{
		const ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
		return ConstitutiveLawPtr->GetScalingFactorCrackAngle();
	}
	catch (NuTo::MechanicsException& e)
	{
		e.AddMessage("[NuTo::StructureBase::ConstitutiveLawGetScalingFactorCrackAngle] error getting scaling factor for crack angle.");
		throw e;
	}
}

//! @brief ... set scaling factor for the crack angle
//! @param rScalingFactorCrackAngle ...  scaling factor
void NuTo::StructureBase::ConstitutiveLawSetScalingFactorCrackAngle(int rIdent, double rScalingFactorCrackAngle)
{
    try
    {
        ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        ConstitutiveLawPtr->SetScalingFactorCrackAngle(rScalingFactorCrackAngle);
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawSetScalingFactorCrackAngle] error setting scaling factor for crack angle.");
        throw e;
    }
}

//! @brief ... get scaling factor for the crack opening
//! @param rIdent ...  identifier
double NuTo::StructureBase::ConstitutiveLawGetScalingFactorCrackOpening(int rIdent)
{
	try
	{
		const ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
		return ConstitutiveLawPtr->GetScalingFactorCrackOpening();
	}
	catch (NuTo::MechanicsException& e)
	{
		e.AddMessage("[NuTo::StructureBase::ConstitutiveLawGetScalingFactorCrackOpening] error getting scaling factor for crack opening.");
		throw e;
	}
}

//! @brief ... set scaling factor for the crack opening
//! @param rScalingFactorCrackAngle ...  scaling factor
void NuTo::StructureBase::ConstitutiveLawSetScalingFactorCrackOpening(int rIdent, double rScalingFactorCrackOpening)
{
    try
    {
        ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        ConstitutiveLawPtr->SetScalingFactorCrackOpening(rScalingFactorCrackOpening);
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawSetScalingFactorCrackOpening] error setting scaling factor for crack opening.");
        throw e;
    }
}

//! @brief ... get scaling factor for the total strain
//! @param rIdent ...  identifier
double NuTo::StructureBase::ConstitutiveLawGetScalingFactorEpsilon(int rIdent)
{
	try
	{
		const ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
		return ConstitutiveLawPtr->GetScalingFactorEpsilon();
	}
	catch (NuTo::MechanicsException& e)
	{
		e.AddMessage("[NuTo::StructureBase::ConstitutiveLawGetScalingFactorEpsilon] error getting scaling factor for strain.");
		throw e;
	}
}

//! @brief ... set PenaltyStiffnessCrackAngle
//! @param rScalingFactorCrackAngle ...  scaling factor
void NuTo::StructureBase::ConstitutiveLawSetScalingFactorEpsilon(int rIdent, double rScalingFactorEpsilon)
{
    try
    {
        ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        ConstitutiveLawPtr->SetScalingFactorEpsilon(rScalingFactorEpsilon);
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawSetScalingFactorEpsilon] error setting scaling factor for strain.");
        throw e;
    }
}

//! @brief ... get result directory for fine scale models in multiscale simulation
//! @param rIdent ...  identifier
std::string NuTo::StructureBase::ConstitutiveLawGetResultDirectory(int rIdent)
{
	try
	{
		const ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
		return ConstitutiveLawPtr->GetResultDirectory();
	}
	catch (NuTo::MechanicsException& e)
	{
		e.AddMessage("[NuTo::StructureBase::ConstitutiveLawGetResultDirectory] error getting result directory.");
		throw e;
	}
}

//! @brief ... set ResultDirectory for fine scale models in multiscale simulation
//! @param rResultDirectory ...  ResultDirectory
void NuTo::StructureBase::ConstitutiveLawSetResultDirectory(int rIdent, std::string rResultDirectory)
{
    try
    {
        ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        ConstitutiveLawPtr->SetResultDirectory(rResultDirectory);
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawSetResultDirectory] error setting result directory.");
        throw e;
    }

}

//! @brief ... get LoadStepMacro for fine scale models in multiscale simulation
//! @param rIdent ...  identifier
int NuTo::StructureBase::ConstitutiveLawGetLoadStepMacro(int rIdent)
{
	try
	{
		const ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
		return ConstitutiveLawPtr->GetLoadStepMacro();
	}
	catch (NuTo::MechanicsException& e)
	{
		e.AddMessage("[NuTo::StructureBase::ConstitutiveLawGetResultDirectory] error getting LoadStepMacro.");
		throw e;
	}
}

//! @brief ... set LoadStepMacro for fine scale models in multiscale simulation
//! @param rLoadStepMacro ...  LoadStepMacro
void NuTo::StructureBase::ConstitutiveLawSetLoadStepMacro(int rIdent, int rLoadStepMacro)
{
    try
    {
        ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        ConstitutiveLawPtr->SetLoadStepMacro(rLoadStepMacro);
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawSetLoadStepMacro] error setting LoadStepMacro.");
        throw e;
    }
}

//! @brief ... get if the fine scale model is to be used with the linear elastic periodic boundary shape functions
//! @return rUseAdditionalPeriodicShapeFunctions
bool NuTo::StructureBase::ConstitutiveLawGetUseAdditionalPeriodicShapeFunctions(int rIdent)const
{
	try
	{
		const ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
		return ConstitutiveLawPtr->GetUseAdditionalPeriodicShapeFunctions();
	}
	catch (NuTo::MechanicsException& e)
	{
		e.AddMessage("[NuTo::StructureBase::ConstitutiveLawGetUseAdditionalPeriodicShapeFunctions] error getting UseAdditionalPeriodicShapeFunctions.");
		throw e;
	}
}

//! @brief ... set if the fine scale model is to be used with the linear elastic periodic boundary shape functions
//! @param rUseAdditionalPeriodicShapeFunctions ...  true or false
void NuTo::StructureBase::ConstitutiveLawSetUseAdditionalPeriodicShapeFunctions(int rIdent, bool rUseAdditionalPeriodicShapeFunctions)
{
    try
    {
        ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        ConstitutiveLawPtr->SetUseAdditionalPeriodicShapeFunctions(rUseAdditionalPeriodicShapeFunctions);
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawSetLoadStepMacro] error setting UseAdditionalPeriodicShapeFunctions.");
        throw e;
    }
}

//! @brief ... get the treshold for crack initiation (transistion from a single fine scale model to a combined cracked/uncracked model)
//! @return treshold
double NuTo::StructureBase::ConstitutiveLawGetDamageTresholdCrackInitiation(int rIdent)const
{
	try
	{
		const ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
		return ConstitutiveLawPtr->GetDamageTresholdCrackInitiation();
	}
	catch (NuTo::MechanicsException& e)
	{
		e.AddMessage("[NuTo::StructureBase::ConstitutiveLawGetDamageTresholdCrackInitiation] error getting DamageTresholdCrackInitiation.");
		throw e;
	}

}

//! @brief ... set the treshold for crack initiation (transistion from a single fine scale model to a combined cracked/uncracked model)
//! @param rDamageTresholdCrackInitiation ...  treshold
void NuTo::StructureBase::ConstitutiveLawSetDamageTresholdCrackInitiation(int rIdent, double rDamageTresholdCrackInitiation)
{
    try
    {
        ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        ConstitutiveLawPtr->SetDamageTresholdCrackInitiation(rDamageTresholdCrackInitiation);
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawSetDamageTresholdCrackInitiation] error setting DamageTresholdCrackInitiation.");
        throw e;
    }
}

//! @brief ... get the number of possible crack angles that are checked when the crack is inserted
//! @return number of crack angles
int NuTo::StructureBase::ConstitutiveLawGetNumPossibleCrackAngles(int rIdent)const
{
	try
	{
		const ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
		return ConstitutiveLawPtr->GetNumPossibleCrackAngles();
	}
	catch (NuTo::MechanicsException& e)
	{
		e.AddMessage("[NuTo::StructureBase::ConstitutiveLawGetNumPossibleCrackAngles] error getting number of possible crack angles.");
		throw e;
	}

}

//! @brief ... set the number of possible crack angles that are checked when the crack is inserted
//! @param rNumPossibleCrackAngles ...  number of crack angles
void NuTo::StructureBase::ConstitutiveLawSetNumPossibleCrackAngles(int rIdent, int rNumPossibleCrackAngles)
{
    try
    {
        ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        ConstitutiveLawPtr->SetNumPossibleCrackAngles(rNumPossibleCrackAngles);
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawSetNumPossibleCrackAngles] error setting number of possible crack angles.");
        throw e;
    }
}

//! @brief ... get the number of possible crack shifts that are checked when the crack is inserted
//! @return number of crack angles
int NuTo::StructureBase::ConstitutiveLawGetNumPossibleCrackShifts(int rIdent)const
{
	try
	{
		const ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
		return ConstitutiveLawPtr->GetNumPossibleCrackShifts();
	}
	catch (NuTo::MechanicsException& e)
	{
		e.AddMessage("[NuTo::StructureBase::ConstitutiveLawGetNumPossibleCrackShifts] error getting number of possible crack shifts.");
		throw e;
	}

}

//! @brief ... set the number of possible crack shifts that are checked when the crack is inserted
//! @param rNumPossibleCrackShifts ...  number of crack shifts
void NuTo::StructureBase::ConstitutiveLawSetNumPossibleCrackShifts(int rIdent, int rNumPossibleCrackShifts)
{
    try
    {
        ConstitutiveBase* ConstitutiveLawPtr = this->ConstitutiveLawGetConstitutiveLawPtr(rIdent);
        ConstitutiveLawPtr->SetNumPossibleCrackShifts(rNumPossibleCrackShifts);
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawSetNumPossibleCrackShifts] error setting number of possible crack shifts.");
        throw e;
    }
}
