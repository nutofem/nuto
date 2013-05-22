// $Id$

#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/constitutive/mechanics/LinearElasticEngineeringStress.h"
#include "nuto/mechanics/constitutive/mechanics/MisesPlasticityEngineeringStress.h"
#include "nuto/mechanics/constitutive/mechanics/GradientDamagePlasticityEngineeringStress.h"
#include "nuto/mechanics/constitutive/mechanics/NonlocalDamagePlasticityEngineeringStress.h"
#include "nuto/mechanics/constitutive/mechanics/StrainGradientDamagePlasticityEngineeringStress.h"
#include "nuto/mechanics/constitutive/thermal/LinearHeatFlux.h"

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
    else if (ConstitutiveLawTypeString == "LINEARHEATFLUX")
    {
        ConstitutiveLawType = Constitutive::LINEAR_HEAT_FLUX;
    }
    else if (ConstitutiveLawTypeString == "STRAINGRADIENTDAMAGEPLASTICITYENGINEERINGSTRESS")
    {
        ConstitutiveLawType = Constitutive::STRAIN_GRADIENT_DAMAGE_PLASTICITY_ENGINEERING_STRESS;
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
        case NuTo::Constitutive::STRAIN_GRADIENT_DAMAGE_PLASTICITY_ENGINEERING_STRESS:
            ConstitutiveLawPtr = new NuTo::StrainGradientDamagePlasticityEngineeringStress();
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
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawSetThermalConductivity] error setting density.");
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
        e.AddMessage("[NuTo::StructureBase::ConstitutiveLawGetDensity] error getting density.");
        throw e;
    }
    return thermalConductivity;
}
