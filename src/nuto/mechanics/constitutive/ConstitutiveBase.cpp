// $Id$
#include "nuto/mechanics/constitutive/ConstitutiveBase.h"
#include "nuto/math/FullMatrix.h"
#include <iostream>

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/base/Logger.h"
#include "nuto/mechanics/MechanicsException.h"

// constructor
NuTo::ConstitutiveBase::ConstitutiveBase()
{
    this->mParametersValid = false;
}

// set density
void NuTo::ConstitutiveBase::SetDensity(double rRho)
{
    throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::SetDensity] The constitutive relationship does not have a parameter density.");
}

// get density
double NuTo::ConstitutiveBase::GetDensity() const
{
    throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::GetDensity] The constitutive relationship does not have a parameter density.");
}


// set Young's modulus
void NuTo::ConstitutiveBase::SetYoungsModulus(double rE)
{
    throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::SetYoungsModulus] The constitutive relationship does not have a parameter Young's modulus.");
}

// get Young's modulus
double NuTo::ConstitutiveBase::GetYoungsModulus() const
{
    throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::GetYoungsModulus] The constitutive relationship does not have a parameter Young's modulus.");
}

//! @brief ... get factor to modify Young's modulus (using random fields)
//! @param rElement ...  element
//! @param rIp ...  integration point
double NuTo::ConstitutiveBase::GetRanfieldFactorYoungsModulus(const ElementBase* rElement,int rIp) const
{
    return 1;
}

// set Poisson's ratio
void NuTo::ConstitutiveBase::SetPoissonsRatio(double rNu)
{
    throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::SetPoissonsRatio] The constitutive relationship does not have a parameter Poisson's ratio.");
}

// get Poisson's ratio
double NuTo::ConstitutiveBase::GetPoissonsRatio() const
{
    throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::GetPoissonsRatio] The constitutive relationship does not have a parameter Poisson's ratio.");
}

//! @brief ... get factor to modify Poisson's ratio (using random fields)
//! @param rElement ...  element
//! @param rIp ...  integration point
double NuTo::ConstitutiveBase::GetRanfieldFactorPoissonsRatio(const ElementBase* rElement,int rIp) const
{
    return 1;
}

//! @brief ... get initial yield strength
//! @return ... yield strength
double NuTo::ConstitutiveBase::GetInitialYieldStrength() const
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::GetInitialYieldStrength] The constitutive relationship does not have a parameter yield strength.");
}

//! @brief ... set initial yield strength
//! @param rSigma ...  yield strength
void NuTo::ConstitutiveBase::SetInitialYieldStrength(double rSigma)
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::SetInitialYieldStrength] The constitutive relationship does not have a parameter yield strength.");
}

//! @brief ... get yield strength for multilinear response
//! @return ... first column: equivalent plastic strain
//! @return ... second column: corresponding yield strength
NuTo::FullMatrix<double> NuTo::ConstitutiveBase::GetYieldStrength() const
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::GetYieldStrength] The constitutive relationship does not have a parameter yield strength.");
}

//! @brief ... add yield strength
//! @param rEpsilon ...  equivalent plastic strain
//! @param rSigma ...  yield strength
void NuTo::ConstitutiveBase::AddYieldStrength(double rEpsilon, double rSigma)
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::AddYieldStrength] The constitutive relationship does not have a parameter yield strength.");
}

//! @brief ... get factor to modify yield strength (using random fields)
//! @param rElement ...  element
//! @param rIp ...  integration point
double NuTo::ConstitutiveBase::GetRanfieldFactorYieldStrength(const ElementBase* rElement,int rIp) const
{
    return 1;
}

//! @brief ... get initial hardening modulus
//! @return ... hardening modulus
double NuTo::ConstitutiveBase::GetInitialHardeningModulus() const
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::GetHardeningModulus] The constitutive relationship does not have a parameter hardening modulus.");
}

//! @brief ... set initial hardening modulus
//! @param rH ...  hardening modulus
void NuTo::ConstitutiveBase::SetInitialHardeningModulus(double rH)
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::SetHardeningModulus] The constitutive relationship does not have a parameter hardening modulus.");
}

//! @brief ... get hardening modulus for multilinear response
//! @return ... first column: equivalent plastic strain
//! @return ... second column: corresponding hardening modulus
NuTo::FullMatrix<double> NuTo::ConstitutiveBase::GetHardeningModulus() const
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::GetYieldStrength] The constitutive relationship does not have a parameter yield strength.");
}

//! @brief ... add hardening modulus
//! @param rEpsilon ...  equivalent plastic strain
//! @param rSigma ...  hardening modulus
void NuTo::ConstitutiveBase::AddHardeningModulus(double rEpsilon, double rH)
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::AddYieldStrength] The constitutive relationship does not have a parameter yield strength.");
}

//! @brief ... get factor to modify hardening modulus (using random fields)
//! @param rElement ...  element
//! @param rIp ...  integration point
double NuTo::ConstitutiveBase::GetRanfieldFactorHardeningModulus(const ElementBase* rElement,int rIp) const
{
    return 1;
}

//! @brief ... get nonlocal radius
//! @return ... nonlocal radius
double NuTo::ConstitutiveBase::GetNonlocalRadius() const
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::GetNonlocalRadius] The constitutive relationship does not have a nonlocal radius.");
}

//! @brief ... set nonlocal radius
//! @param rRadius ...  nonlocal radius
void NuTo::ConstitutiveBase::SetNonlocalRadius(double rH)
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::SetNonlocalRadius] The constitutive relationship does not have a nonlocal radius.");
}
//! @brief ... get tensile strength
//! @return ... tensile strength
double NuTo::ConstitutiveBase::GetTensileStrength() const
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::GetTensileStrength] The constitutive relationship does not have a tensile strength.");
}

//! @brief ... set tensile strength
//! @param rTensileStrength...  tensile strength
void NuTo::ConstitutiveBase::SetTensileStrength(double rTensileStrength)
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::SetTensileStrength] The constitutive relationship does not have a tensile strength.");
}

//! @brief ... get compressive strength
//! @return ... compressive strength
double NuTo::ConstitutiveBase::GetCompressiveStrength() const
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::GetCompressiveStrength] The constitutive relationship does not have a compressive strength.");
}

//! @brief ... set compressive strength
//! @param rCompressiveStrength...  compressive strength
void NuTo::ConstitutiveBase::SetCompressiveStrength(double rCompressiveStrength)
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::SetCompressiveStrength] The constitutive relationship does not have a compressive strength.");
}

//! @brief ... get biaxial compressive strength
//! @return ... biaxial compressive strength
double NuTo::ConstitutiveBase::GetBiaxialCompressiveStrength() const
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::GetBiaxialCompressiveStrength] The constitutive relationship does not have a biaxial compressive strength.");
}

//! @brief ... set biaxial compressive strength
//! @param rBiaxialCompressiveStrength...  biaxial compressive strength
void NuTo::ConstitutiveBase::SetBiaxialCompressiveStrength(double rBiaxialCompressiveStrength)
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::SetBiaxialCompressiveStrength] The constitutive relationship does not have a biaxial compressive strength.");
}

//! @brief ... get fracture energy
//! @return ... fracture energy
double NuTo::ConstitutiveBase::GetFractureEnergy() const
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::GetFractureEnergy] The constitutive relationship does not have a fracture energy.");
}

//! @brief ... set fracture energy
//! @param rFractureEnergy... fracture energy
void NuTo::ConstitutiveBase::SetFractureEnergy(double rFractureEnergy)
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::SetFractureEnergy] The constitutive relationship does not have a fracture energy.");
}

//! @brief ... set the elastic matrix
//! @param rElasticStiffness... elastic matrix
NuTo::FullMatrix<double> NuTo::ConstitutiveBase::GetElasticStiffness()const
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::GetElasticStiffness] The constitutive relationship does not have an elastic stiffness as material parameter.");
}

//! @brief ... set the elastic matrix
//! @param rElasticStiffness... elastic matrix
void NuTo::ConstitutiveBase::SetElasticStiffness(NuTo::FullMatrix<double> rElasticStiffness)
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::SetElasticStiffness] The constitutive relationship does not have an elastic stiffness as material parameter.");
}

//! @brief ... return the binary file from which the fine scale model is eventually deserialized
//! @return name of the file
std::string NuTo::ConstitutiveBase::GetMultiscaleFile()const
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::GetMultiscaleFile] The constitutive relationship does not have a multiscale file.");
}

//! @brief ... set the binary file from which the fine scale model is eventually deserialized
//! @param rFileName... name of the file
void NuTo::ConstitutiveBase::SetMultiscaleFile(std::string rFileName)
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::SetMultiscaleFile] The constitutive relationship does not have a multiscale file.");
}

//! @brief ... return crack transition radius to smooth the Heaviside function in the multiscale model
//! @return crack transition radius
double NuTo::ConstitutiveBase::GetCrackTransitionRadius()const
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::GetCrackTransitionRadius] The constitutive relationship does not have a crack transition radius.");
}

//! @brief ... crack transition radius to smooth the Heaviside function in the multiscale model
//! @param rCrackTransitionRadius... crack transition radius
void NuTo::ConstitutiveBase::SetCrackTransitionRadius(double rCrackTransitionRadius)
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::SetCrackTransitionRadius] The constitutive relationship does not have a crack transition radius.");
}

//! @brief ... get scaling factor for the dofs of the crack angle
//! @return ... scaling factor
double NuTo::ConstitutiveBase::GetScalingFactorCrackAngle() const
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::GetScalingFactorCrackAngle] The constitutive relationship does not have a scaling factor for the crack angle.");
}

//! @brief ... set scaling factor for the dofs of the crack angle
//! @param rScalingFactor...  scaling factor
void NuTo::ConstitutiveBase::SetScalingFactorCrackAngle(double rScalingFactorCrackAngle)
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::SetScalingFactorCrackAngle] The constitutive relationship does not have a scaling factor for the crack angle.");
}

//! @brief ... get scaling factor for the dofs of the crack opening
//! @return ... scaling factor
double NuTo::ConstitutiveBase::GetScalingFactorCrackOpening() const
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::GetScalingFactorCrackOpening] The constitutive relationship does not have a scaling factor for the crack opening.");
}

//! @brief ... set scaling factor for the dofs of the crack opening
//! @param rScalingFactor...  scaling factor
void NuTo::ConstitutiveBase::SetScalingFactorCrackOpening(double rScalingFactorCrackOpening)
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::SetScalingFactorCrackOpening] The constitutive relationship does not have a scaling factor for the crack opening.");
}

//! @brief ... get scaling factor for the dofs of total strain
//! @return ... scaling factor
double NuTo::ConstitutiveBase::GetScalingFactorEpsilon() const
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::GetScalingFactorEpsilon] The constitutive relationship does not have a scaling factor for the strain.");
}

//! @brief ... set scaling factor for the dofs of the total strain
//! @param rScalingFactor...  scaling factor
void NuTo::ConstitutiveBase::SetScalingFactorEpsilon(double rScalingFactorEpsilon)
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::SetScalingFactorEpsilon] The constitutive relationship does not have a scaling factor for the strain.");
}


//! @brief ... get AugmentedLagrangeStiffnessCrackOpening
//! @return ...AugmentedLagrangeStiffnessCrackOpening
double NuTo::ConstitutiveBase::GetAugmentedLagrangeStiffnessCrackOpening() const
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::GetAugmentedLagrangeStiffnessCrackOpening] The constitutive relationship does not have this parameter.");
}

//! @brief ... set AugmentedLagrangeStiffnessCrackOpening
//! @param rAugmentedLagrangeStiffnessCrackOpening...AugmentedLagrangeStiffnessCrackOpening
void NuTo::ConstitutiveBase::SetAugmentedLagrangeStiffnessCrackOpening(double rAugmentedLagrangeStiffnessCrackOpening)
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::SetAugmentedLagrangeStiffnessCrackOpening] The constitutive relationship does not have this parameter.");
}

//! @brief ... get ToleranceResidualForce in Newton iteration (for multiscale constitutive model)
//! @return ...ToleranceResidualForce
double NuTo::ConstitutiveBase::GetToleranceResidualForce() const
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::GetToleranceResidualForce] The constitutive relationship does not have this parameter.");
}

//! @brief ... set ToleranceResidualForce in Newton iteration (for multiscale constitutive model)
//! @param rToleranceResidualForce... ToleranceResidualForce
void NuTo::ConstitutiveBase::SetToleranceResidualForce(double rToleranceResidualForce)
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::SetToleranceResidualForce] The constitutive relationship does not have this parameter.");
}

//! @brief ... get MaxNumNewtonIterations in Newton iteration (for multiscale constitutive model)
//! @return ...MaxNumNewtonIterations
int NuTo::ConstitutiveBase::GetMaxNumNewtonIterations() const
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::GetMaxNumNewtonIterations] The constitutive relationship does not have this parameter.");
}

//! @brief ... MaxNumNewtonIterations in Newton iteration (for multiscale constitutive model)
//! @param rMaxNumNewtonIterations...MaxNumNewtonIterations
void NuTo::ConstitutiveBase::SetMaxNumNewtonIterations(int rMaxNumNewtonIterations)
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::SetMaxNumNewtonIterations] The constitutive relationship does not have this parameter.");
}

//! @brief ... get MaxDeltaLoadFactor in Newton iteration (for multiscale constitutive model)
//! @return ...MaxDeltaLoadFactor
double NuTo::ConstitutiveBase::GetMaxDeltaLoadFactor() const
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::GetMaxDeltaLoadFactor] The constitutive relationship does not have this parameter.");
}

//! @brief ... set MaxDeltaLoadFactor in Newton iteration (for multiscale constitutive model)
//! @param rMaxDeltaLoadFactor...MaxDeltaLoadFactor
void NuTo::ConstitutiveBase::SetMaxDeltaLoadFactor(double rMaxDeltaLoadFactor)
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::SetMaxDeltaLoadFactor] The constitutive relationship does not have this parameter.");
}

//! @brief ... get DecreaseFactor in Newton iteration (for multiscale constitutive model)
//! @return ...DecreaseFactor
double NuTo::ConstitutiveBase::GetDecreaseFactor() const
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::GetDecreaseFactor] The constitutive relationship does not have this parameter.");
}

//! @brief ... set DecreaseFactor in Newton iteration (for multiscale constitutive model)
//! @param rDecreaseFactor...DecreaseFactor
void NuTo::ConstitutiveBase::SetDecreaseFactor(double rDecreaseFactor)
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::SetDecreaseFactor] The constitutive relationship does not have this parameter.");
}

//! @brief ... get MinNumNewtonIterations in Newton iteration (for multiscale constitutive model)
//! @return ...MinNumNewtonIterations
int NuTo::ConstitutiveBase::GetMinNumNewtonIterations() const
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::GetMinNumNewtonIterations] The constitutive relationship does not have this parameter.");
}

//! @brief ... set MinNumNewtonIterations in Newton iteration (for multiscale constitutive model)
//! @param rMinNumNewtonIterations...MinNumNewtonIterations
void NuTo::ConstitutiveBase::SetMinNumNewtonIterations(int rMinNumNewtonIterations)
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::SetMinNumNewtonIterations] The constitutive relationship does not have this parameter.");
}

//! @brief ... get IncreaseFactor in Newton iteration (for multiscale constitutive model)
//! @return ...IncreaseFactor
double NuTo::ConstitutiveBase::GetIncreaseFactor() const
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::GetIncreaseFactor] The constitutive relationship does not have this parameter.");
}

//! @brief ... set IncreaseFactor in Newton iteration (for multiscale constitutive model)
//! @param rIncreaseFactor...
void NuTo::ConstitutiveBase::SetIncreaseFactor(double rIncreaseFactor)
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::SetIncreaseFactor] The constitutive relationship does not have this parameter.");
}

//! @brief ... get MinLoadFactor in Newton iteration (for multiscale constitutive model)
//! @return ...MinLoadFactor
double NuTo::ConstitutiveBase::GetMinLoadFactor() const
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::GetMinLoadFactor] The constitutive relationship does not have this parameter.");
}

//! @brief ... set MinLoadFactor in Newton iteration (for multiscale constitutive model)
//! @param rMinLoadFactor...MinLoadFactor
void NuTo::ConstitutiveBase::SetMinLoadFactor(double rMinLoadFactor)
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::SetMinLoadFactor] The constitutive relationship does not have this parameter.");
}

//! @brief ... get MinLineSearchFactorFactor in Newton iteration (for multiscale constitutive model)
//! @return ... MinLineSearchFactor
double NuTo::ConstitutiveBase::GetMinLineSearchFactor() const
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::GetMinLineSearchFactor] The constitutive relationship does not have this parameter.");
}

//! @brief ... set MinLineSearchFactorFactor in Newton iteration (for multiscale constitutive model)
//! @param rMinLineSearchFactor...MinLineSearchFactor
void NuTo::ConstitutiveBase::SetMinLineSearchFactor(double rMinLineSearchFactor)
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::SetMinLineSearchFactor] The constitutive relationship does not have this parameter.");
}

//! @brief ... get result directory (for results of finescale in multiscale simulations)
//! @return ... ResultDirectory
const std::string& NuTo::ConstitutiveBase::GetResultDirectory() const
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::GetResultDirectory] The constitutive relationship does not have this parameter.");
}

//! @brief ... set ResultDirectory (for results of finescale in multiscale simulations)
//! @param rResultDirectory...ResultDirectory
void NuTo::ConstitutiveBase::SetResultDirectory(const std::string& rResultDirectory)
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::SetResultDirectory] The constitutive relationship does not have this parameter.");
}

//! @brief ... get load step macro
//! @return ... LoadStepMacro
int NuTo::ConstitutiveBase::GetLoadStepMacro() const
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::GetLoadStepMacro] The constitutive relationship does not have this parameter.");
}

//! @brief ... set LoadStepMacro
//! @param LoadStepMacro...LoadStepMacro
void NuTo::ConstitutiveBase::SetLoadStepMacro(int rLoadStepMacro)
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::GetLoadStepMacro] The constitutive relationship does not have this parameter.");
}

//! @brief ... get if additional periodic shape functions are used
//! @return ... true (periodic) or false (fixed displacements)
bool NuTo::ConstitutiveBase::GetUseAdditionalPeriodicShapeFunctions() const
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::GetUseAddPeriodicShapeFunctions] The constitutive relationship does not have this parameter.");
}

//! @brief ... set to use additional periodic shape functions
//! @param rUseAddPeriodicShapeFunctions...rUseAddPeriodicShapeFunctions
void NuTo::ConstitutiveBase::SetUseAdditionalPeriodicShapeFunctions(bool rUseAddPeriodicShapeFunctions)
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::SetUseAddPeriodicShapeFunctions] The constitutive relationship does not have this parameter.");
}

//! @brief ... get threshold for crack initiation based on the maximum damage value within the structure
//! @return ... mDamageTresholdCrackInitiation
double NuTo::ConstitutiveBase::GetDamageTresholdCrackInitiation() const
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::GetDamageTresholdCrackInitiation] The constitutive relationship does not have this parameter.");
}

//! @brief ... set DamageTresholdCrackInitiation
//! @param rDamageTresholdCrackInitiation...DamageTresholdCrackInitiation
void NuTo::ConstitutiveBase::SetDamageTresholdCrackInitiation(double rDamageTresholdCrackInitiation)
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::SetDamageTresholdCrackInitiation] The constitutive relationship does not have this parameter.");
}

//! @brief ... get number of possible crack shifts that are checked when the crack is inserted
//! @return ... NumPossibleCrackAngles
int NuTo::ConstitutiveBase::GetNumPossibleCrackAngles() const
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::GetNumPossibleCrackAngles] The constitutive relationship does not have this parameter.");
}

//! @brief ... set number of possible crack shifts that are checked when the crack is inserted
//! @param rNumPossibleCrackAngles...NumPossibleCrackAngles
void NuTo::ConstitutiveBase::SetNumPossibleCrackAngles(int rNumPossibleCrackAngles)
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::SetNumPossibleCrackAngles] The constitutive relationship does not have this parameter.");
}

//! @brief ... get number of possible crack orientations that are checked when the crack is inserted
//! @return ... mNumPossibleCrackShifts
int NuTo::ConstitutiveBase::GetNumPossibleCrackShifts() const
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::GetNumPossibleCrackShifts] The constitutive relationship does not have this parameter.");
}

//! @brief ... set number of possible crack orientations that are checked when the crack is inserted
//! @param mNumPossibleCrackShifts...mNumPossibleCrackShifts
void NuTo::ConstitutiveBase::SetNumPossibleCrackShifts(int rNumPossibleCrackShifts)
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::SetNumPossibleCrackShifts] The constitutive relationship does not have this parameter.");
}


// modify parameter validity flag
void NuTo::ConstitutiveBase::SetParametersValid()
{
    try
    {
        this->CheckParameters();
    }
    catch (NuTo::MechanicsException)
    {
        this->mParametersValid = false;
        return;
    }
    catch (...)
    {
        throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::SetParametersValid] Unhandled exception");
    }
    this->mParametersValid = true;
}

//! @brief ... avoid dynamic cast
//! @return ... see brief explanation
NuTo::ConstitutiveEngineeringStressStrain* NuTo::ConstitutiveBase::AsConstitutiveEngineeringStressStrain()
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::AsConstitutiveEngineeringStressStrain] Constitutive Law is not of type EngineeringStressStrain.");
}

//! @brief ... avoid dynamic cast
//! @return ... see brief explanation
const NuTo::ConstitutiveEngineeringStressStrain* NuTo::ConstitutiveBase::AsConstitutiveEngineeringStressStrain()const
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::AsConstitutiveEngineeringStressStrain] Constitutive Law is not of type EngineeringStressStrain.");
}

// info routine
void NuTo::ConstitutiveBase::Info(unsigned short rVerboseLevel, Logger& rLogger) const
{
    std::cout << "    parameter validity flag: " << this->mParametersValid << std::endl;
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::ConstitutiveBase::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveBase::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveBase::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveBase::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveBase::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveBase::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstitutiveBase::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize constitutive Base" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_NVP(mParametersValid);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize Constitutive Base" << std::endl;
#endif
}
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::ConstitutiveBase)
#endif // ENABLE_SERIALIZATION
