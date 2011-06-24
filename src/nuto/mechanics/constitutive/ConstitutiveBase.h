// $Id$

#ifndef CONSTITUTIVEBASE_H_
#define CONSTITUTIVEBASE_H_

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include <string>

#include "nuto/mechanics/elements/ElementEnum.h"
#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"

namespace NuTo
{
// forward declarations
class ConstitutiveEngineeringStressStrain;
class ConstitutiveTangentLocal1x1;
class ConstitutiveTangentLocal3x3;
class ConstitutiveTangentLocal6x6;
class DeformationGradient1D;
class DeformationGradient2D;
class DeformationGradient3D;
class ElementBase;
class EngineeringStrain3D;
class EngineeringStress1D;
class EngineeringStress2D;
class EngineeringStress3D;
template<class T>
class FullMatrix;
class Logger;
class SecondPiolaKirchhoffStress1D;
class SecondPiolaKirchhoffStress2D;
class SecondPiolaKirchhoffStress3D;
class StructureBase;

//! @brief ... base class for the constitutive relationship, e.g. material laws
//! @author Stefan Eckardt, ISM
//! @date November 2009
class ConstitutiveBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    //! @brief ... constructor
    ConstitutiveBase();

    //! @brief ... constructor
    virtual ~ConstitutiveBase()
    {}

    // parameters /////////////////////////////////////////////////////////////
    //! @brief ... get density
    //! @return ... density
    virtual double GetDensity() const;

    //! @brief ... set density
    //! @param rRho ... density
    virtual void SetDensity(double rRho);

    //! @brief ... get Young's modulus
    //! @return ... Young's modulus
    virtual double GetYoungsModulus() const;

    //! @brief ... get factor to modify Youngs modulus (using random fields)
    //! @param rElement ...  element
    //! @param rIp ...  integration point
    double GetRanfieldFactorYoungsModulus(const ElementBase* rElement,int rIp) const;

    //! @brief ... set Young's modulus
    //! @param rE ... Young's modulus
    virtual void SetYoungsModulus(double rE);

    //! @brief ... get Poisson's ratio
    //! @return ... Poisson's ratio
    virtual double GetPoissonsRatio() const;

    //! @brief ... set Poisson's ratio
    //! @param rNu ... Poisson's ratio
    virtual void SetPoissonsRatio(double rNu);

    //! @brief ... get factor to modify Poisson's ratio (using random fields)
    //! @param rElement ...  element
    //! @param rIp ...  integration point
    double GetRanfieldFactorPoissonsRatio(const ElementBase* rElement,int rIp) const;

    //! @brief ... get initial yield strength
    //! @return ... yield strength
    virtual double GetInitialYieldStrength() const;

    //! @brief ... set initial yield strength
    //! @param rSigma ...  yield strength
    virtual void SetInitialYieldStrength(double rSigma);

    //! @brief ... get yield strength for multilinear response
    //! @return ... first column: equivalent plastic strain
    //! @return ... second column: corresponding yield strength
    virtual NuTo::FullMatrix<double> GetYieldStrength() const;

    //! @brief ... add yield strength
    //! @param rEpsilon ...  equivalent plastic strain
    //! @param rSigma ...  yield strength
    virtual void AddYieldStrength(double rEpsilon, double rSigma);

    //! @brief ... get factor to modify yield strength (using random fields)
    //! @param rElement ...  element
    //! @param rIp ...  integration point
    double GetRanfieldFactorYieldStrength(const ElementBase* rElement,int rIp) const;

    //! @brief ... get initial hardening modulus
    //! @return ... hardening modulus
    virtual double GetInitialHardeningModulus() const;

    //! @brief ... set hardening modulus
    //! @param rH ...  hardening modulus
    virtual void SetInitialHardeningModulus(double rH);

    //! @brief ... get hardening modulus for multilinear response
    //! @return ... first column: equivalent plastic strain
    //! @return ... second column: corresponding hardening modulus
    virtual NuTo::FullMatrix<double> GetHardeningModulus() const;

    //! @brief ... add hardening modulus
    //! @param rEpsilon ...  equivalent plastic strain
    //! @param rSigma ...  hardening modulus
    virtual void AddHardeningModulus(double rEpsilon, double rH);

    //! @brief ... get factor to modify hardening modulus (using random fields)
    //! @param rElement ...  element
    //! @param rIp ...  integration point
    double GetRanfieldFactorHardeningModulus(const ElementBase* rElement,int rIp) const;

    //! @brief ... get nonlocal radius
    //! @return ... nonlocal radius
    virtual double GetNonlocalRadius() const;

    //! @brief ... set nonlocal radius
    //! @param rRadius...  nonlocal radius
    virtual void SetNonlocalRadius(double rRadius);

    //! @brief ... get tensile strength
    //! @return ... tensile strength
    virtual double GetTensileStrength() const;

    //! @brief ... set tensile strength
    //! @param rTensileStrength...  tensile strength
    virtual void SetTensileStrength(double rTensileStrength);

    //! @brief ... get compressive strength
    //! @return ... compressive strength
    virtual double GetCompressiveStrength() const;

    //! @brief ... set compressive strength
    //! @param rCompressiveStrength...  compressive strength
    virtual void SetCompressiveStrength(double rCompressiveStrength);

    //! @brief ... get biaxial compressive strength
    //! @return ... biaxial compressive strength
    virtual double GetBiaxialCompressiveStrength() const;

    //! @brief ... set biaxial compressive strength
    //! @param rBiaxialCompressiveStrength...  biaxial compressive strength
    virtual void SetBiaxialCompressiveStrength(double rBiaxialCompressiveStrength);

    //! @brief ... get fracture energy
    //! @return ... fracture energy
    virtual double GetFractureEnergy() const;

    //! @brief ... set fracture energy
    //! @param rFractureEnergy... fracture energy
    virtual void SetFractureEnergy(double rFractureEnergy);

    //! @brief ... set the elastic matrix
    //! @param rElasticStiffness... elastic matrix
    virtual NuTo::FullMatrix<double> GetElasticStiffness()const;

    //! @brief ... set the elastic matrix
    //! @param rElasticStiffness... elastic matrix
    virtual void SetElasticStiffness(NuTo::FullMatrix<double> rElasticStiffness);

    //! @brief ... return the binary file from which the fine scale model is eventually deserialized
    //! @return name of the file
    virtual std::string GetMultiscaleFile()const;

    //! @brief ... set the binary file from which the fine scale model is eventually deserialized
    //! @param rFileName... name of the file
    virtual void SetMultiscaleFile(std::string rFileName);

    //! @brief ... return crack transition radius to smooth the Heaviside function in the multiscale model
    //! @return crack transition radius
    virtual double GetCrackTransitionRadius()const;

    //! @brief ... crack transition radius to smooth the Heaviside function in the multiscale model
    //! @param rCrackTransitionRadius... crack transition radius
    virtual void SetCrackTransitionRadius(double rCrackTransitionRadius);

    //! @brief ... get penalty stiffness crack angle
    //! @return ... penalty stiffness crack angle
    virtual double GetPenaltyStiffnessCrackAngle() const;

    //! @brief ... set PenaltyStiffnessCrackAngle
    //! @param rPenaltyStiffnessCrackAngle...  penalty stiffness crack angle
    virtual void SetPenaltyStiffnessCrackAngle(double rPenaltyStiffnessCrackAngle);

    //! @brief ... get scaling factor for the dofs of the crack angle
    //! @return ... scaling factor
    virtual double GetScalingFactorCrackAngle() const;

    //! @brief ... set scaling factor for the dofs of the crack angle
    //! @param rScalingFactor...  scaling factor
    virtual void SetScalingFactorCrackAngle(double rScalingFactorCrackAngle);

    //! @brief ... get scaling factor for the dofs of the crack opening
    //! @return ... scaling factor
    virtual double GetScalingFactorCrackOpening() const;

    //! @brief ... set scaling factor for the dofs of the crack opening
    //! @param rScalingFactor...  scaling factor
    virtual void SetScalingFactorCrackOpening(double rScalingFactorCrackOpening);

    //! @brief ... get scaling factor for the dofs of total strain
    //! @return ... scaling factor
    virtual double GetScalingFactorEpsilon() const;

    //! @brief ... set scaling factor for the dofs of the total strain
    //! @param rScalingFactor...  scaling factor
    virtual void SetScalingFactorEpsilon(double rScalingFactorEpsilon);

    //! @brief ... get AugmentedLagrangeStiffnessCrackOpening
    //! @return ...AugmentedLagrangeStiffnessCrackOpening
    virtual double GetAugmentedLagrangeStiffnessCrackOpening() const;

    //! @brief ... set AugmentedLagrangeStiffnessCrackOpening
    //! @param rAugmentedLagrangeStiffnessCrackOpening...AugmentedLagrangeStiffnessCrackOpening
    virtual void SetAugmentedLagrangeStiffnessCrackOpening(double rAugmentedLagrangeStiffnessCrackOpening);

    //! @brief ... get ToleranceResidualForce in Newton iteration (for multiscale constitutive model)
    //! @return ...ToleranceResidualForce
    virtual double GetToleranceResidualForce() const;

    //! @brief ... set ToleranceResidualForce in Newton iteration (for multiscale constitutive model)
    //! @param rToleranceResidualForce... ToleranceResidualForce
    virtual void SetToleranceResidualForce(double rToleranceResidualForce);

    //! @brief ... get MaxNumNewtonIterations in Newton iteration (for multiscale constitutive model)
    //! @return ...MaxNumNewtonIterations
    virtual int GetMaxNumNewtonIterations() const;

    //! @brief ... MaxNumNewtonIterations in Newton iteration (for multiscale constitutive model)
    //! @param rMaxNumNewtonIterations...MaxNumNewtonIterations
    virtual void SetMaxNumNewtonIterations(int rMaxNumNewtonIterations);

    //! @brief ... get MaxDeltaLoadFactor in Newton iteration (for multiscale constitutive model)
    //! @return ...MaxDeltaLoadFactor
    virtual double GetMaxDeltaLoadFactor() const;

    //! @brief ... set MaxDeltaLoadFactor in Newton iteration (for multiscale constitutive model)
    //! @param rMaxDeltaLoadFactor...MaxDeltaLoadFactor
    virtual void SetMaxDeltaLoadFactor(double rMaxDeltaLoadFactor);

    //! @brief ... get DecreaseFactor in Newton iteration (for multiscale constitutive model)
    //! @return ...DecreaseFactor
    virtual double GetDecreaseFactor() const;

    //! @brief ... set DecreaseFactor in Newton iteration (for multiscale constitutive model)
    //! @param rDecreaseFactor...DecreaseFactor
    virtual void SetDecreaseFactor(double rDecreaseFactor);

    //! @brief ... get MinNumNewtonIterations in Newton iteration (for multiscale constitutive model)
    //! @return ...MinNumNewtonIterations
    virtual int GetMinNumNewtonIterations() const;

    //! @brief ... set MinNumNewtonIterations in Newton iteration (for multiscale constitutive model)
    //! @param rMinNumNewtonIterations...MinNumNewtonIterations
    virtual void SetMinNumNewtonIterations(int rMinNumNewtonIterations);

    //! @brief ... get IncreaseFactor in Newton iteration (for multiscale constitutive model)
    //! @return ...IncreaseFactor
    virtual double GetIncreaseFactor() const;

    //! @brief ... set IncreaseFactor in Newton iteration (for multiscale constitutive model)
    //! @param rIncreaseFactor...
    virtual void SetIncreaseFactor(double rIncreaseFactor);

    //! @brief ... get MinLoadFactor in Newton iteration (for multiscale constitutive model)
    //! @return ...MinLoadFactor
    virtual double GetMinLoadFactor() const;

    //! @brief ... set MinLoadFactor in Newton iteration (for multiscale constitutive model)
    //! @param rMinLoadFactor...MinLoadFactor
    virtual void SetMinLoadFactor(double rMinLoadFactor);

    //! @brief ... get MinLineSearchFactorFactor in Newton iteration (for multiscale constitutive model)
    //! @return ... MinLineSearchFactorFactor
    virtual double GetMinLineSearchFactor() const;

    //! @brief ... set MinLineSearchFactorFactor in Newton iteration (for multiscale constitutive model)
    //! @param rMinLineSearchFactorFactor...MinLineSearchFactorFactor
    virtual void SetMinLineSearchFactor(double rMinLineSearchFactorFactor);

    //! @brief ... get result directory (for results of finescale in multiscale simulations)
    //! @return ... ResultDirectory
    virtual const std::string& GetResultDirectory() const;

    //! @brief ... set ResultDirectory (for results of finescale in multiscale simulations)
    //! @param rResultDirectory...ResultDirectory
    virtual void SetResultDirectory(const std::string& rResultDirectory);

    //! @brief ... get load step macro
    //! @return ... LoadStepMacro
    virtual int GetLoadStepMacro() const;

    //! @brief ... set LoadStepMacro
    //! @param LoadStepMacro...LoadStepMacro
    virtual void SetLoadStepMacro(int rLoadStepMacro);

    ///////////////////////////////////////////////////////////////////////////

    //! @brief ... get type of constitutive relationship
    //! @return ... type of constitutive relationship
    //! @sa eConstitutiveType
    virtual Constitutive::eConstitutiveType GetType() const = 0;

    //! @brief ... check compatibility between element type and type of constitutive relationship
    //! @param rElementType ... element type
    //! @return ... <B>true</B> if the element is compatible with the constitutive relationship, <B>false</B> otherwise.
    virtual bool CheckElementCompatibility(NuTo::Element::eElementType rElementType) const = 0;

    //! @brief ... returns whether the parameters of the constitutive relationship are valid or not
    //! @return ...  <B>true</B> if all parameters of the constitutive relationship are valid and <B>false</B> otherwise
    inline bool AreParametersValid() const
    {
        return this->mParametersValid;
    }

    //! @brief ... check if all parameters are valid and modify parameter validity flag
    //! @sa mParametersValid
    void SetParametersValid();

    //! @brief ... print information about the object
    //! @param rVerboseLevel ... verbosity of the information
    virtual void Info(unsigned short rVerboseLevel, Logger& rLogger) const;

    //! @brief ... returns true, if a material model has tmp static data (which has to be updated before stress or stiffness are calculated)
    //! @return ... see brief explanation
    virtual bool HaveTmpStaticData() const=0;

    //! @brief ... avoid dynamic cast
    //! @return ... see brief explanation
    virtual const ConstitutiveEngineeringStressStrain* AsConstitutiveEngineeringStressStrain()const;

    //! @brief ... avoid dynamic cast
    //! @return ... see brief explanation
    virtual ConstitutiveEngineeringStressStrain* AsConstitutiveEngineeringStressStrain();


#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION
protected:
    //! @brief ... check parameters of the constitutive relationship
    //! if one check fails, an exception is thrwon
    virtual void CheckParameters() const = 0;

    //! @brief ... flag which is <B>true</B> if all parameters of the constitutive relationship are valid and <B>false</B> otherwise
    bool mParametersValid;


};

}

#endif // CONSTITUTIVEBASE_H_ 
