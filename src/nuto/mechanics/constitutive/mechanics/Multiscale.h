// $Id: Multiscale.h 342 2010-10-18 12:39:08Z arnold2 $
#ifndef CONSTITUTIVEMULTISCALE_H_
#define CONSTITUTIVEMULTISCALE_H_

#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/elements/ElementEnum.h"

#include "nuto/mechanics/constitutive/mechanics/ConstitutiveEngineeringStressStrain.h"

namespace NuTo
{
class StructureMultiscale;
class ConstitutiveStaticDataMultiscale2DPlaneStrain;
//! @author Joerg F. Unger
//! @date Apr 26, 2010
//! @brief ...
class Multiscale : public ConstitutiveEngineeringStressStrain
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    Multiscale();

    //  Engineering strain /////////////////////////////////////
    //! @brief ... calculate engineering plastic strain from deformation gradient in 3D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rEngineeringStrain ... engineering strain
    void GetEngineeringPlasticStrain(const ElementBase* rElement, int rIp,
                                      const DeformationGradient1D& rDeformationGradient, EngineeringStrain3D& rEngineeringPlasticStrain) const;

    //  Engineering strain /////////////////////////////////////
    //! @brief ... calculate engineering plastic strain from deformation gradient in 3D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rEngineeringStrain ... engineering strain
    void GetEngineeringPlasticStrain(const ElementBase* rElement, int rIp,
                                      const DeformationGradient2D& rDeformationGradient, EngineeringStrain3D& rEngineeringPlasticStrain) const;

    //  Engineering strain /////////////////////////////////////
    //! @brief ... calculate engineering plastic strain from deformation gradient in 3D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rEngineeringStrain ... engineering strain
    void GetEngineeringPlasticStrain(const ElementBase* rElement, int rIp,
                                      const DeformationGradient3D& rDeformationGradient, EngineeringStrain3D& rEngineeringPlasticStrain) const;


    // Engineering stress - Engineering strain /////////////////////////////////////
    //! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rEngineeringStress ... engineering stress
    void GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
    		      const DeformationGradient1D& rDeformationGradient, EngineeringStress1D& rEngineeringStress) const;

    // Engineering stress - Engineering strain /////////////////////////////////////
    //! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rEngineeringStress ... engineering stress
    void GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
    		      const DeformationGradient1D& rDeformationGradient, EngineeringStress3D& rEngineeringStress) const;

    // Engineering stress - Engineering strain /////////////////////////////////////
    //! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rEngineeringStress ... engineering stress
    void GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
    		      const DeformationGradient2D& rDeformationGradient, EngineeringStress2D& rEngineeringStress) const;

    // Engineering stress - Engineering strain /////////////////////////////////////
    //! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rEngineeringStress ... engineering stress
    void GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
    		      const DeformationGradient2D& rDeformationGradient, EngineeringStress3D& rEngineeringStress) const;

    // Engineering stress - Engineering strain /////////////////////////////////////
    //! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rCauchyStress ... Cauchy stress
    void GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
    		      const DeformationGradient3D& rDeformationGradient, EngineeringStress3D& rEngineeringStress) const;

    //! @brief ... calculate the tangent (derivative of the Engineering stresses with respect to the engineering strains) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rTangent ... tangent
    void GetTangent_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient1D& rDeformationGradient,
            ConstitutiveTangentBase* rTangent) const;

    //! @brief ... calculate the tangent (derivative of the Engineering stresses with respect to the engineering strains) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rTangent ... tangent
    void GetTangent_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient2D& rDeformationGradient,
            ConstitutiveTangentBase* rTangent) const;

    //! @brief ... calculate the tangent (derivative of the Engineering stresses with respect to the engineering strains) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rTangent ... tangent
    void GetTangent_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient3D& rDeformationGradient,
            ConstitutiveTangentBase* rTangent) const;

    //  Damage /////////////////////////////////////
    //! @brief ... calculate isotropic damage from deformation gradient in 1D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rDamage ... damage variable
    void GetDamage(const ElementBase* rElement, int rIp,
                                      const DeformationGradient1D& rDeformationGradient, double& rDamage) const;

    //  Damage /////////////////////////////////////
    //! @brief ... calculate isotropic damage from deformation gradient in 2D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rDamage ... damage variable
    void GetDamage(const ElementBase* rElement, int rIp,
                                      const DeformationGradient2D& rDeformationGradient, double& rDamage) const;

    //  Damage /////////////////////////////////////
    //! @brief ... calculate isotropic damage from deformation gradient in 3D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rDamage ... damage variable
    void GetDamage(const ElementBase* rElement, int rIp,
                                      const DeformationGradient3D& rDeformationGradient, double& rDamage) const;

    //! @brief ... update static data (history variables) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    void UpdateStaticData_EngineeringStress_EngineeringStrain( ElementBase* rElement, int rIp,
    		const DeformationGradient1D& rDeformationGradient) const;

    //! @brief ... update static data (history variables) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    void UpdateStaticData_EngineeringStress_EngineeringStrain( ElementBase* rElement, int rIp,
    		const DeformationGradient2D& rDeformationGradient) const;

    //! @brief ... update static data (history variables) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    void UpdateStaticData_EngineeringStress_EngineeringStrain( ElementBase* rElement, int rIp,
    		const DeformationGradient3D& rDeformationGradient) const;

    //! @brief ... update tmp static data (history variables) of the constitutive relationship
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    void UpdateTmpStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
            const DeformationGradient1D& rDeformationGradient) const;

    //! @brief ... update tmp static data (history variables) of the constitutive relationship
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    void UpdateTmpStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
            const DeformationGradient2D& rDeformationGradient) const;

    //! @brief ... update tmp static data (history variables) of the constitutive relationship
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    void UpdateTmpStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
            const DeformationGradient3D& rDeformationGradient) const;

    //! @brief ... checks, if a model has to be switched from linear to nonlinear, and then performs the adaption
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    void MultiscaleSwitchToNonlinear(ElementBase* rElement, int rIp,
            const DeformationGradient2D& rDeformationGradient) const;

    //! @brief ... create new static data object for an integration point
    //! @return ... pointer to static data object
    ConstitutiveStaticDataBase* AllocateStaticDataEngineeringStress_EngineeringStrain1D(const ElementBase* rElement) const;

    //! @brief ... create new static data object for an integration point
    //! @return ... pointer to static data object
    ConstitutiveStaticDataBase* AllocateStaticDataEngineeringStress_EngineeringStrain2D(const ElementBase* rElement) const;

    //! @brief ... create new static data object for an integration point
    //! @return ... pointer to static data object
    ConstitutiveStaticDataBase* AllocateStaticDataEngineeringStress_EngineeringStrain3D(const ElementBase* rElement) const;

    //! @brief ... calculate the total energy density
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    double GetTotalEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient1D& rDeformationGradient) const;

    //! @brief ... calculate the total energy density
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    double GetTotalEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient2D& rDeformationGradient) const;

    //! @brief ... calculate the total energy density
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    double GetTotalEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient3D& rDeformationGradient) const;

    //! @brief ... calculate the elastic energy density
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    double GetElasticEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient1D& rDeformationGradient) const;

    //! @brief ... calculate the elastic energy density
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    double GetElasticEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient2D& rDeformationGradient) const;

    //! @brief ... calculate the elastic energy density
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    double GetElasticEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient3D& rDeformationGradient) const;

    //! @brief ... calculates the difference of the elastic strain between the current state and the previous update
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rDeltaElasticEngineeringStrain ... delta elastic engineering strain (return value)
    void GetDeltaElasticEngineeringStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient1D& rDeformationGradient, EngineeringStrain1D& rDeltaElasticEngineeringStrain) const;

    //! @brief ... calculates the difference of the elastic strain between the current state and the previous update
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rDeltaElasticEngineeringStrain ... delta elastic engineering strain (return value)
    void GetDeltaElasticEngineeringStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient2D& rDeformationGradient, EngineeringStrain2D& rDeltaElasticEngineeringStrain) const;

    //! @brief ... calculates the difference of the elastic strain between the current state and the previous update
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rDeltaElasticEngineeringStrain ... delta elastic engineering strain (return value)
    void GetDeltaElasticEngineeringStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient3D& rDeformationGradient, EngineeringStrain3D& rDeltaElasticEngineeringStrain) const;

///////////////////////////////////////////////////////////////////////////
    //! @brief ... get dimension of the constitutive relationship
    //! @return ... dimension of the constitutive relationship (1, 2 or 3)
    int GetGlobalDimension() const;

    //! @brief ... get type of constitutive relationship
    //! @return ... type of constitutive relationship
    //! @sa eConstitutiveType
    Constitutive::eConstitutiveType GetType() const;

    //! @brief ... set the elastic matrix
    //! @param rElasticStiffness... elastic matrix
    NuTo::FullMatrix<double> GetElasticStiffness()const;

    //! @brief ... set the elastic matrix
    //! @param rElasticStiffness... elastic matrix
    void SetElasticStiffness(NuTo::FullMatrix<double> rElasticStiffness);

    //! @brief ... check elastic stiffness
    //! @param rElasticStiffness ... crack transition radius
    void CheckElasticStiffness(const NuTo::FullMatrix<double>& rElasticStiffness) const;

    //! @brief ... return the binary file from which the fine scale model is eventually deserialized
    //! @return name of the file
    std::string GetMultiscaleFile()const;

    //! @brief ... set the binary file from which the fine scale model is eventually deserialized
    //! @param rFileName... name of the file
    void SetMultiscaleFile(std::string rFileName);

    //! @brief ... return crack transition radius to smooth the Heaviside function in the multiscale model
    //! @return crack transition radius
    double GetCrackTransitionRadius()const;

    //! @brief ... crack transition radius to smooth the Heaviside function in the multiscale model
    //! @param rCrackTransitionRadius... crack transition radius
    void SetCrackTransitionRadius(double rCrackTransitionRadius);

    //! @brief ... check if crack transition radius is positive
    //! @param rCrackTransitionRadius ... crack transition radius
    void CheckCrackTransitionRadius(double rCrackTransitionRadius) const;

    //! @brief ... get tensile strength
    //! @return ... tensile strength
    double GetTensileStrength() const;

    //! @brief ... set tensile strength
    //! @param rTensileStrength...  tensile strength
    void SetTensileStrength(double rTensileStrength);

    //! @brief ... check if tensile strength is positive
    //! @param rTensileStrength ... tensile strength
    void CheckTensileStrength(double rTensileStrength) const;

    //! @brief ... get penalty stiffness crack angle
    //! @return ... penalty stiffness crack angle
    double GetPenaltyStiffnessCrackAngle() const;

    //! @brief ... set PenaltyStiffnessCrackAngle
    //! @param rPenaltyStiffnessCrackAngle...  penalty stiffness crack angle
    void SetPenaltyStiffnessCrackAngle(double rPenaltyStiffnessCrackAngle);

    //! @brief ... check if penalty stiffness crack angle is positive
    //! @param rPenaltyStiffnessCrackAngle ... PenaltyStiffnessCrackAngle
    void CheckPenaltyStiffnessCrackAngle(double rPenaltyStiffnessCrackAngle) const;

    //! @brief ... get scaling factor for the dofs of the crack angle
    //! @return ... scaling factor
    double GetScalingFactorCrackAngle() const;

    //! @brief ... set scaling factor for the dofs of the crack angle
    //! @param rScalingFactor...  scaling factor
    void SetScalingFactorCrackAngle(double rScalingFactorCrackAngle);

    //! @brief ... check if ScalingFactorCrackAngle is positive
    //! @param rScalingFactorCrackAngle ... ScalingFactorCrackAngle
    void CheckScalingFactorCrackAngle(double rScalingFactorCrackAngle) const;

    //! @brief ... get scaling factor for the dofs of the crack opening
    //! @return ... scaling factor
    double GetScalingFactorCrackOpening() const;

    //! @brief ... set scaling factor for the dofs of the crack opening
    //! @param rScalingFactor...  scaling factor
    void SetScalingFactorCrackOpening(double rScalingFactorCrackOpening);

    //! @brief ... check if ScalingFactorCrackOpening is positive
    //! @param rScalingFactorCrackOpening ... ScalingFactorCrackOpening
    void CheckScalingFactorCrackOpening(double rScalingFactorCrackOpening) const;

    //! @brief ... get scaling factor for the dofs of total strain
    //! @return ... scaling factor
    double GetScalingFactorEpsilon() const;

    //! @brief ... set scaling factor for the dofs of the total strain
    //! @param rScalingFactor...  scaling factor
    void SetScalingFactorEpsilon(double rScalingFactorEpsilon);

    //! @brief ... check if ScalingFactorEpsilon is positive
    //! @param rScalingFactorEpsilon ... ScalingFactorEpsilon
    void CheckScalingFactorEpsilon(double rScalingFactorEpsilon) const;

    //! @brief ... get AugmentedLagrangeStiffnessCrackOpening
    //! @return ...AugmentedLagrangeStiffnessCrackOpening
    double GetAugmentedLagrangeStiffnessCrackOpening() const;

    //! @brief ... set AugmentedLagrangeStiffnessCrackOpening
    //! @param rAugmentedLagrangeStiffnessCrackOpening...AugmentedLagrangeStiffnessCrackOpening
    void SetAugmentedLagrangeStiffnessCrackOpening(double rAugmentedLagrangeStiffnessCrackOpening);

    //! @brief ... check AugmentedLagrangeStiffnessCrackOpening
    //! @param rAugmentedLagrangeStiffnessCrackOpening ...AugmentedLagrangeStiffnessCrackOpening
    void CheckAugmentedLagrangeStiffnessCrackOpening(double rAugmentedLagrangeStiffnessCrackOpening) const;

    //! @brief ... get ToleranceResidualForce
    //! @return ...ToleranceResidualForce
    double GetToleranceResidualForce() const;

    //! @brief ... set ToleranceResidualForce
    //! @param rToleranceResidualForce... ToleranceResidualForce
    void SetToleranceResidualForce(double rToleranceResidualForce);

    //! @brief ... check ToleranceResidualForce
    //! @param r ...
    void CheckToleranceResidualForce(double rToleranceResidualForce) const;

    //! @brief ... get MaxNumNewtonIterations
    //! @return ...MaxNumNewtonIterations
    int GetMaxNumNewtonIterations() const;

    //! @brief ... MaxNumNewtonIterations
    //! @param rMaxNumNewtonIterations...MaxNumNewtonIterations
    void SetMaxNumNewtonIterations(int rMaxNumNewtonIterations);

    //! @brief ... check MaxNumNewtonIterations
    //! @param rMaxNumNewtonIterations ...MaxNumNewtonIterations
    void CheckMaxNumNewtonIterations(int rMaxNumNewtonIterations) const;

    //! @brief ... get MaxDeltaLoadFactor
    //! @return ...MaxDeltaLoadFactor
    double GetMaxDeltaLoadFactor() const;

    //! @brief ... set MaxDeltaLoadFactor
    //! @param rMaxDeltaLoadFactor...MaxDeltaLoadFactor
    void SetMaxDeltaLoadFactor(double rMaxDeltaLoadFactor);

    //! @brief ... check MaxDeltaLoadFactor
    //! @param rMaxDeltaLoadFactor ...MaxDeltaLoadFactor
    void CheckMaxDeltaLoadFactor(double rMaxDeltaLoadFactor) const;

    //! @brief ... get DecreaseFactor
    //! @return ...DecreaseFactor
    double GetDecreaseFactor() const;

    //! @brief ... set DecreaseFactor
    //! @param rDecreaseFactor...DecreaseFactor
    void SetDecreaseFactor(double rDecreaseFactor);

    //! @brief ... check DecreaseFactor
    //! @param rDecreaseFactor ...DecreaseFactor
    void CheckDecreaseFactor(double r) const;

    //! @brief ... get MinNumNewtonIterations
    //! @return ...MinNumNewtonIterations
    int GetMinNumNewtonIterations() const;

    //! @brief ... set MinNumNewtonIterations
    //! @param rMinNumNewtonIterations...MinNumNewtonIterations
    void SetMinNumNewtonIterations(double rMinNumNewtonIterations);

    //! @brief ... check MinNumNewtonIterations
    //! @param rMinNumNewtonIterations ...
    void CheckMinNumNewtonIterations(int rMinNumNewtonIterations) const;

    //! @brief ... get IncreaseFactor
    //! @return ...IncreaseFactor
    double GetIncreaseFactor() const;

    //! @brief ... set IncreaseFactor
    //! @param rIncreaseFactor...
    void SetIncreaseFactor(double rIncreaseFactor);

    //! @brief ... check IncreaseFactor
    //! @param rIncreaseFactor ...IncreaseFactor
    void CheckIncreaseFactor(double rIncreaseFactor) const;

    //! @brief ... get MinLoadFactor
    //! @return ...MinLoadFactor
    double GetMinLoadFactor() const;

    //! @brief ... set MinLoadFactor
    //! @param rMinLoadFactor...MinLoadFactor
    void SetMinLoadFactor(double rMinLoadFactor);

    //! @brief ... check MinLoadFactor
    //! @param rMinLoadFactor ...MinLoadFactor
    void CheckMinLoadFactor(double rMinLoadFactor) const;

    //! @brief ... get MinLineSearchFactor
    //! @return ... MinLineSearchFactor
    double GetMinLineSearchFactor() const;

    //! @brief ... set MinLineSearchFactor
    //! @param rMinLineSearchFactor...MinLineSearchFactor
    void SetMinLineSearchFactor(double rMinLineSearchFactor);

    //! @brief ... check MinLineSearchFactor
    //! @param rMinLineSearchFactor ...MinLineSearchFactor
    void CheckMinLineSearchFactor(double rMinLineSearchFactor) const;

    //! @brief ... get result directory
    //! @return ... ResultDirectory
    const std::string& GetResultDirectory() const;

    //! @brief ... set ResultDirectory
    //! @param rResultDirectory...ResultDirectory
    void SetResultDirectory(const std::string& rResultDirectory);

    //! @brief ... check ResultDirectory
    //! @param ResultDirectory ...ResultDirectory
    void CheckResultDirectory(const std::string& rResultDirectory) const;

    //! @brief ... get load step macro
    //! @return ... LoadStepMacro
    int GetLoadStepMacro() const;

    //! @brief ... set LoadStepMacro
    //! @param LoadStepMacro...LoadStepMacro
    void SetLoadStepMacro(int rLoadStepMacro);

    //! @brief ... check LoadStepMacro
    //! @param LoadStepMacro ...LoadStepMacro
    void CheckLoadStepMacro(int rLoadStepMacro) const;

    //! @brief ... get load step macro
    //! @return ... LoadStepMacro
    bool GetSquareCoarseScaleModel() const;

    //! @brief ... set LoadStepMacro
    //! @param LoadStepMacro...LoadStepMacro
    void SetSquareCoarseScaleModel(bool rSquareCoarseScaleModel);

    //! @brief ... check LoadStepMacro
    //! @param LoadStepMacro ...LoadStepMacro
    void CheckSquareCoarseScaleModel(bool rSquareCoarseScaleModel) const;

    //! @brief ... check parameters of the constitutive relationship
    void CheckParameters()const;

    //! @brief ... check compatibility between element type and type of constitutive relationship
    //! @param rElementType ... element type
    //! @return ... <B>true</B> if the element is compatible with the constitutive relationship, <B>false</B> otherwise.
    bool CheckElementCompatibility(Element::eElementType rElementType) const;

    //! @brief ... print information about the object
    //! @param rVerboseLevel ... verbosity of the information
    void Info(unsigned short rVerboseLevel) const;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION


    //! @brief ... returns true, if a material model has tmp static data (which has to be updated before stress or stiffness are calculated)
    //! @return ... see brief explanation
    bool HaveTmpStaticData() const;

    //! @brief ... returns true, if a material model has is nonlocal (stiffness is of dynamic size, nonlocal averaging)
    //! @return ... see brief explanation
    bool IsNonlocalModel()const;

protected:
    // calculate coefficients of the linear elastic material matrix
    void CalculateCoefficients3D(double& C11, double& C12, double& C33) const;

    // this is just for debugging purposes
    bool CheckStiffness(NuTo::StructureMultiscale* rFineScaleStructure)const;

    // this is just for debugging purposes
    bool CheckGradient(NuTo::StructureMultiscale* rFineScaleStructure)const;

    //move from elastic to nonlinear solution
    void SwitchToNonlinear(ConstitutiveStaticDataMultiscale2DPlaneStrain *rStaticData, ElementBase* rElement, int rIp, EngineeringStrain2D& rEngineeringStrain)const;

    //! @brief elastic stiffness (before the model has to be transfered to the fine scale)
    NuTo::FullMatrix<double> mElasticStiffness;

    //! @brief stores the file to read in the fine scale model
    std::string mFileName;

    //! @brief stiffness of the augmented Lagrangian to prevent crack opening
    double mAugmentedLagrangeStiffnessCrackOpening;

    //! @brief if true, scale the crack length of the fine scale model based on the crack angle
    bool mSquareCoarseScaleModel;

    //! @brief tensile strength, this parameter should be chosen smaller than the actual value
    //! it is used to determine the transition from the linear elastic model to the full mesoscale model
    double mTensileStrength;

    //! @brief scaling factof for the dofs in order to reflect the different dimensions (remember, the standard displacements are not scaled)
    double mScalingFactorCrackAngle;
    double mScalingFactorCrackOpening;
    double mScalingFactorEpsilon;

    //! @brief scaling factor m for the penalty constraint of the crack angle pot = m/2*(alpha-alpha_e)^2
    double mPenaltyStiffnessCrackAngle;

    //! @brief crack transition radius
    double mCrackTransitionRadius;

    //! @brief parameters for the Newton Raphson iteration on the fine scale
    double mToleranceResidualForce;
    int mMaxNumNewtonIterations;
    double mMaxDeltaLoadFactor;
    double mDecreaseFactor;
    int mMinNumNewtonIterations;
    double mIncreaseFactor;
    double mMinLoadFactor;
    double mMinLineSearchFactor;

    //! @brief criterion to add the crack enrichment function
    double mDamageCrackInitiation;

    //directory, where all the results for the fine scale solutions are stored
    std::string mResultDirectory;
    int mLoadStepMacro;

    //! @brief if true, rotate strain into crack coordinate system, solve and rotate back
    bool mRotateCrack;

};
}

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::Multiscale)
//this is due to the diamond structure (virtual public derivative)
namespace boost{
template<>
struct is_virtual_base_of<NuTo::ConstitutiveEngineeringStressStrain, NuTo::Multiscale>: public mpl::true_ {};
}
#endif // ENABLE_SERIALIZATION

#endif /* CONSTITUTIVEMULTISCALE_H_ */
