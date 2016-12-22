/*
 * GradientDamageEngineeringStressFatigue.h
 *
 *  Created on: 10 Nov 2014
 *      Author: ttitsche
 */

#ifndef CONSTITUTIVEGRADIENTDAMAGEENGINEERINGSTRESSFATIGUE_H_
#define CONSTITUTIVEGRADIENTDAMAGEENGINEERINGSTRESSFATIGUE_H_

#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/elements/ElementEnum.h"

#include "mechanics/constitutive/ConstitutiveBase.h"

namespace NuTo
{
class ConstitutiveStaticDataGradientDamagePlasticity1D;
class ConstitutiveStaticDataGradientDamage1DFatigue;
class ConstitutiveStaticDataGradientDamage2DFatigue;
class Logger;
class ConstitutiveTangentBase;
//! @author Joerg F. Unger
//! @date Apr 26, 2010
//! @brief ...
class GradientDamageEngineeringStressFatigue : public ConstitutiveBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    GradientDamageEngineeringStressFatigue();

    //! @brief ... evaluate the constitutive relation in 1D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
    //! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
    NuTo::Error::eError Evaluate1D(
            ElementBase* rElement, int rIp,
            const std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*>& rConstitutiveInput,
            std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>& rConstitutiveOutput) override;

    //! @brief ... evaluate the constitutive relation in 2D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
    //! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
    NuTo::Error::eError Evaluate2D(
            ElementBase* rElement, int rIp,
            const std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*>& rConstitutiveInput,
            std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>& rConstitutiveOutput) override;

    //! @brief ... evaluate the constitutive relation in 3D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
    //! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
    NuTo::Error::eError Evaluate3D(
            ElementBase* rElement, int rIp,
            const std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*>& rConstitutiveInput,
            std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>& rConstitutiveOutput) override;

    //! @brief ... create new static data object for an integration point
    //! @return ... pointer to static data object
    ConstitutiveStaticDataBase* AllocateStaticDataEngineeringStress_EngineeringStrain1D(const ElementBase* rElement) const override;

    //! @brief ... create new static data object for an integration point
    //! @return ... pointer to static data object
    ConstitutiveStaticDataBase* AllocateStaticDataEngineeringStress_EngineeringStrain2D(const ElementBase* rElement) const override;

    //! @brief ... create new static data object for an integration point
    //! @return ... pointer to static data object
    ConstitutiveStaticDataBase* AllocateStaticDataEngineeringStress_EngineeringStrain3D(const ElementBase* rElement) const override;

    //! @brief ... performs the return mapping procedure in 1D
    //! @param rElement ... structure
    //! @param rIp ... integration point
    //! @param rNonlocalEqStrain ... nonlocal eq strain
    //! @param rEngineeringStrain ... engineering strain
    //! @param rNewStress ... new stress (if a 0-pointer is given, no values are written)
    //! @param rNewTangent ... new tangent matrix (if a 0-pointer is given, no values are written)
    //! @param rNewTangentNonlocal ... new tangent matrix (if a 0-pointer is given, no values are written) derivative of stress by the nonlocal eq strain
    //! @param rNewStaticData ... new static data (if a 0-pointer is given, no values are written)
    NuTo::Error::eError ReturnMapping1D(const ElementBase* rElement,int rIp,
    		const double rNonlocalEqStrain,
    		const EngineeringStrain1D& rEngineeringStrain,
    		EngineeringStress1D* rNewStress,
    		ConstitutiveTangentLocal<1,1>* rNewTangent,
    		ConstitutiveTangentLocal<1,1>* rNewTangentNonlocal,
    		ConstitutiveStaticDataGradientDamage1DFatigue* rNewStaticData,
    		Logger& rLogger) const;

    //! @brief ... performs the return mapping procedure in 2D
    //! @param rElement ... structure
    //! @param rIp ... integration point
    //! @param rNonlocalEqStrain ... nonlocal eq strain
    //! @param rEngineeringStrain ... engineering strain
    //! @param rNewStress2D, rNewStress3D ... new stresses in 2D and 3D (if a 0-pointer is given, no values are written)
    //! @param rNewTangent ... new tangent matrix (if a 0-pointer is given, no values are written)
    //! @param rNewTangentNonlocal ... new tangent matrix (if a 0-pointer is given, no values are written) derivative of stress by the nonlocal eq strain
    //! @param rNewStaticData ... new static data (if a 0-pointer is given, no values are written)
	NuTo::Error::eError ReturnMapping2D(const ElementBase* rElement, int rIp,
			const double rNonlocalEqStrain,
			const EngineeringStrain2D& rEngineeringStrain,
			EngineeringStress2D* rNewStress2D,
			EngineeringStress3D* rNewStress3D,
			ConstitutiveTangentLocal<3,3>* rNewTangent,
			ConstitutiveTangentLocal<3,1>* rNewTangentNonlocal,
			ConstitutiveStaticDataGradientDamage2DFatigue* rNewStaticData,
			Logger& rLogger) const;

    // calculate coefficients of the material matrix
    void CalculateCoefficients3D(double& C11, double& C12, double& C44) const;

    // calculate coefficients of the material matrix
    void CalculateCoefficients2DPlaneStress(double& C11, double& C12, double& C33) const;

    // parameters /////////////////////////////////////////////////////////////

    //! @brief ... gets a parameter of the constitutive law which is selected by an enum
    //! @param rIdentifier ... Enum to identify the requested parameter
    //! @return ... value of the requested variable
    virtual double GetParameterDouble(Constitutive::eConstitutiveParameter rIdentifier) const override;

    //! @brief ... sets a parameter of the constitutive law which is selected by an enum
    //! @param rIdentifier ... Enum to identify the requested parameter
    //! @param rValue ... new value for requested variable
    virtual void SetParameterDouble(Constitutive::eConstitutiveParameter rIdentifier, double rValue) override;


    //! @brief ... gets a parameter of the constitutive law which is selected by an enum
    //! @param rIdentifier ... Enum to identify the requested parameter
    //! @return ... value of the requested variable
    virtual NuTo::FullVector<double,Eigen::Dynamic> GetParameterFullVectorDouble(Constitutive::eConstitutiveParameter rIdentifier) const;

    //! @brief ... sets a parameter of the constitutive law which is selected by an enum
    //! @param rIdentifier ... Enum to identify the requested parameter
    //! @param rValue ... new value for requested variable
    virtual void SetParameterFullVectorDouble(Constitutive::eConstitutiveParameter rIdentifier, NuTo::FullVector<double,Eigen::Dynamic> rValue);

///////////////////////////////////////////////////////////////////////////

    //! @brief ... get type of constitutive relationship
    //! @return ... type of constitutive relationship
    //! @sa eConstitutiveType
    Constitutive::eConstitutiveType GetType() const override;

    //! @brief ... check parameters of the constitutive relationship
    void CheckParameters()const override;

    //! @brief ... check compatibility between element type and type of constitutive relationship
    //! @param rElementType ... element type
    //! @return ... <B>true</B> if the element is compatible with the constitutive relationship, <B>false</B> otherwise.
    bool CheckElementCompatibility(Element::eElementType rElementType) const override;

    //! @brief ... print information about the object
    //! @param rVerboseLevel ... verbosity of the information
    void Info(unsigned short rVerboseLevel, Logger& rLogger) const override;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

    //! @brief ... returns true, if a material model has tmp static data (which has to be updated before stress or stiffness are calculated)
    //! @return ... see brief explanation
    bool HaveTmpStaticData() const override;

    //! @brief ... calculate the isotropic damage variable from the nonlocal eq strain history variable
    //! kappa and the damage law parameters (class members)
    //! \f$ \omega = 1 - \frac{\varepsilon_0}{\kappa} \exp \left(\frac{\varepsilon_0 - \kappa}{\varepsilon_f} \right) \f$
    //! @param rKappa ... history variable
    //! @return isotropic damage variable
    double CalculateDamage(double rKappa) const;

    //! @brief ... calculate the D_Omega_D_Kappa from the nonlocal eq strain history variable
    //! kappa and the damage law parameters (class members)
    //! \f$ \frac{\partial \omega}{\partial \kappa} = \frac{\varepsilon_0}{\kappa} \left(\frac{1}{\kappa} + \frac{1}{\varepsilon_f} \right) \exp \left(\frac{\varepsilon_0 - \kappa}{\varepsilon_f} \right) \f$
    //! @param rKappa ... history variable
    //! @return ... isotropic damage variable
    double CalculateDerivativeDamage(double rKappa) const;

    //! @brief ... calculate the second derivative DD_Omega_DD_Kappa from the nonlocal eq strain history variable
    //! kappa and the damage law parameters (class members)
    //! \f$ \frac{\partial \omega}{\partial \kappa} = \frac{\varepsilon_0}{\kappa} \left(\frac{1}{\kappa} + \frac{1}{\varepsilon_f} \right) \exp \left(\frac{\varepsilon_0 - \kappa}{\varepsilon_f} \right) \f$
    //! @param rKappa ... history variable
    //! @return ... isotropic damage variable
    double CalculateSecondDerivativeDamage(double rKappa) const;

    //! @brief returns the local eq strain calculated from 2d strain
    //! @param rStrain2D ... 2d strain
    NuTo::LocalEqStrain CalculateLocalEqStrain2D(const NuTo::EngineeringStrain2D& rStrain2D) const;

    //! @brief returns the local eq strain calculated from 3d strain
    //! @param rStrain3D ... 3d strain
    NuTo::LocalEqStrain CalculateLocalEqStrain3D(const NuTo::EngineeringStrain3D& rStrain3D) const;

    //! @brief returns the derivative of the local eq strain with respect to the 2d strain
    //! @param rStrain2D ... 2d strain
    NuTo::ConstitutiveTangentLocal<3,1> CalculateLocalEqStrainTangent2D(const NuTo::EngineeringStrain2D& rStrain2D) const;

    //! @brief returns the derivative of the local eq strain with respect to the 3d strain
    //! @param rStrain3D ... 3d strain
    NuTo::ConstitutiveTangentLocal<6,1> CalculateLocalEqStrainTangent3D(const NuTo::EngineeringStrain3D& rStrain3D) const;

    //! @brief calculates nonlocal eq strain and its derivatives for the modified mises model
    //! @param rStrain2D ... 2d strain
    //! @param rLocalEqStrain ... local eq strain
    //! @param rLocalEqStrainTangent ... d local eq strain / d strain 2D
    void CalculateLocalEqStrainAndDerivativeModifiedMises2DPlaneStrain(const NuTo::EngineeringStrain2D& rStrain2D, NuTo::LocalEqStrain& rLocalEqStrain, NuTo::ConstitutiveTangentLocal<3,1>& rLocalEqStrainTangent) const;

    //! @brief calculates the local eq strain and its derivatives for the modified mises model
    //! @param rStrain3D ... 3d strain
    //! @param rLocalEqStrain ... local eq strain
    //! @param rLocalEqStrainTangent ... d local eq strain / d strain 3D
    void CalculateLocalEqStrainAndDerivativeModifiedMises3D(const NuTo::EngineeringStrain3D& rStrain3D, NuTo::LocalEqStrain& rLocalEqStrain, NuTo::ConstitutiveTangentLocal<6, 1>& rLocalEqStrainTangent) const;

    //! @brief calculates the local eq strain and its derivatives for the modified mises model in plane stress state
    //! @param rStrain2D ... 2d strain
    //! @param rLocalEqStrain ... local eq strain
    //! @param rLocalEqStrainTangent ... d local eq strain / d strain 2D
    void CalculateLocalEqStrainAndDerivativeModifiedMises2DPlaneStress(const NuTo::EngineeringStrain2D& rStrain2D, NuTo::LocalEqStrain& rLocalEqStrain, NuTo::ConstitutiveTangentLocal<3, 1>& rLocalEqStrainTangent) const;


protected:
    //! @brief ... density
    double mRho;

    //! @brief ... Young's modulus \f$ E \f$
    double mE;

    //! @brief ... Poisson's ratio \f$ \nu \f$
    double mNu;

    //! @brief ... nonlocal radius
    double mNonlocalRadius;

    //! @brief ... nonlocal radius parameter for non constant c-paramter
    double mNonlocalRadiusParameter;

    //! @brief ... thermal expansion coefficient \f$ \alpha \f$
    double mThermalExpansionCoefficient;

    //! @brief ... tensile strength
    double mTensileStrength;

    //! @brief ... uniaxial compressive strength
    double mCompressiveStrength;

    //! @brief ... fracture energy
    double mFractureEnergy;

    //! @brief ... exponent for the damage evolution law
    double mViscosityExponent;

    //! @brief ... damage law type
    int mDamageLawType;

    //! @brief ... damage law parameters
    FullVector<double, Eigen::Dynamic> mDamageLawParameters;

    //! @brief ... check if density is non negative
    //! @param rRho ... Young's modulus
    void CheckDensity(double rRho) const;

    //! @brief ... check if Young's modulus is positive
    //! @param rE ... Young's modulus
    void CheckYoungsModulus(double rE) const;

    //! @brief ... check if Poisson's ratio is valid \f$ (-1.0 < \nu < 0.5) \f$
    //! @param rNu ... Poisson's ratio
    void CheckPoissonsRatio(double rNu) const;


    //! @brief ... check if nonlocal radius is positive
    //! @param rRadius ... nonlocal radius
    void CheckNonlocalRadius(double rNonlocalRadius) const;

    //! @brief ... check if the nonlocal radius parameter is greater than 1
    //! @param rRadiusParameter ... nonlocal radius parameter
    void CheckNonlocalRadiusParameter(double rNonlocalRadiusParameter) const;

    //! @brief ... check thermal expansion coefficient
    //! @param rAlpha ... thermal expansion coefficient
    void CheckThermalExpansionCoefficient(double rAlpha) const;

    //! @brief ... check if tensile strength is positive
    //! @param rTensileStrength ... nonlocal radius
    void CheckTensileStrength(double rTensileStrength) const;

    //! @brief ... check if compressive strength is positive
    //! @param rRadius ... compressive strength
    void CheckCompressiveStrength(double rCompressiveStrength) const;

    //! @brief ... check if fracture energy is positive
    //! @param rFractureEnergy ... fracture energy
    void CheckFractureEnergy(double rFractureEnergy) const;

    //! @brief ... check if viscosity exponent is equal or larger then 1
    //! @param rViscosityExponent ... viscosity exponent
    void CheckViscosityExponent(double rViscosityExponent) const;

    //! @brief ... check damage law
    //! @param rDamageLaw ... damage law
    void CheckDamageLaw(const NuTo::FullVector<double, Eigen::Dynamic>& rDamageLaw) const;

private:
#define Plus(rVP) ((rVP >= 0.)?rVP:0.)				// define McCauley brackets
#define Sgn(rVP) ((rVP >= 0.)?1.:0.)				// define positive sgn function
#define DELTA(i, j) ((i==j)?1:0)					// define Kronecker Symbol
//#define mTOLF 1.0e-12								// define tolerance for Newton

    //! @brief ... calculates the residual vector of an incremental formulation
    NuTo::FullVector<double,Eigen::Dynamic> Residual(
    		const NuTo::FullVector<double,Eigen::Dynamic> &rParameter,
    		NuTo::FullVector<double,Eigen::Dynamic> rUnknown) const
		{
    		if (rParameter.size()!=4) {
    			throw MechanicsException("[NuTo::GradientDamageEngineeringStressFatigue::Residual] the parameter vector expects 4 components.");
    		}

			NuTo::FullVector<double,Eigen::Dynamic> residual(rUnknown.GetNumRows());

			double rPrevOmega, rPrevKappa, rPrevNonlocalEqStrain, rNonlocalEqStrain;
			double rDeltaOmega, rDeltaKappa;

			// initialization of physical values from the vectors rParameter and rUnknown
			rPrevOmega = rParameter[0];					// rParameter(0) omega at the beginning of the time increment
			rPrevKappa = rParameter[1];					// rParameter(1) kappa at the beginning of the time increment
			rPrevNonlocalEqStrain = rParameter[2];		// rParameter(2) nonlocal eq strain at the beginning of the time increment
			rNonlocalEqStrain     = rParameter[3];		// current nonlocal eq strain

			rDeltaOmega = rUnknown[0];					// delta omega
			rDeltaKappa = rUnknown[1];					// delta kappa

			// define current omega and kappa
			double omega(rPrevOmega + rDeltaOmega);
			double kappa(rPrevKappa + rDeltaKappa);

			// calculate the residual vector
			// damage equation
			residual[0] = rDeltaOmega - CalculateDerivativeDamage(kappa) * pow(rNonlocalEqStrain/kappa,mViscosityExponent) * Plus(rNonlocalEqStrain - rPrevNonlocalEqStrain);

			// kappa equation
			residual[1] = omega - CalculateDamage(kappa);

			return residual;
		}

    //! @brief ... calculates the analytical matrix:= derivative of the Residual by nonlocal eq strain
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> DResidualDEpsNonlocalAn(
    		const NuTo::FullVector<double,Eigen::Dynamic> &rParameter,
    		NuTo::FullVector<double,Eigen::Dynamic> rUnknown) const
		{
			if (rParameter.size()!=4) {
				throw MechanicsException("[NuTo::GradientDamageEngineeringStressFatigue::DResidualDEpsNonlocalAn] the parameter vector expects 4 components.");
			}
			NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> deriv(rUnknown.GetNumRows(),1);

			double rPrevKappa, rPrevNonlocalEqStrain, rNonlocalEqStrain;
			double rDeltaKappa;

			// initialization of physical values from the vectors rParameter and rUnknown
			rPrevKappa = rParameter[1];					// rParameter(1) kappa at the beginning of the time increment
			rPrevNonlocalEqStrain = rParameter[2];		// rParameter(2) nonlocal eq strain at the beginning of the time increment
			rNonlocalEqStrain     = rParameter[3];		// current nonlocal eq strain

			rDeltaKappa = rUnknown[1];					// delta kappa

			// define kappa
			double kappa(rPrevKappa + rDeltaKappa);

			// compute the components
			deriv(0,0)  = mViscosityExponent * pow(rNonlocalEqStrain/kappa,mViscosityExponent - 1.) * Plus(rNonlocalEqStrain - rPrevNonlocalEqStrain) / kappa;
			deriv(0,0) += pow(rNonlocalEqStrain/kappa,mViscosityExponent) * Sgn(rNonlocalEqStrain - rPrevNonlocalEqStrain);
			deriv(0,0) *= -CalculateDerivativeDamage(kappa);

			deriv(1,0) = 0.;

			return deriv;
		}

    //! @brief ... calculates the analytical Jacobi matrix:= derivative of the Residual
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> DResidualAn(
    		const NuTo::FullVector<double,Eigen::Dynamic> &rParameter,
    		NuTo::FullVector<double,Eigen::Dynamic> rUnknown) const
		{
			if (rParameter.size()!=4) {
				throw MechanicsException("[NuTo::GradientDamageEngineeringStressFatigue::DResidualAn] the parameter vector expects 4 components.");
			}

    		NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> deriv(rUnknown.GetNumRows(),rUnknown.GetNumRows());

			double rPrevKappa, rPrevNonlocalEqStrain, rNonlocalEqStrain;
			double rDeltaKappa;

			// initialization of physical values from the vectors rParameter and rUnknown
			rPrevKappa = rParameter[1];					// rParameter(1) kappa at the beginning of the time increment
			rPrevNonlocalEqStrain = rParameter[2];		// rParameter(2) nonlocal eq strain at the beginning of the time increment
			rNonlocalEqStrain     = rParameter[3];		// current nonlocal eq strain

			rDeltaKappa = rUnknown[1];					// delta kappa

			// define current kappa
			double kappa(rPrevKappa + rDeltaKappa);

			deriv.Zero(2,2);

			// compute the components
			deriv(0,0) = 1.;

			deriv(0,1)  = -CalculateSecondDerivativeDamage(kappa) + mViscosityExponent * CalculateDerivativeDamage(kappa) / kappa;
			deriv(0,1) *= pow(rNonlocalEqStrain/kappa,mViscosityExponent) * Plus(rNonlocalEqStrain - rPrevNonlocalEqStrain);

			deriv(1,0) = 1.;

			deriv(1,1) = -CalculateDerivativeDamage(kappa);

			return deriv;
		}

};
}

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::GradientDamageEngineeringStressFatigue)
#endif // ENABLE_SERIALIZATION

#endif /* CONSTITUTIVEGRADIENTDAMAGEENGINEERINGSTRESSFATIGUE_H_ */
