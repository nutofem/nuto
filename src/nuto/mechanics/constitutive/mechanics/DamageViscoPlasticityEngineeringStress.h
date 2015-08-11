// $Id$

#ifndef DAMAGEVISCOPLASTICITYENGINEERINGSTRESS_H_
#define DAMAGEVISCOPLASTICITYENGINEERINGSTRESS_H_

#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/elements/ElementEnum.h"
#include "nuto/mechanics/constitutive/ConstitutiveBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress3D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain3D.h"

//#include "nuto/optimize/NewtonRaphson.h"

using namespace std;

namespace NuTo
{
//! @brief ... visco-plastic material model with damage
class ConstitutiveStaticDataDamageViscoPlasticity3D;
class Logger;
template <int TNumRows,int TNumColumns>
class ConstitutiveTangentLocal;

class DamageViscoPlasticityEngineeringStress: public ConstitutiveBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    DamageViscoPlasticityEngineeringStress();

    //! @brief ... evaluate the constitutive relation in 2D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rUpdateHistory ... update history variables after leaving the routine
    //! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
    //! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
    NuTo::Error::eError Evaluate1D(ElementBase* rElement, int rIp,
    		const std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*>& rConstitutiveInput,
    		std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>& rConstitutiveOutput) override;

    //! @brief ... evaluate the constitutive relation in 2D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rUpdateHistory ... update history variables after leaving the routine
    //! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
    //! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
    NuTo::Error::eError Evaluate2D(ElementBase* rElement, int rIp,
    		const std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*>& rConstitutiveInput,
    		std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>& rConstitutiveOutput) override;

    //! @brief ... evaluate the constitutive relation in 3D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rUpdateHistory ... update history variables after leaving the routine
    //! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
    //! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
    NuTo::Error::eError Evaluate3D(ElementBase* rElement, int rIp,
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

    //! @brief ... performs the return mapping procedure in 3D
    //! @param rElement ... structure
    //! @param rIp ... integration point
    //! @param rEngineeringStrain ... engineering strain
    //! @param rNewStress ... new stress (if a 0-pointer is given, no values are written)
    //! @param rNewTangent ... new tangent matrix (if a 0-pointer is given, no values are written)
    //! @param rNewStaticData ... new static data (if a 0-pointer is given, no values are written)
    NuTo::Error::eError ReturnMapping3D(const ElementBase* rElement,int rIp,
    		const EngineeringStrain3D& rEngineeringStrain,
    		EngineeringStress3D* rNewStress,
    		ConstitutiveTangentLocal<6,6>* rNewTangent,
    		ConstitutiveStaticDataDamageViscoPlasticity3D* rNewStaticData,
    		Logger& rLogger) const;

    ///////////////////////////////////////////////////////////////////////////


    // calculate coefficients of the material matrix
    void CalculateCoefficients2DPlainStress(double& C11, double& C12, double& C33) const;
    void CalculateCoefficients3D(double& C11, double& C12, double& C44) const;

    //! @brief ... calculate the ductility parameter
    //! @brief ... return ductility
    double CalculateDuctility()const;

    // parameters /////////////////////////////////////////////////////////////

    //! @brief ... gets a variable of the constitutive law which is selected by an enum
    //! @param rIdentifier ... Enum to identify the requested variable
    //! @return ... value of the requested variable
    virtual double GetParameterDouble(Constitutive::eConstitutiveParameter rIdentifier) const override;

    //! @brief ... sets a variable of the constitutive law which is selected by an enum
    //! @param rIdentifier ... Enum to identify the requested variable
    //! @param rValue ... new value for requested variable
    virtual void SetParameterDouble(Constitutive::eConstitutiveParameter rIdentifier, double rValue) override;


    //! @brief ... get fracture energy
    //! @return ... fracture energy
    double GetFractureEnergy() const;

    //! @brief ... set fracture energy
    //! @param rFractureEnergy... fracture energy
    void SetFractureEnergy(double rFractureEnergy);

    ///////////////////////////////////////////////////////////////////////////

    //! @brief ... get type of constitutive relationship
    //! @return ... type of constitutive relationship
    //! @sa eConstitutiveType
    Constitutive::eConstitutiveType GetType() const override;

    //! @brief ... check parameters of the constitutive relationship
    void CheckParameters()const;

    //! @brief ... check compatibility between element type and type of constitutive relationship
    //! @param rElementType ... element type
    //! @return ... <B>true</B> if the element is compatible with the constitutive relationship, <B>false</B> otherwise.
    bool CheckElementCompatibility(Element::eElementType rElementType) const override;

    //! @brief ... print information about the object
    //! @param rVerboseLevel ... verbosity of the information
    //! @param rLogger stream for the output
    void Info(unsigned short rVerboseLevel, Logger& rLogger) const override;

    //! @brief ... returns true, if a material model has tmp static data (which has to be updated before stress or stiffness are calculated)
    //! @return ... see brief explanation
    bool HaveTmpStaticData() const override
    {
    	return false;
    }

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

protected:
    //! @brief ... Young's modulus \f$ E \f$
    double mE;

    //! @brief ... Poisson's ratio \f$ \nu \f$
    double mNu;

    //! @brief ... density \f$ \rho \f$
    double mRho;

    //! @brief ... thermal expansion coefficient \f$ \alpha \f$
    double mThermalExpansionCoefficient;

    //! @param ... Tensile Strength
    double mTensileStrength;

    //! @param ... Compressive Strength
    double mCompressiveStrength;

    //! @param ... Biaxial Compressive Strength
    double mBiaxialCompressiveStrength;

    //! @param ... Viscosity
    double mViscosity;

    //! @param ... Damage Distribution (determines the portion of damage via viscoplasticity and plasticity)
    double mDamageDistribution;

    //! @param ... Viscoplastic Yield Surface Offset with respect to the plastic yield surface
    double mViscoplasticYieldSurfaceOffset;

    //! @param ... Fracture energy
    double mFractureEnergy;

    //! @brief ... either use a pure plasticity model (false) or add softening using the damage model (true)
    bool mDamage;

    //! @brief ... check if density is positive
    //! @param rRho ... density
    void CheckDensity(double rRho) const;

    //! @brief ... check if Young's modulus is positive
    //! @param rE ... Young's modulus
    void CheckYoungsModulus(double rE) const;

    //! @brief ... check if Poisson's ratio is valid
    //! @param rNu ... Poisson's ratio
    void CheckPoissonsRatio(double rNu) const;

    //! @brief ... check thermal expansion coefficient
    //! @param rAlpha ... thermal expansion coefficient
    void CheckThermalExpansionCoefficient(double rAlpha) const;

    //! @brief ... check if tensile strength is positive
    //! @param rTensileStrength ...
    void CheckTensileStrength(double rTensileStrength) const;

    //! @brief ... check if compressive strength is positive
    //! @param rRadius ... compressive strength
    void CheckCompressiveStrength(double rCompressiveStrength) const;

    //! @brief ... check if biaxial compressive strength is positive
    //! @param rBiaxialCompressiveStrength ... biaxial compressive strength
    void CheckBiaxialCompressiveStrength(double rBiaxialCompressiveStrength) const;

    //! @brief ... check if viscosity is positive
    //! @param rViscosity ... viscosity
    void CheckViscosity(double rViscosity) const;

    //! @brief ... check whether the damage distribution ranges between 0 and 1
    //! @param rDamageDistribution ... damage distribution
    void CheckDamageDistribution(double rDamageDistribution) const;

    //! @brief ... check whether the viscoplastic yield surface offset is negative
    //! @param rViscoplasticYieldSurfaceOffset ... viscoplastic yield surface offset
    void CheckViscoplasticYieldSurfaceOffset(double rViscoplasticYieldSurfaceOffset) const;

    //! @brief ... check whether fracture energy is positive
    //! @param rFractureEnergy ... fracture energy
    void CheckFractureEnergy(double rFractureEnergy) const;

private:

#define Plus(rVP) ((rVP >= 0.)?rVP:0.)				// define McCauley brackets
#define Sgn(rVP) ((rVP >= 0.)?1.:0.)				// define positive sgn function
#define DELTA(i, j) ((i==j)?1:0)					// define Kronecker Symbol
#define mTOLF 1.0e-12*this->mE						// define tolerance for Newton


    //! @brief ... calculates the residual vector of an incremental formulation
    //! @brief ... rElasticEngineeringStrain ... elastic engineering strain at the end of time increment
    NuTo::FullVector<double,Eigen::Dynamic> Residual(
    		const NuTo::FullVector<double,Eigen::Dynamic> &rParameter,
    		NuTo::FullVector<double,Eigen::Dynamic> rUnknown) const
		{
			NuTo::FullVector<double,Eigen::Dynamic> residual(rUnknown.GetNumRows());

			EngineeringStress3D rPrevStress;      // stress at the beginning of the time increment
			EngineeringStrain3D rDeltaStrain;     // mechanical strain increment (strain increment without thermal component)
			double rDeltaTime, rBeta, rHP, rHVP, rViscosity;  // time increment, rBeta, rHP, rHVP, rViscosity

			EngineeringStress3D rDeltaStress, rDeltaProjStress;		// stress increment, projection stress increment
			double rVP, rVviscoP;					// state variable, its positive counterpart is cumulative plastic and viscoplastic
													// strain rate respectively

			// initialization of physical values from the vectors rParameter and rUnknown
			for (int i = 0; i < 6; i++) {
				rDeltaStrain[i] = rParameter[i];	// rParameter(0:5) increment of mechanical strain
				rPrevStress[i] = rParameter[i+6];	// rParameter(6:11) stress at the beginning of the time increment
				rDeltaStress[i] = rUnknown[i];		// rUnknown(0:5) stress increment
				rDeltaProjStress[i] = rUnknown[i+6];// rUnknown(6:11) projection stress increment
			}

			rDeltaTime = rParameter[12]; 			// rParameter(12) time increment itself
			rBeta = rParameter[13];					// rParameter(13) rBeta
			rHP = rParameter[14];					// rParameter(14) rHP
			rHVP = rParameter[15];					// rParameter(15) rHVP
			rViscosity = rParameter[16];			// rParameter(16) rViscosity

			rVP = rUnknown[12];						// rUnknown(12) is the state variable associated with cumulative plastic flow rate
			rVviscoP = rUnknown[13];				// rUnknown(13) is the state variable associated with cumulative viscoplastic flow rate

			// elastic stiffness
			NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> ElasticStiffness(6, 6);
			double C11, C12, C44;
			this->CalculateCoefficients3D(C11, C12, C44);
			ElasticStiffness << C11, C12, C12,  0., 0.,  0.,
								C12, C11, C12,  0., 0.,  0.,
								C12, C12, C11,  0., 0.,  0.,
								 0.,  0.,  0., C44, 0.,  0.,
								 0.,  0.,  0.,  0., C44, 0.,
								 0.,  0.,  0.,  0., 0., C44;

			NuTo::FullMatrix<double,6,6> d2F_d2Sigma;	// The second derivative of the flow surface by the stress
			NuTo::FullVector<double,6> dF_dSigma;		// The firtst derivatife of the flow surface by the stress
			EngineeringStress3D Stress, ProjStress;	// Stress and projection stress at the end of the time increment

			// calculate Stress and derivatives of the flow surface by the Stress
			Stress = rPrevStress + rDeltaStress;
			ProjStress = rPrevStress + rDeltaProjStress;
			Stress.YieldSurfaceDruckerPrager3DDerivatives(dF_dSigma,d2F_d2Sigma,rBeta);

			// calculate residual with respect to the stress increments, that is residual(0:5)
			residual.segment<6>(0) = rDeltaStress + rDeltaTime*Sgn(Stress.YieldSurfaceDruckerPrager3D(rBeta, rHVP))*
					(rDeltaStress - rDeltaProjStress)/rViscosity - ElasticStiffness*
					(rDeltaStrain - rDeltaTime*Plus(rVP)*dF_dSigma);

			// calculate residual with respect to the projection stress increments, that is residual(6:11)
			ProjStress.YieldSurfaceDruckerPrager3DDerivatives(dF_dSigma,d2F_d2Sigma,rBeta);
			residual.segment<6>(6) = rDeltaProjStress - ElasticStiffness*
					(rDeltaStrain - rDeltaTime*Plus(rVviscoP)*dF_dSigma);

			// calculate residual with respect to the cumulative plastic strain rate, that is residual(12)
			residual[12] = this->mE*rDeltaTime*(rVP - Plus(rVP)) -
					Stress.YieldSurfaceDruckerPrager3D(rBeta, rHP);

			// calculate residual with respect to the cumulative viscoplastic strain rate, that is residual(13)
			residual[13] = this->mE*rDeltaTime*(rVviscoP - Plus(rVviscoP)) -
					ProjStress.YieldSurfaceDruckerPrager3D(rBeta, rHVP);


//			residual[0] = 2.*rUnknown[1]*rUnknown[2] + 2.*rUnknown[0]*rUnknown[1];
//			residual[1] = rUnknown[1]*rUnknown[1] + 2.*rUnknown[0]*rUnknown[2];
//			residual[2] = rUnknown[1]*rUnknown[1] + rUnknown[2]*rUnknown[2] -3.;

			return residual;
		}

    //! @brief ... calculates the analytical matrix:= derivative of the Residual by dEngineeringStress3D
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> DResidualDEpsAn(NuTo::FullVector<double,Eigen::Dynamic> rUnknown) const
		{
    	NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> deriv(rUnknown.GetNumRows(),6);
    	NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> ElasticStiffness(6, 6);
    	double C11, C12, C44;
    	this->CalculateCoefficients3D(C11, C12, C44);
    	ElasticStiffness << C11, C12, C12,  0., 0.,  0.,
    						C12, C11, C12,  0., 0.,  0.,
    						C12, C12, C11,  0., 0.,  0.,
    						 0.,  0.,  0., C44, 0.,  0.,
    						 0.,  0.,  0.,  0., C44, 0.,
    						 0.,  0.,  0.,  0., 0., C44;

    	for (int i = 0; i < 6; ++i) {
			for (int j = 0; j < 6; ++j) {
				deriv(i,j) = -ElasticStiffness(i,j);
				deriv(i+6,j) = -ElasticStiffness(i,j);
			}
		}

		for (int j = 0; j < 6; ++j) {
			deriv(12,j) = 0.;
			deriv(13,j) = 0.;
		}

        return deriv;
		}

    //! @brief ... calculates the analytical Jacobi matrix:= derivative of the Residual
    //! @brief ... rElasticEngineeringStrain ... elastic engineering strain at the end of time increment
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> DResidualAn(
    		const NuTo::FullVector<double,Eigen::Dynamic> &rParameter,
    		NuTo::FullVector<double,Eigen::Dynamic> rUnknown) const
		{
    		NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> deriv(rUnknown.GetNumRows(),rUnknown.GetNumRows());

    		EngineeringStress3D rPrevStress;      // stress at the beginning of the time increment
    		EngineeringStrain3D rDeltaStrain;     // mechanical strain increment (strain increment without thermal component)
            double rDeltaTime, rBeta, rHVP, rViscosity;  // time increment, rBeta, rHP, rHVP, rViscosity

    		EngineeringStress3D rDeltaStress, rDeltaProjStress;		// stress increment, projection stress increment
    		double rVP, rVviscoP;					// state variable, its positive counterpart is cumulative plastic and viscoplastic
												// strain rate respectively

    		// initialization of physical values from the vectors rParameter and rUnknown
    		for (int i = 0; i < 6; i++) {
    			rDeltaStrain[i] = rParameter[i];	// rParameter(0:5) increment of mechanical strain
    			rPrevStress[i] = rParameter[i+6];	// rParameter(6:11) stress at the beginning of the time increment
    			rDeltaStress[i] = rUnknown[i];		// rUnknown(0:5) stress increment
    			rDeltaProjStress[i] = rUnknown[i+6];// rUnknown(6:11) projection stress increment
    		}

    		rDeltaTime = rParameter[12]; 			// rParameter(12) time increment itself
    		rBeta = rParameter[13];					// rParameter(13) rBeta
            //rHP = rParameter[14];					// rParameter(14) rHP
    		rHVP = rParameter[15];					// rParameter(15) rHVP
    		rViscosity = rParameter[16];			// rParameter(16) rViscosity

    		rVP = rUnknown[12];						// rUnknown(12) is the state variable associated with cumulative plastic flow rate
    		rVviscoP = rUnknown[13];				// rUnknown(13) is the state variable associated with cumulative viscoplastic flow rate

    		// elastic stiffness
    		NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> ElasticStiffness(6, 6);
    		double C11, C12, C44;
    		this->CalculateCoefficients3D(C11, C12, C44);
    		ElasticStiffness << C11, C12, C12,  0., 0.,  0.,
								C12, C11, C12,  0., 0.,  0.,
								C12, C12, C11,  0., 0.,  0.,
								0.,  0.,  0., C44, 0.,  0.,
								0.,  0.,  0.,  0., C44, 0.,
								0.,  0.,  0.,  0., 0., C44;

    		NuTo::FullMatrix<double,6,6> d2F_d2Sigma;	 	// The second derivative of the flow surface by the stress
    		NuTo::FullVector<double,6> dF_dSigma, Temp;	// The firtst derivatife of the flow surface by the stress
    		EngineeringStress3D Stress, ProjStress;		// Stress and projection stress at the end of the time increment

    		// calculate Stress and derivatives of the flow surface by the Stress
    		Stress = rPrevStress + rDeltaStress;
    		ProjStress = rPrevStress + rDeltaProjStress;
    		Stress.YieldSurfaceDruckerPrager3DDerivatives(dF_dSigma,d2F_d2Sigma,rBeta);

			deriv.Zero(14,14);
    		Temp.Zero();

			// components DResidualAn(0:5,0:5)
			for (int j = 0; j < 6; j++) {
		    	Temp = ElasticStiffness*d2F_d2Sigma.col(j);
		    	for (int i = 0; i < 6; i++) {
		    		deriv(i,j) = DELTA(i, j) *(1. + rDeltaTime*Sgn(Stress.YieldSurfaceDruckerPrager3D(rBeta, rHVP))/rViscosity) +
		    				rDeltaTime*Plus(rVP)*Temp[i];
				}
			}
			// components DResidualAn(0:5,6:11)
			for (int j = 6; j < 12; j++) {
		    	for (int i = 0; i < 6; i++) {
		    		deriv(i,j) = - DELTA(i, j-6) * rDeltaTime*Sgn(Stress.YieldSurfaceDruckerPrager3D(rBeta, rHVP))/rViscosity;
				}
			}
			// components DResidualAn(0:5,12)
			Temp = ElasticStiffness*dF_dSigma;
			for (int i = 0; i < 6; i++) {
				deriv(i,12) = rDeltaTime*Sgn(rVP)*Temp[i];
			}
			// components DResidualAn(12,0:5)
			for (int j = 0; j < 6; j++) {
				deriv(12,j) = -dF_dSigma[j];
			}
			// components DResidualAn(0:5,13) are already zeroed
			// components DResidualAn(6:11,0:5) are already zeroed
			// components DResidualAn(6:11,6:11)
			ProjStress.YieldSurfaceDruckerPrager3DDerivatives(dF_dSigma,d2F_d2Sigma,rBeta);
			for (int j = 6; j < 12; j++) {
		    	Temp = ElasticStiffness*d2F_d2Sigma.col(j-6);
		    	for (int i = 6; i < 12; i++) {
		    		deriv(i,j) = DELTA(i, j) + rDeltaTime*Plus(rVviscoP)*Temp[i-6];
				}
			}
			// components DResidualAn(6:11,12) are already zeroed
			// components DResidualAn(12,6:11) are already zeroed
			// components DResidualAn(6:11,13)
			Temp = ElasticStiffness*dF_dSigma;
			for (int i = 6; i < 12; i++) {
				deriv(i,13) = rDeltaTime*Sgn(rVviscoP)*Temp[i-6];
			}
			// components DResidualAn(12,12)
			deriv(12,12) = this->mE*rDeltaTime*(1. - Sgn(rVP));
			// components DResidualAn(12,13) is already zeroed
			// components DResidualAn(13,12) is already zeroed
			// components DResidualAn(13,0:5) are already zeroed
			// components DResidualAn(13,6:11)
			for (int j = 6; j < 12; j++) {
				deriv(13,j) = -dF_dSigma[j-6];
			}
			// components DResidualAn(13,13)
			deriv(13,13) = this->mE*rDeltaTime*(1. - Sgn(rVviscoP));

//			deriv(0,0) = 2*rUnknown[1];
//			deriv(0,1) = 2*rUnknown[0] + 2*rUnknown[2];
//			deriv(0,2) = 2*rUnknown[1];
//
//			deriv(1,0) = 2*rUnknown[2];
//			deriv(1,1) = 2*rUnknown[1];
//			deriv(1,2) = 2*rUnknown[0];
//
//			deriv(2,0) = 0.;
//			deriv(2,1) = 2*rUnknown[1];
//			deriv(2,2) = 2*rUnknown[2];

			return deriv;
		}

    //! @brief ... calculates 0.5*Residual^2 and updates fvec = Residual(rz)
    double Fmin(const NuTo::FullVector<double,Eigen::Dynamic> &rParameter,
    		NuTo::FullVector<double,Eigen::Dynamic> rz, NuTo::FullVector<double,Eigen::Dynamic> &fvec) const {
    	fvec = this->Residual(rParameter,rz);
    	double sum=0;
		sum = fvec.dot(fvec);
		return 0.5*sum;
    }

    //! @brief ... calculates the numerical Jacobi matrix:= derivative of the Residual
    //! @brief ... rElasticEngineeringStrain ... elastic engineering strain at the end of time increment
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> DResidualNum(const NuTo::FullVector<double,Eigen::Dynamic> &rParameter,
    		NuTo::FullVector<double,Eigen::Dynamic> rUnknown,
    		NuTo::FullVector<double,Eigen::Dynamic> &fvec) const {
    	const double EPS = 1.0e-8;
    	int n=rUnknown.GetNumRows();
    	NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> deriv(n,n);
    	NuTo::FullVector<double,Eigen::Dynamic> xh=rUnknown;
    	for (int j=0;j<n;j++) {
    		double temp=xh[j];
    		double h=EPS*fabs(temp);
    	   	if (h == 0.0) h=EPS;
    	   		xh[j]=temp+h;
    	   		h=xh[j]-temp;
    	   		NuTo::FullVector<double,Eigen::Dynamic> f=this->Residual(rParameter,xh);
    	   		xh[j]=temp;
    	   		for (int i=0;i<n;i++)
    	   			deriv(i,j)=(f[i]-fvec[i])/h;
    	   		}
    	return deriv;
    }

    //! @brief ... the routine performs line search correction of the Newton step
// NR    template <class T>
    void LineSearch(const NuTo::FullVector<double,Eigen::Dynamic> &rParameter,
    		NuTo::FullVector<double,Eigen::Dynamic> &xold,
    		const double fold,
    		NuTo::FullVector<double,Eigen::Dynamic> &g,
    		NuTo::FullVector<double,Eigen::Dynamic> &p,
    		NuTo::FullVector<double,Eigen::Dynamic> &x,
// NR    		double &f, const double stpmax, bool &check, T &func) {
		double &f, const double stpmax, bool &check, NuTo::FullVector<double,Eigen::Dynamic> &fvec) const {  // AnstattNR
    	const double ALF=1.0e-4, TOLX=numeric_limits<double>::epsilon();
    	double a,alam,alam2=0.0,alamin,b,disc,f2=0.0;
    	double rhs1,rhs2,slope=0.0,sum=0.0,test,tmplam;

    	check = false;
    //                                              // OPTIMIZED
    //	for (i=0;i<n;i++) sum += p[i]*p[i];   		// OPTIMIZED
    //	sum=sqrt(sum);                        		// OPTIMIZED
    	sum = p.lpNorm<2>();                  		// OPTIMIZED

    	if (sum > stpmax)
    //  											// OPTIMIZED
    //		for (i=0;i<n;i++)                       // OPTIMIZED
    //			p[i] *= stpmax/sum;                 // OPTIMIZED
    		p *= stpmax/sum;                        // OPTIMIZED

    //                                              // OPTIMIZED
    //	for (i=0;i<n;i++)                           // OPTIMIZED
    //		slope += g[i]*p[i];                     // OPTIMIZED
    	slope = g.dot(p);	                        // OPTIMIZED

    	if (slope >= 0.0){
    		cout << "Roundoff problem in lnsrch." << endl;
    		throw("Roundoff problem in lnsrch.");
    	}
    	test=0.0;
    //                 								// OPTIMIZED
    //	for (i=0;i<n;i++) {
    //		temp=fabs(p[i])/max(fabs(xold[i]),1.0); // OPTIMIZED
    //		if (temp > test) test=temp;             // OPTIMIZED
    //	}    										// OPTIMIZED
    	test = ( p.array().abs() / xold.array().abs().max(1.0) ).maxCoeff();  // OPTIMIZED

    	alamin=TOLX/test;
    	alam=1.0;
    	for (;;) {
    // 												// OPTIMIZED
    //		for (i=0;i<n;i++) x[i]=xold[i]+alam*p[i];  // OPTIMIZED
    		x = xold + alam*p; 						// OPTIMIZED
    // NR		f=func(x);
    		f = this->Fmin(rParameter,x,fvec);
    		if (alam < alamin) {
    //												// OPTIMIZED
    //			for (i=0;i<n;i++) x[i]=xold[i]; 	// OPTIMIZED
    			x = xold;							// OPTIMIZED
    			check=true;
    			return;
    		} else if (f <= fold+ALF*alam*slope) return;
    		else {
    			if (alam == 1.0)
    				tmplam = -slope/(2.0*(f-fold-slope));
    			else {
    				rhs1=f-fold-alam*slope;
    				rhs2=f2-fold-alam2*slope;
    				a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
    				b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
    				if (a == 0.0) tmplam = -slope/(2.0*b);
    				else {
    					disc=b*b-3.0*a*slope;
    					if (disc < 0.0) tmplam=0.5*alam;
    					else if (b <= 0.0) tmplam=(-b+sqrt(disc))/(3.0*a);
    					else tmplam=-slope/(b+sqrt(disc));
    				}
    				if (tmplam>0.5*alam)
    					tmplam=0.5*alam;
    			}
    		}
    		alam2=alam;
    		f2 = f;
    		alam=max(tmplam,0.1*alam);
    	}
    }

    //! @brief ... the routine performs Newton-Raphson integration
//NR    template <class T>
    void Newton(const NuTo::FullVector<double,Eigen::Dynamic>& rParameter,
    		NuTo::FullVector<double,Eigen::Dynamic> &x, bool &check,
//NR    		T &vecfunc,
   		NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>
    	(DamageViscoPlasticityEngineeringStress::*fdjacAn)(const NuTo::FullVector<double,Eigen::Dynamic>&,
    			NuTo::FullVector<double,Eigen::Dynamic>) const = 0) const{
    	const int MAXITS=200;
    	const double TOLF = mTOLF, TOLMIN=1.0e-12, STPMX=100.0;
    	const double TOLX=numeric_limits<double>::epsilon();
    	int its,n=x.GetNumRows();
    	double den,f,fold,stpmax,test;
    	NuTo::FullVector<double,Eigen::Dynamic> g(n),p(n),xold(n);
    	NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> fjac(n,n);
//NR    	NRfmin<T> fmin(vecfunc);
//NR    	NRfdjac<T> fdjac(vecfunc);
//NR    	NuTo::FullVector<double,Eigen::Dynamic> &fvec=fmin.fvec;
    	NuTo::FullVector<double,Eigen::Dynamic> fvec;   // AnstattNR
//NR    	f=fmin(x);
    	f = this->Fmin(rParameter, x, fvec); // AnstattNR
    //												// OPTIMIZED
    //	test=0.0;
    //	for (i=0;i<n;i++)							// OPTIMIZED
    //		if (fabs(fvec[i]) > test) test=fabs(fvec[i]);   // OPTIMIZED
    	test = fvec.array().abs().maxCoeff();		// OPTIMIZED
    	if (test < 0.01*TOLF) {
    		check=false;
    		return;
    	}
    //												// OPTIMIZED
    //	sum=0.0;									// OPTIMIZED
    //	for (i=0;i<n;i++) sum += x[i]*x[i];			// OPTIMIZED
    //	stpmax=STPMX*max(sqrt(sum),double(n));		// OPTIMIZED
    	stpmax=STPMX*max(x.lpNorm<2>(),double(n));	// OPTIMIZED

    	for (its=0;its<MAXITS;its++) {
    		cout<< "===== Iteration =====" << its <<endl;   // Test
    		if (fdjacAn != 0) {         			// if analytical Jacobi is there
    			fjac = (this->*fdjacAn)(rParameter,x);     	// analytical Jacobi matrix
//    			cout<<"*** Analytical ***"<<endl;  	// Test
//    			cout << fjac << endl;             	// Test
    		} else {                   				// if not analytic -> take numeric
    			fjac=this->DResidualNum(rParameter,x,fvec); // numerical Jacobi matrix
//    			cout<<"*** Numerical ***"<<endl;   	// Test
//    			cout << fjac << endl;	           	// Test
    		}
    //												// OPTIMIZED
    //		for (i=0;i<n;i++) {						// OPTIMIZED
    //			sum=0.0;							// OPTIMIZED
    //			for (j=0;j<n;j++) sum += fjac(j,i)*fvec[j];   // OPTIMIZED
    //			g[i]=sum;							// OPTIMIZED
    //		}										// OPTIMIZED
    		g = fjac.transpose() * fvec;			// OPTIMIZED
    //												// OPTIMIZED
    //		for (i=0;i<n;i++) xold[i]=x[i];         // OPTIMIZED
    		xold = x;
    		fold=f;
    // 												// OPTIMIZED
    //		for (i=0;i<n;i++) p[i] = -fvec[i];      // OPTIMIZED
    		p = -fvec;								// OPTIMIZED
    //												// Test
    		cout << "fvec before solve= " << fvec.transpose()<<endl;   // Test
    		cout << "p = " << fjac.fullPivLu().solve(-fvec).transpose()<<endl;   // Test

    //												// OPTIMIZED
    //		LUdcmp alu(fjac);						// OPTIMIZED
    //		alu.solve(p,p)							// OPTIMIZED
    		p =	fjac.fullPivLu().solve(-fvec).transpose();		// LU SOLVER of fjac * p = -fvec
    //															// SVD SOLVER
    //		p = fjac.jacobiSvd().solve(-fvec).transpose();      // SVD SOLVER
//NR    		LineSearch(xold,fold,g,p,x,f,stpmax,check,fmin);
    		LineSearch(rParameter,xold,fold,g,p,x,f,stpmax,check,fvec);    // AnstattNR
    		cout << "fvec after solve= " << fvec.transpose()<<endl;   // Test
    //												// OPTIMIZED
    //		test=0.0;								// OPTIMIZED
    //		for (i=0;i<n;i++)						// OPTIMIZED
    //			if (fabs(fvec[i]) > test) test=fabs(fvec[i]);   // OPTIMIZED
    		test = fvec.array().abs().maxCoeff();	// OPTIMIZED

    		if (test < TOLF) {
    			check=false;
    			cout << "iter = " << its << " ORDENARY RETURN, test = " << test << " TOLF = " << TOLF << endl;
    			return;
    		}

    		if (check) {
    			den=max(f,0.5*n);
    //												// OPTIMIZED
    //			test=0.0;							// OPTIMIZED
    //			for (i=0;i<n;i++) {					// OPTIMIZED
    //				temp=fabs(g[i])*max(fabs(x[i]),1.0)/den;   // OPTIMIZED
    //				if (temp > test) test=temp;		// OPTIMIZED
    //			}									// OPTIMIZED
    			test = ( g.array().abs() * x.array().abs().max(1.0) ).maxCoeff() / den;   // OPTIMIZED
    			check=(test < TOLMIN);
    			return;
    		}
    //												// OPTIMIZED
    //		test=0.0;								// OPTIMIZED
    //		for (i=0;i<n;i++) {						// OPTIMIZED
    //			temp=(fabs(x[i]-xold[i]))/max(fabs(x[i]),1.0);   // OPTIMIZED
    //			if (temp > test) test=temp;			// OPTIMIZED
    //		}										// OPTIMIZED
    		test = ( (x.array()-xold.array()).abs() / x.array().abs().max(1.0) ).maxCoeff();   // OPTIMIZED

    		if (test < TOLX)
    			return;
    	}
    	throw("MAXITS exceeded in newt");
    }

    //! @brief ... calculates the numerical algorithmic tangent:= derivative of stress by the strain
    //! @brief ... rParameter the list of parameters, containing rEngineeringStrain engineering strain at the end of time increment
    //! @brief ... rUnknown, state variables (unknowns), which correspond to rParameter. rUnknown contains stress at the end of the time increment
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> DStressDEpsNum(const NuTo::FullVector<double,Eigen::Dynamic> &rParameter,
    		NuTo::FullVector<double,Eigen::Dynamic> rUnknown, const int rDimension) const {
    	const double EPS = 1.0e-8;
    	int n=rDimension;  // 3D n = 6, 2D n = 4, 1D n = 1;
		EngineeringStrain3D rDeltaStrain;   // mechanical strain increment (strain increment without thermal component)
		EngineeringStress3D rDeltaStress;	// stress increment respective to rDeltaStrain
    	NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> deriv(n,n);

		// initialization of rDeltaStrain and rDeltaStress from the vectors rParameter and rUnknown respectively
		for (int i = 0; i < n; i++) {
			rDeltaStrain[i] = rParameter[i];	// rParameter(0:n-1) increment of mechanical strain
			rDeltaStress[i] = rUnknown[i];		// rUnknown(0:n-1) stress increment
		}

		NuTo::FullVector<double,Eigen::Dynamic> rParameterPert(rParameter); // perturbated strain increment
		NuTo::FullVector<double,Eigen::Dynamic> rUnknownPert(rUnknown);    	// respective perturbated stress increment
		EngineeringStress3D rDeltaStressPert;								// respective perturbated stress increment

    	for (int j=0;j<n;j++) {
    		double temp = rDeltaStrain[j];
    	   	double h = EPS;

    	   	rDeltaStrain[j] = temp+h;				// perturbation of strain
    	   	h = rDeltaStrain[j] - temp;				// remove the roundoff error
    	   	rParameterPert[j] = rDeltaStrain[j];	// update rParameter to rParameterPert

    	   	// prepare starting Newton
    	   	bool check;
            NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> (DamageViscoPlasticityEngineeringStress::*fdjacAn)
            	(const NuTo::FullVector<double,Eigen::Dynamic>&,NuTo::FullVector<double,Eigen::Dynamic>) const;

            // set Jacobi to analytical Jacobi
            fdjacAn = &DamageViscoPlasticityEngineeringStress::DResidualAn;

            this->Newton(rParameterPert,rUnknownPert,check,fdjacAn);

    		for (int i = 0; i < n; i++) {
    			rDeltaStressPert[i] = rUnknownPert[i];		// rUnknownPert(0:n-1), perturbated stress increment
    		}

    	   	for (int i=0;i<n;i++) {
    	   		deriv(i,j)=(rDeltaStressPert[i]-rDeltaStress[i])/h;
    	   	}

    	   	rDeltaStrain[j] = temp;
    	   	rParameterPert[j] = rDeltaStrain[j];
    	   	rUnknownPert = rUnknown;
    	}
    	return deriv;
    }

};

}

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::DamageViscoPlasticityEngineeringStress)
#endif //ENABLE_SERIALIZATION

#endif // DamageViscoPlasticityEngineeringStress_H_
