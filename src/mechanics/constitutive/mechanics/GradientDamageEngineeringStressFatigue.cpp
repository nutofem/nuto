// $Id: GradientDamageEngineeringStressFatigue.cpp 612 2012-08-13 07:31:23Z unger3 $
// GradientDamageEngineeringStressFatigue.cpp
// created Apr 26, 2010 by Joerg F. Unger
#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include <eigen3/Eigen/LU>

#include "base/Logger.h"
#include "mechanics/MechanicsException.h"
#include "mechanics/structures/StructureBase.h"
#include "mechanics/constitutive/ConstitutiveBase.h"
#include "mechanics/constitutive/ConstitutiveStaticDataBase.h"
#include "mechanics/constitutive/ConstitutiveTangentNonlocal.h"
#include "mechanics/constitutive/mechanics/ConstitutiveStaticDataGradientDamage1DFatigue.h"
#include "mechanics/constitutive/mechanics/ConstitutiveStaticDataGradientDamage2DFatigue.h"
#include "mechanics/constitutive/mechanics/Damage.h"
#include "mechanics/constitutive/mechanics/DeformationGradient1D.h"
#include "mechanics/constitutive/mechanics/DeformationGradient2D.h"
#include "mechanics/constitutive/mechanics/DeformationGradient3D.h"
#include "mechanics/constitutive/mechanics/GradientDamageEngineeringStressFatigue.h"
#include "mechanics/constitutive/mechanics/EngineeringStress1D.h"
#include "mechanics/constitutive/mechanics/EngineeringStress2D.h"
#include "mechanics/constitutive/mechanics/EngineeringStress3D.h"
#include "mechanics/constitutive/mechanics/EngineeringStrain1D.h"
#include "mechanics/constitutive/mechanics/EngineeringStrain2D.h"
#include "mechanics/constitutive/mechanics/EngineeringStrain3D.h"
#include "mechanics/constitutive/mechanics/LocalEqStrain.h"
#include "mechanics/constitutive/mechanics/NonlocalEqStrain.h"
#include "mechanics/elements/ElementBase.h"
#include "mechanics/sections/SectionBase.h"
#include "mechanics/sections/SectionEnum.h"
#include "optimize/NewtonRaphson.h"


#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <functional>


#define MAX_OMEGA 0.999
//#define ENABLE_DEBUG

NuTo::GradientDamageEngineeringStressFatigue::GradientDamageEngineeringStressFatigue() :
        ConstitutiveBase(), mRho(0.), mE(0.), mNu(0.), mNonlocalRadius(0.), mNonlocalRadiusParameter(0.), mThermalExpansionCoefficient(0.), mTensileStrength(0.),
		mCompressiveStrength(0.), mFractureEnergy(0.), mViscosityExponent(1.), mDamageLawType(Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING)
{
    SetParametersValid();

#ifdef ENABLE_DEBUG
    std::cout << "NuTo::GradientDamageEngineeringStressFatigue::GradientDamageEngineeringStressFatigue debug active" << std::endl;
#else
    std::cout << "NuTo::GradientDamageEngineeringStressFatigue::GradientDamageEngineeringStressFatigue debug inactive" << std::endl;
#endif

}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template<class Archive>
void NuTo::GradientDamageEngineeringStressFatigue::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize GradientDamageEngineeringStressFatigue" << "\n";
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveBase)
    & BOOST_SERIALIZATION_NVP(mRho)
    & BOOST_SERIALIZATION_NVP(mE)
    & BOOST_SERIALIZATION_NVP(mNu)
    & BOOST_SERIALIZATION_NVP(mThermalExpansionCoefficient)
    & BOOST_SERIALIZATION_NVP(mNonlocalRadius)
    & BOOST_SERIALIZATION_NVP(mNonlocalRadiusParameter)
    & BOOST_SERIALIZATION_NVP(mThermalExpansionCoefficient)
    & BOOST_SERIALIZATION_NVP(mTensileStrength)
    & BOOST_SERIALIZATION_NVP(mCompressiveStrength)
    & BOOST_SERIALIZATION_NVP(mViscosityExponent)
    & BOOST_SERIALIZATION_NVP(mFractureEnergy)
    & BOOST_SERIALIZATION_NVP(mDamageLawType)
    & BOOST_SERIALIZATION_NVP(mDamageLawParameters);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize GradientDamageEngineeringStressFatigue" << "\n";
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::GradientDamageEngineeringStressFatigue)
#endif // ENABLE_SERIALIZATION

//! @brief ... evaluate the constitutive relation in 1D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
//! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
NuTo::Error::eError NuTo::GradientDamageEngineeringStressFatigue::Evaluate1D(ElementBase* rElement, int rIp,
        const std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*>& rConstitutiveInput,
        std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>& rConstitutiveOutput)
{
    // get section information determining which input on the constitutive level should be used
    const SectionBase* section(rElement->GetSection());
	const StructureBase* structure(rElement->GetStructure());

    // check if parameters are valid
    if (this->mParametersValid == false)
    {
        //throw an exception giving information related to the wrong parameter
        CheckParameters();
    }

    if (section->GetType() != Section::TRUSS)
        throw MechanicsException("[NuTo::GradientDamageEngineeringStressFatigue::Evaluate1D] only truss sections are implemented.");

    // get engineering strain
    EngineeringStrain1D engineeringStrain1D;
	if (rConstitutiveInput.find(NuTo::Constitutive::Input::ENGINEERING_STRAIN_1D)==rConstitutiveInput.end()) {
		if (rConstitutiveInput.find(NuTo::Constitutive::Input::DEFORMATION_GRADIENT_1D)==rConstitutiveInput.end()) {
			throw MechanicsException("[NuTo::GradientDamageEngineeringStressFatigue::Evaluate1D] engineering strain 1d or deformation gradient 1d needed.");
		} else {
			const DeformationGradient1D& deformationGradient(rConstitutiveInput.find(NuTo::Constitutive::Input::DEFORMATION_GRADIENT_1D)->second->GetDeformationGradient1D());
			deformationGradient.GetEngineeringStrain(engineeringStrain1D);
		}
	} else {
		engineeringStrain1D = rConstitutiveInput.find(NuTo::Constitutive::Input::ENGINEERING_STRAIN_1D)->second->GetEngineeringStrain1D();
	}

    // get nonlocal eq strain
    if (rConstitutiveInput.find(NuTo::Constitutive::Input::NONLOCAL_EQ_STRAIN) == rConstitutiveInput.end())
        throw MechanicsException("[NuTo::GradientDamageEngineeringStressFatigue::Evaluate1D] nonlocal eq strain needed to evaluate engineering stress 1d.");
    const double nonlocalEqStrain = rConstitutiveInput.find(NuTo::Constitutive::Input::NONLOCAL_EQ_STRAIN)->second->GetNonlocalEqStrain().GetValue(0);

	// check whether Return Mapping should be done: check output
	bool performReturnMapping(false);

	if (rConstitutiveOutput.find(NuTo::Constitutive::Output::ENGINEERING_STRESS_1D)!=rConstitutiveOutput.end())
	{
		performReturnMapping = true;
	}
	if (rConstitutiveOutput.find(NuTo::Constitutive::Output::ENGINEERING_STRESS_3D)!=rConstitutiveOutput.end())
	{
		performReturnMapping = true;
	}
	if (rConstitutiveOutput.find(NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_1D)!=rConstitutiveOutput.end())
	{
		performReturnMapping = true;
	}
	if (rConstitutiveOutput.find(NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_NONLOCAL_EQ_STRAIN_1D)!=rConstitutiveOutput.end())
	{
		performReturnMapping = true;
	}
	if (rConstitutiveOutput.find(NuTo::Constitutive::Output::UPDATE_STATIC_DATA)!=rConstitutiveOutput.end())
	{
		performReturnMapping = true;
	}

	// perform a nonlinear iteration of Return Mapping
	// declare output of Return Mapping: stress, algorithmic tangent and static data
	EngineeringStress1D engineeringStress1D;
	ConstitutiveTangentLocal<1, 1> tangent,tangentNonlocal;
	ConstitutiveStaticDataGradientDamage1DFatigue newStaticData;

    if (performReturnMapping)
    {
    	// perform return mapping
    	NuTo::Error::eError errorReturnMapping = ReturnMapping1D(rElement, rIp, nonlocalEqStrain, engineeringStrain1D, &engineeringStress1D, &tangent, &tangentNonlocal, &newStaticData ,rElement->GetStructure()->GetLogger());
    	if (errorReturnMapping!=Error::SUCCESSFUL)
    		return errorReturnMapping;
    }

    // calculate output
    for (std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>::iterator itOutput = rConstitutiveOutput.begin(); itOutput != rConstitutiveOutput.end(); itOutput++)
    {
        switch (itOutput->first)
        {
        case NuTo::Constitutive::Output::ENGINEERING_STRESS_1D:
    	{
   			itOutput->second->GetEngineeringStress1D() = (1. - newStaticData.mOmega)*engineeringStress1D;
    	}
    		break;
        case NuTo::Constitutive::Output::ENGINEERING_STRESS_3D:
        {
            //this is for the visualize routines
            EngineeringStress3D& engineeringStress3D(itOutput->second->GetEngineeringStress3D());

            engineeringStress3D[0] = (1. - newStaticData.mOmega)*engineeringStress1D[0];
            engineeringStress3D[1] = 0.;
            engineeringStress3D[2] = 0.;
            engineeringStress3D[3] = 0.;
            engineeringStress3D[4] = 0.;
            engineeringStress3D[5] = 0.;
        }
            break;
        case NuTo::Constitutive::Output::ENGINEERING_STRAIN_3D:
        {
            //this is for the visualize routines
            EngineeringStrain3D& engineeringStrain3D(itOutput->second->GetEngineeringStrain3D());
            engineeringStrain3D[0] = engineeringStrain1D[0];
            engineeringStrain3D[1] = -mNu * engineeringStrain1D[0];
            engineeringStrain3D[2] = -mNu * engineeringStrain1D[0];
            engineeringStrain3D[3] = 0.;
            engineeringStrain3D[4] = 0.;
            engineeringStrain3D[5] = 0.;
        }
            break;
        case NuTo::Constitutive::Output::LOCAL_EQ_STRAIN:
        {
            itOutput->second->GetLocalEqStrain() = engineeringStrain1D.Norm();
        }
            break;
        case NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_1D:
        {
    	    itOutput->second->AsConstitutiveTangentLocal_1x1() = tangent;
    	    itOutput->second->AsConstitutiveTangentLocal_1x1().SetSymmetry(true);
        }
            break;
        case NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_NONLOCAL_EQ_STRAIN_1D:
        {
    	    itOutput->second->AsConstitutiveTangentLocal_1x1() = tangentNonlocal;
    	    itOutput->second->AsConstitutiveTangentLocal_1x1().SetSymmetry(true);
        }
            break;
        case NuTo::Constitutive::Output::D_LOCAL_EQ_STRAIN_D_STRAIN_1D:
        {
        	itOutput->second->AsConstitutiveTangentLocal_1x1() = engineeringStrain1D.DNorm();
        }
            break;
        case NuTo::Constitutive::Output::D_LOCAL_EQ_STRAIN_XI_D_STRAIN_1D:
        {
            ConstitutiveTangentLocal<1, 1>& tangent(itOutput->second->AsConstitutiveTangentLocal_1x1());
            tangent.SetSymmetry(true);

            double xi = mNonlocalRadius;
            double dXi = 0;

            if (mNonlocalRadiusParameter != 0.)
            {
                double e_xi = mTensileStrength / mE * mNonlocalRadiusParameter;
                if (engineeringStrain1D.Norm()[0] <= e_xi)
                {
                    double c0 = 0.01;
                    xi = c0 + (mNonlocalRadius - c0) * (engineeringStrain1D.Norm()[0] / e_xi);
                    dXi = (mNonlocalRadius - c0) / e_xi;
                }
            }

            double factor = 1. / xi - (engineeringStrain1D.Norm()[0] - nonlocalEqStrain) / (xi * xi) * dXi;

            tangent = factor * engineeringStrain1D.DNorm();
        }
            break;
        case NuTo::Constitutive::Output::NONLOCAL_PARAMETER_XI:
        {
            ConstitutiveTangentLocal<1, 1>& variableNonlocalRadius(itOutput->second->AsConstitutiveTangentLocal_1x1());
            variableNonlocalRadius[0] = mNonlocalRadius;
            if (mNonlocalRadiusParameter != 0.)
            {
                double e_xi = mTensileStrength / mE * mNonlocalRadiusParameter;
                if (engineeringStrain1D.Norm()[0] <= e_xi)
                {
                    double c0 = 0.01;
                    variableNonlocalRadius[0] = c0 + (mNonlocalRadius - c0) * (engineeringStrain1D.Norm()[0] / e_xi);
                }
            }
        }
            break;
        case NuTo::Constitutive::Output::DAMAGE:
        {
            itOutput->second->GetDamage().SetDamage(newStaticData.mOmega);
        }
            break;
        case NuTo::Constitutive::Output::UPDATE_TMP_STATIC_DATA:
        {
            throw MechanicsException("[NuTo::GradientDamageEngineeringStressFatigue::Evaluate1D] tmp_static_data has to be updated without any other outputs, call it separately.");
        }
            break;
        case NuTo::Constitutive::Output::UPDATE_STATIC_DATA:
        {
    	    *(rElement->GetStaticData(rIp)->AsGradientDamage1DFatigue()) = newStaticData;
    	    rElement->GetStaticData(rIp)->AsGradientDamage1DFatigue()->mPrevStrain = engineeringStrain1D;
    	    rElement->GetStaticData(rIp)->AsGradientDamage1DFatigue()->mPrevSigma = engineeringStress1D;

    	    // Damage
    	    std::ofstream DamageFile;
    	    DamageFile.open("Damage.txt", std::ios::app);

//    	    DamageFile << " time     omega     kappa     nonLocalEqStrain" << std::endl;
            if (rElement->ElementGetId() == 49 && rIp == 0) {
            	DamageFile << structure->GetTime() << " " << newStaticData.mOmega << " " << newStaticData.mKappa <<
            		" " << newStaticData.mPrevNonlocalEqStrain.GetValue(0) << std::endl;
            }
        }
            break;
        default:
            throw MechanicsException(
                    std::string("[NuTo::GradientDamageEngineeringStressFatigue::Evaluate1D] output object ") + NuTo::Constitutive::OutputToString(itOutput->first)
                            + std::string(" could not be calculated, check the allocated material law and the section behavior."));
        }
    }

    return Error::SUCCESSFUL;
}

//! @brief ... performs the return mapping procedure in 1D
//! @param rElement ... structure
//! @param rIp ... integration point
//! @param rNonlocalEqStrain ... nonlocal eq strain
//! @param rEngineeringStrain ... engineering strain
//! @param rNewStress ... new stress (if a 0-pointer is given, no values are written)
//! @param rNewTangent ... new tangent matrix (if a 0-pointer is given, no values are written)
//! @param rNewTangentNonlocal ... new tangent matrix (if a 0-pointer is given, no values are written) derivative of stress by the nonlocal eq strain
//! @param rNewStaticData ... new static data (if a 0-pointer is given, no values are written)
NuTo::Error::eError NuTo::GradientDamageEngineeringStressFatigue::ReturnMapping1D(const ElementBase* rElement,int rIp,
		const double rNonlocalEqStrain,
		const EngineeringStrain1D& rEngineeringStrain,
		EngineeringStress1D* rNewStress,
		ConstitutiveTangentLocal<1,1>* rNewTangent,
		ConstitutiveTangentLocal<1,1>* rNewTangentNonlocal,
		ConstitutiveStaticDataGradientDamage1DFatigue* rNewStaticData,
		Logger& rLogger) const
	{
	    // get state variables to begin of the time increment
    	const ConstitutiveStaticDataGradientDamage1DFatigue *oldStaticData = rElement->GetStaticData(rIp)->AsGradientDamage1DFatigue();

	    // extract total mechanical strain to begin of the time increment
    	const double prevNonlocalEqStrain(oldStaticData->mPrevNonlocalEqStrain.GetValue(0));
    	const EngineeringStress1D prevStress(oldStaticData->GetPrevStress());
	    const EngineeringStrain1D prevStrain(oldStaticData->GetPrevStrain());  // prevStrain is the mechanical strain, without thermal strain

	    // calculate damage and kappa
	    double e_0 = mTensileStrength / mE;
	    double kappa(std::max(oldStaticData->mKappa, e_0 + e_0*1.0e-5)), omega;


	    if (rNonlocalEqStrain - kappa >= 0.) {
	    	// calculate omega and kappa
	    	kappa = rNonlocalEqStrain;
	    	omega = CalculateDamage(kappa);

	    	// calculate nonlocal eq tangent
	    	if (rNewTangentNonlocal!=0) {
	    		(*rNewTangentNonlocal)(0, 0) = - mE * rEngineeringStrain[0] * CalculateDerivativeDamage(kappa) * Sgn(rNonlocalEqStrain - prevNonlocalEqStrain);
	    	}
	    } else {
	    	// call Newton and calculate omega and kappa

	    	// initialize the vector of unknowns
	    	Eigen::VectorXd Unknown(2);
	    	double deltaOmega;

	    	deltaOmega = CalculateDerivativeDamage(kappa) * pow(rNonlocalEqStrain/kappa,mViscosityExponent) * Plus(rNonlocalEqStrain - prevNonlocalEqStrain);

	    	Unknown[0] = deltaOmega;									// delta omega
	    	Unknown[1] = deltaOmega/CalculateDerivativeDamage(kappa); 	// delta kappa

	    	// compose vector of known Parameter, which are necessary for ResidualAn
	    	Eigen::VectorXd Parameter(4);

	    	Parameter[0] = oldStaticData->mOmega;									// omega at the beginning of the time increment
	    	Parameter[1] = std::max(oldStaticData->mKappa, e_0 + e_0*1.0e-5);		// kappa at the beginning of the time increment
	    	Parameter[2] = prevNonlocalEqStrain;									// nonlocal eq strain at the beginning of the time increment
	    	Parameter[3] = rNonlocalEqStrain;										// current nonlocal eq strain

	    	const Eigen::VectorXd ParameterList(Parameter);

	        // prepare starting Newton solver with respect to the "Unknown"
//	        NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> (GradientDamageEngineeringStressFatigue::*fdjacAn)
//	        	(const Eigen::VectorXd&,Eigen::VectorXd) const;
//
//	        Eigen::VectorXd (GradientDamageEngineeringStressFatigue::*residual)
//	        		(const Eigen::VectorXd&,Eigen::VectorXd) const;
//
//	        // set Jacobi to analytical Jacobi
//	        fdjacAn = &GradientDamageEngineeringStressFatigue::DResidualAn;
//
//	        // set residual
//	        residual = &GradientDamageEngineeringStressFatigue::Residual;

	        boost::function<Eigen::VectorXd (const Eigen::VectorXd&,Eigen::VectorXd)> residualFunction;
	        residualFunction = boost::bind( &GradientDamageEngineeringStressFatigue::Residual, this, _1, _2 );

	        boost::function<NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>
	        	(const Eigen::VectorXd&,Eigen::VectorXd)> jacobiFunction;
	        jacobiFunction = boost::bind( &GradientDamageEngineeringStressFatigue::DResidualAn, this, _1, _2 );

	        // start Newton solver
	        try
	        {
#ifdef ENABLE_OPTIMIZE
	        	NuTo::NewtonRaphson myNonlinearSolver;
	        	// provide the solver with the pointers to the (!) residual, (!!) jacobi and (!!!) parameters necessary to evaluate the  both
	        	myNonlinearSolver.SetParameters(Parameter);
	        	myNonlinearSolver.SetResidualFunction(residualFunction);
	        	myNonlinearSolver.SetResidualDerivativeFunction(jacobiFunction);
//	        	myNonlinearSolver.SetResidualFunction(residual);
//	        	myNonlinearSolver.SetResidualDerivativeFunction(fdjacAn);

	        	// solve
	        	myNonlinearSolver.Solve(Unknown);
#else
    std::cout << "OPTIMIZE module not enabled - NewtonRaphson solver not available." << std::endl;
#endif // ENABLE_OPTIMIZE

//	        	// is the solution spurious?
//		        bool check;
//	        	check = myNonlinearSolver.GetCheckNewtonRaphson();
	        }
	        catch (...)
	        {
	            rLogger << "[NuTo::GradientDamageEngineeringStressFatigue::ReturnMapping1D] No convergence in Newton." << "\n";
	        }

	        // provide omega and kappa
	        omega  = oldStaticData->mOmega + Unknown[0];
	        kappa += Unknown[1];

	    	// calculate nonlocal eq tangent
	        if (rNewTangentNonlocal!=0) {
	        	NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> matrixMultipl;

	        	matrixMultipl  = ((GradientDamageEngineeringStressFatigue::DResidualAn(ParameterList,Unknown)).fullPivLu().
	        			solve(GradientDamageEngineeringStressFatigue::DResidualDEpsNonlocalAn(ParameterList,Unknown))).eval();
	        	matrixMultipl *= mE * rEngineeringStrain[0];
	        	(*rNewTangentNonlocal) = matrixMultipl.block<1,1>(0,0);
	        }
	    }

	    // calculate stress
	    if (rNewStress!=0) {
	    	(*rNewStress) = mE * rEngineeringStrain;
	    }
	    // calculate local tangent
	    if (rNewTangent!=0) {
	    	(*rNewTangent)(0,0) = (1. - omega) * mE;
	    }
	    // update statevs
	    if (rNewStaticData!=0) {
	    	rNewStaticData->mKappa = kappa;
	    	rNewStaticData->mOmega = omega;
	    	rNewStaticData->mPrevNonlocalEqStrain.SetValue(0,rNonlocalEqStrain);
	    	rNewStaticData->mPrevStrain = rEngineeringStrain;
	    	rNewStaticData->mPrevSigma = mE * rEngineeringStrain;
	    }

		return Error::SUCCESSFUL;
	}

//! @brief ... evaluate the constitutive relation in 2D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
//! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
NuTo::Error::eError NuTo::GradientDamageEngineeringStressFatigue::Evaluate2D(ElementBase* rElement, int rIp,
        const std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*>& rConstitutiveInput,
        std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>& rConstitutiveOutput)
{

    // check if parameters are valid
    if (this->mParametersValid == false)
    {
        //throw an exception giving information related to the wrong parameter
        CheckParameters();
    }

    // get engineering strain
    EngineeringStrain2D engineeringStrain2D;
	if (rConstitutiveInput.find(NuTo::Constitutive::Input::ENGINEERING_STRAIN_2D)==rConstitutiveInput.end()) {
		if (rConstitutiveInput.find(NuTo::Constitutive::Input::DEFORMATION_GRADIENT_2D)==rConstitutiveInput.end()) {
			throw MechanicsException("[NuTo::GradientDamageEngineeringStressFatigue::Evaluate2D] engineering strain 2d or deformation gradient 2d needed.");
		} else {
			const DeformationGradient2D& deformationGradient(rConstitutiveInput.find(NuTo::Constitutive::Input::DEFORMATION_GRADIENT_2D)->second->GetDeformationGradient2D());
			deformationGradient.GetEngineeringStrain(engineeringStrain2D);
		}
	} else {
		engineeringStrain2D = rConstitutiveInput.find(NuTo::Constitutive::Input::ENGINEERING_STRAIN_2D)->second->GetEngineeringStrain2D();
	}

    // nonlocal eq strain
    if (rConstitutiveInput.find(NuTo::Constitutive::Input::NONLOCAL_EQ_STRAIN) == rConstitutiveInput.end())
        throw MechanicsException("[NuTo::GradientDamageEngineeringStressFatigue::Evaluate2D] nonlocal eq strain needed to evaluate engineering stress 2d.");
    const double nonlocalEqStrain = rConstitutiveInput.find(NuTo::Constitutive::Input::NONLOCAL_EQ_STRAIN)->second->GetNonlocalEqStrain().GetValue(0);

	// check whether Return Mapping should be done: check output
	bool performReturnMapping(false);

	if (rConstitutiveOutput.find(NuTo::Constitutive::Output::ENGINEERING_STRESS_2D)!=rConstitutiveOutput.end())
	{
		performReturnMapping = true;
	}
	if (rConstitutiveOutput.find(NuTo::Constitutive::Output::ENGINEERING_STRESS_3D)!=rConstitutiveOutput.end())
	{
		performReturnMapping = true;
	}
	if (rConstitutiveOutput.find(NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_2D)!=rConstitutiveOutput.end())
	{
		performReturnMapping = true;
	}
	// D_ENGINEERING_STRESS_D_NONLOCAL_EQ_STRAIN_2D is performed simultaneously with:
//				D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_2D			-> requires ReturnMapping
//				or
//				D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_ELASTIC_2D	-> does not require ReturnMapping
//		for the extrapolation D_ENGINEERING_STRESS_D_NONLOCAL_EQ_STRAIN_2D is zero; therefore it is important to check whether the extrapolation does occur
//	if (rConstitutiveOutput.find(NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_NONLOCAL_EQ_STRAIN_2D)!=rConstitutiveOutput.end())
//	{
//		performReturnMapping = true;
//	}
	if (rConstitutiveOutput.find(NuTo::Constitutive::Output::UPDATE_STATIC_DATA)!=rConstitutiveOutput.end())
	{
		performReturnMapping = true;
	}

	// perform a nonlinear iteration of Return Mapping
	// declare output of Return Mapping: stress, algorithmic tangent and static data
	EngineeringStress2D engineeringStress2D;
	EngineeringStress3D engineeringStress3D;
	ConstitutiveTangentLocal<3, 3> tangent;
	ConstitutiveTangentLocal<3, 1> tangentNonlocal;
	ConstitutiveStaticDataGradientDamage2DFatigue newStaticData;


    if (performReturnMapping)
    {
    	// perform return mapping
    	NuTo::Error::eError errorReturnMapping = ReturnMapping2D(rElement, rIp, nonlocalEqStrain, engineeringStrain2D, &engineeringStress2D, &engineeringStress3D,
    			&tangent, &tangentNonlocal, &newStaticData, rElement->GetStructure()->GetLogger());
    	if (errorReturnMapping!=Error::SUCCESSFUL)
    		return errorReturnMapping;
    }


    // calculate local eq strain and tangent
    LocalEqStrain localEqStrain;
    ConstitutiveTangentLocal<3, 1> localEqStrainTangent;

    // Modified Mises Strain
    switch (rElement->GetSection()->GetType())
    {
    case Section::PLANE_STRAIN:
        CalculateLocalEqStrainAndDerivativeModifiedMises2DPlaneStrain(engineeringStrain2D, localEqStrain, localEqStrainTangent);
        break;
    case Section::PLANE_STRESS:
        CalculateLocalEqStrainAndDerivativeModifiedMises2DPlaneStress(engineeringStrain2D, localEqStrain, localEqStrainTangent);
        break;

    default:
        throw MechanicsException("[NuTo::GradientDamageEngineeringStressFatigue::Evaluate2D] Invalid type of 2D section behaviour found!!!");
    }

    // get elastic matrix
	// calculate coefficients of the linear elastic material matrix
	double C11, C12, C33;

	switch (rElement->GetSection()->GetType())
	{
	case Section::PLANE_STRAIN:
		// calculate coefficients of the material matrix
		this->CalculateCoefficients3D(C11, C12, C33);
		break;
	case Section::PLANE_STRESS:
    // calculate coefficients of the material matrix
		this->CalculateCoefficients2DPlaneStress(C11, C12, C33);
		break;
	default:
		throw MechanicsException("[NuTo::GradientDamageEngineeringStressFatigue::Evaluate2D] Invalid type of 2D section behaviour found!!!");
	}

	NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> ElasticStiffness(3,3);
	ElasticStiffness << C11, C12,  0.,
						C12, C11,  0.,
						0.,  0.,  C33;

//    // Standard equivalent strain
//    localEqStrain = CalculateLocalEqStrain2D(engineeringStrain2D);
//    localEqStrainTangent = CalculateLocalEqStrainTangent2D(engineeringStrain2D);

    // calculate output
    for (std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>::iterator itOutput = rConstitutiveOutput.begin(); itOutput != rConstitutiveOutput.end(); itOutput++)
    {
        switch (itOutput->first)
        {
        case NuTo::Constitutive::Output::ENGINEERING_STRESS_2D:
        {
   			itOutput->second->GetEngineeringStress2D() = (1. - newStaticData.mOmega)*engineeringStress2D;
        }
        	break;
        case NuTo::Constitutive::Output::ENGINEERING_STRESS_3D:
        {
            //this is for the visualize routines
            itOutput->second->GetEngineeringStress3D() = (1. - newStaticData.mOmega)*engineeringStress3D;
        }
            break;
        case NuTo::Constitutive::Output::ENGINEERING_STRAIN_3D:
        {
            EngineeringStrain3D& engineeringStrain3D(itOutput->second->GetEngineeringStrain3D());
            switch (rElement->GetSection()->GetType())
            {
            case Section::PLANE_STRAIN:
                engineeringStrain3D[0] = engineeringStrain2D[0];
                engineeringStrain3D[1] = engineeringStrain2D[1];
                engineeringStrain3D[2] = 0;
                engineeringStrain3D[3] = engineeringStrain2D[2];
                engineeringStrain3D[4] = 0.;
                engineeringStrain3D[5] = 0.;
                break;
            case Section::PLANE_STRESS:
                engineeringStrain3D[0] = engineeringStrain2D[0];
                engineeringStrain3D[1] = engineeringStrain2D[1];
                engineeringStrain3D[2] = mNu / (mNu - 1.) * (engineeringStrain2D[0] + engineeringStrain2D[1]);
                engineeringStrain3D[3] = engineeringStrain2D[2];
                engineeringStrain3D[4] = 0.;
                engineeringStrain3D[5] = 0.;
                break;
            default:
                throw MechanicsException("[NuTo::GradientDamageEngineeringStressFatigue::ENGINEERING_STRAIN_3D] Invalid type of 2D section behaviour found!!!");
            }
        }
            break;
        case NuTo::Constitutive::Output::LOCAL_EQ_STRAIN:
        {
            itOutput->second->GetLocalEqStrain() = localEqStrain;
        }
            break;
        case NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_2D:
        {
    	    itOutput->second->AsConstitutiveTangentLocal_3x3() = tangent;
    	    itOutput->second->AsConstitutiveTangentLocal_3x3().SetSymmetry(true);
        }
            break;
        case NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_NONLOCAL_EQ_STRAIN_2D:
        {
        	// check whether an extrapolation occurs
        	if (rConstitutiveOutput.find(NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_ELASTIC_2D)!=rConstitutiveOutput.end()) {
        		// extrapolation, put zero
        		ConstitutiveTangentLocal<3, 1> zeroTangentNonlocal;
        		zeroTangentNonlocal << 0., 0.,  0.;
        		itOutput->second->AsConstitutiveTangentLocal_3x1() = zeroTangentNonlocal;
        		if (rElement->ElementGetId() == 189 && rIp == 0) {
        			std::cout << "tangentNonLocal ===============" << std::endl;
        			std::cout << itOutput->second->AsConstitutiveTangentLocal_3x1() << std::endl;
        		}
        	} else {
        		// conventional incremental integration, take from ReturnMapping
        		itOutput->second->AsConstitutiveTangentLocal_3x1() = tangentNonlocal;
        	}
    	    itOutput->second->AsConstitutiveTangentLocal_3x1().SetSymmetry(false);
        }
            break;
        case NuTo::Constitutive::Output::D_LOCAL_EQ_STRAIN_D_STRAIN_2D:
        {
            itOutput->second->AsConstitutiveTangentLocal_3x1() = localEqStrainTangent;
        }
            break;
        case NuTo::Constitutive::Output::D_LOCAL_EQ_STRAIN_XI_D_STRAIN_2D:
        {
            ConstitutiveTangentLocal<3, 1>& tangent(itOutput->second->AsConstitutiveTangentLocal_3x1());

            double xi = mNonlocalRadius;
            double dXi = 0;

            if (mNonlocalRadiusParameter != 0.)
            {
                double e_xi = mTensileStrength / mE * mNonlocalRadiusParameter;
                if (localEqStrain[0] <= e_xi)
                {
                    double c0 = 0.01;
                    xi = c0 + (mNonlocalRadius - c0) * (localEqStrain[0] / e_xi);
                    dXi = (mNonlocalRadius - c0) / e_xi;
                }
            }

            double factor = 1. / xi - (localEqStrain[0] - nonlocalEqStrain) / (xi * xi) * dXi;

            tangent = factor * localEqStrainTangent;
        }
            break;
        case NuTo::Constitutive::Output::NONLOCAL_PARAMETER_XI:
        {
            ConstitutiveTangentLocal<1, 1>& variableNonlocalRadius(itOutput->second->AsConstitutiveTangentLocal_1x1());
            variableNonlocalRadius[0] = mNonlocalRadius;
            if (mNonlocalRadiusParameter != 0.)
            {
                double e_xi = mTensileStrength / mE * mNonlocalRadiusParameter;
                if (localEqStrain[0] <= e_xi)
                {
                    double c0 = 0.01;
                    variableNonlocalRadius[0] = c0 + (mNonlocalRadius - c0) * (localEqStrain[0] / e_xi);
                }
            }
        }
            break;
        case NuTo::Constitutive::Output::DAMAGE:
        {
            itOutput->second->GetDamage().SetDamage(newStaticData.mOmega);
        }
        	break;
        case NuTo::Constitutive::Output::UPDATE_TMP_STATIC_DATA:
        {
            throw MechanicsException("[NuTo::GradientDamageEngineeringStressFatigue::Evaluate2D] tmp_static_data has to be updated without any other outputs, call it separately.");
        }
            break;
        case NuTo::Constitutive::Output::UPDATE_STATIC_DATA:
        {	// the statevs for storing (like mOmegaFatigue ...) should not be updated, cause they will be set to zero, what is wrong
//    	    *(rElement->GetStaticData(rIp)->AsGradientDamage2DFatigue()) = newStaticData;
    	    rElement->GetStaticData(rIp)->AsGradientDamage2DFatigue()->mPrevStrain = engineeringStrain2D;
    	    rElement->GetStaticData(rIp)->AsGradientDamage2DFatigue()->mPrevSigma = engineeringStress2D;
    	    rElement->GetStaticData(rIp)->AsGradientDamage2DFatigue()->mPrevNonlocalEqStrain.SetValue(0,nonlocalEqStrain);
    	    rElement->GetStaticData(rIp)->AsGradientDamage2DFatigue()->mOmega = newStaticData.mOmega;
    	    rElement->GetStaticData(rIp)->AsGradientDamage2DFatigue()->mKappa = newStaticData.mKappa;

    	    // Damage
//    	    std::ofstream DamageFile;
//    	    DamageFile.open("Damage.txt", std::ios::app);
//
//    	    DamageFile << " time     omega     kappa     nonLocalEqStrain" << std::endl;
//            if (rElement->ElementGetId() == 29566 && rIp == 0) {
//            	DamageFile << structure->GetTime() << " " << newStaticData.mOmega << " " << newStaticData.mKappa <<
//            		" " << newStaticData.mPrevNonlocalEqStrain.GetValue(0) << std::endl;
//            }
        }
            break;
    	case NuTo::Constitutive::Output::ENGINEERING_STRESS_ELASTIC_2D:
    	{
    		// get damage from static data
    	    const ConstitutiveStaticDataGradientDamage2DFatigue *ExtrapolatedStaticData = rElement->GetStaticData(rIp)->AsGradientDamage2DFatigue();

    		itOutput->second->GetEngineeringStress2D() = (1. - ExtrapolatedStaticData->mOmega) * ElasticStiffness * engineeringStrain2D;
    		if (rElement->ElementGetId() == 189 && rIp == 0) {
    			std::cout << "Elastic Stress ===============" << std::endl;
    			std::cout << itOutput->second->GetEngineeringStress2D() << std::endl;
    		}
    	}
    		break;
    	case NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_ELASTIC_2D:
    	{
    		// get damage from static data
    	    const ConstitutiveStaticDataGradientDamage2DFatigue *ExtrapolatedStaticData = rElement->GetStaticData(rIp)->AsGradientDamage2DFatigue();

    	    itOutput->second->AsConstitutiveTangentLocal_3x3() = (1. - ExtrapolatedStaticData->mOmega) * ElasticStiffness;
    	    itOutput->second->AsConstitutiveTangentLocal_3x3().SetSymmetry(true);
    		if (rElement->ElementGetId() == 189 && rIp == 0) {
    			std::cout << "Elastic Stiffness ===============" << std::endl;
    			std::cout << itOutput->second->AsConstitutiveTangentLocal_3x3() << std::endl;
    		}
    	}

    		break;
    	case NuTo::Constitutive::Output::FATIGUE_SAVE_STATIC_DATA:
    	{
    		rElement->GetStaticData(rIp)->AsGradientDamage2DFatigue()->FatigueSaveStaticData();
//    		std::cout<< "AsDVP3D = " << rElement->GetStaticData(rIp)->AsDamageViscoPlasticity3D() <<
//    		   		", damage value = " << rElement->GetStaticData(rIp)->AsDamageViscoPlasticity3D()->mPrevStrain <<
//    		   		"AsDVP3DFatigue = " << rElement->GetStaticData(rIp)->AsDamageViscoPlasticity3DFatigue() <<
//    		   		", damage value = " << rElement->GetStaticData(rIp)->AsDamageViscoPlasticity3DFatigue()->mPrevStrain <<
//    		  		", damageFatigue value = " << rElement->GetStaticData(rIp)->AsDamageViscoPlasticity3DFatigue()->mPrevStrainFatigue << std::endl;
    	}
		break;
    	case NuTo::Constitutive::Output::FATIGUE_RESTORE_STATIC_DATA:
    	{
    		rElement->GetStaticData(rIp)->AsGradientDamage2DFatigue()->FatigueRestoreStaticData();
    	}
		break;
    	case NuTo::Constitutive::Output::FATIGUE_EXTRAPOLATE_STATIC_DATA:
    	{
    		// get static data
    		ConstitutiveStaticDataGradientDamage2DFatigue *ExtrapolatedStaticData = rElement->GetStaticData(rIp)->AsGradientDamage2DFatigue();

    		if (rElement->GetStructure()->GetNumExtrapolatedCycles()[0] > 1) {
    			// make extrapolation of state variables, except of mPrevStrain, mPrevSigma and mPrevNonlocalStrain, because this should be calculated after the equilibrium is found
        		ExtrapolatedStaticData->FatigueExtrapolateStaticData(rElement->GetStructure()->GetNumExtrapolatedCycles());	// linear extrapolation of kappa
        		ExtrapolatedStaticData->SetOmega(this->CalculateDamage(ExtrapolatedStaticData->mKappa));					// calculate omega from kappa
    		} else {
				throw MechanicsException("[NuTo::DamageViscoPlasticityHardeningEngineeringStress::Evaluate3D] the number of cycles to be extrapolated mNumExtrapolatedCycles <= 1!");
			}
    	}
    	break;
        default:
            throw MechanicsException(
                    std::string("[NuTo::GradientDamageEngineeringStressFatigue::Evaluate2D] output object ") + NuTo::Constitutive::OutputToString(itOutput->first)
                            + std::string(" could not be calculated, check the allocated material law and the section behavior."));
        }
    }

    return Error::SUCCESSFUL;
}

//! @brief ... performs the return mapping procedure in 2D
//! @param rElement ... structure
//! @param rIp ... integration point
//! @param rNonlocalEqStrain ... nonlocal eq strain
//! @param rEngineeringStrain ... engineering strain
//! @param rNewStress2D, rNewStress3D ... new stress in 2D and 3D  (if a 0-pointer is given, no values are written)
//! @param rNewTangent ... new tangent matrix (if a 0-pointer is given, no values are written)
//! @param rNewTangentNonlocal ... new tangent matrix (if a 0-pointer is given, no values are written) derivative of stress by the nonlocal eq strain
//! @param rNewStaticData ... new static data (if a 0-pointer is given, no values are written)
NuTo::Error::eError NuTo::GradientDamageEngineeringStressFatigue::ReturnMapping2D(const ElementBase* rElement,int rIp,
		const double rNonlocalEqStrain,
		const EngineeringStrain2D& rEngineeringStrain,
		EngineeringStress2D* rNewStress2D,
		EngineeringStress3D* rNewStress3D,
		ConstitutiveTangentLocal<3,3>* rNewTangent,
		ConstitutiveTangentLocal<3,1>* rNewTangentNonlocal,
		ConstitutiveStaticDataGradientDamage2DFatigue* rNewStaticData,
		Logger& rLogger) const
	{
		// get elastic stiffness matrix
		double C11, C12, C33;

		switch (rElement->GetSection()->GetType())
    	{
    	case Section::PLANE_STRAIN:
    		// calculate coefficients of the material matrix
    		this->CalculateCoefficients3D(C11, C12, C33);
    		break;
    	case Section::PLANE_STRESS:
        // calculate coefficients of the material matrix
    		this->CalculateCoefficients2DPlaneStress(C11, C12, C33);
    		break;
    	default:
    		throw MechanicsException("[NuTo::GradientDamageEngineeringStressFatigue::ReturnMapping2D] Invalid type of 2D section behaviour found!!!");
    	}

		NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> ElasticStiffness(3,3);
		ElasticStiffness << C11, C12,  0.,
							C12, C11,  0.,
							0.,  0.,  C33;

	    // get state variables to begin of the time increment
    	const ConstitutiveStaticDataGradientDamage2DFatigue *oldStaticData = rElement->GetStaticData(rIp)->AsGradientDamage2DFatigue();

	    // extract total mechanical strain to begin of the time increment
    	const double prevNonlocalEqStrain(oldStaticData->mPrevNonlocalEqStrain.GetValue(0));
    	const EngineeringStress2D prevStress(oldStaticData->GetPrevStress());
	    const EngineeringStrain2D prevStrain(oldStaticData->GetPrevStrain());  // prevStrain is the mechanical strain, without thermal strain

	    // calculate damage and kappa
	    double e_0 = mTensileStrength / mE;
	    double kappa(std::max(oldStaticData->mKappa, e_0 + e_0*1.0e-5)), omega;

	    if (rNonlocalEqStrain - kappa >= 0.) {
	    	// calculate omega and kappa
	    	kappa = rNonlocalEqStrain;
	    	omega = CalculateDamage(kappa);

	    	// calculate nonlocal eq tangent
	    	if (rNewTangentNonlocal!=0) {
	    		(*rNewTangentNonlocal) = - CalculateDerivativeDamage(kappa) * Sgn(rNonlocalEqStrain - prevNonlocalEqStrain) * ElasticStiffness * rEngineeringStrain;
	    	}
	    } else {
	    	// call Newton and calculate omega and kappa

	    	// initialize the vector of unknowns
	    	Eigen::VectorXd Unknown(2);
	    	double deltaOmega;

	    	deltaOmega = CalculateDerivativeDamage(kappa) * pow(rNonlocalEqStrain/kappa,mViscosityExponent) * Plus(rNonlocalEqStrain - prevNonlocalEqStrain);

	    	Unknown[0] = deltaOmega;									// delta omega
	    	Unknown[1] = deltaOmega/CalculateDerivativeDamage(kappa); 	// delta kappa

	    	// compose vector of known Parameter, which are necessary for ResidualAn
	    	Eigen::VectorXd Parameter(4);

	    	Parameter[0] = oldStaticData->mOmega;									// omega at the beginning of the time increment
	    	Parameter[1] = std::max(oldStaticData->mKappa, e_0 + e_0*1.0e-5);		// kappa at the beginning of the time increment
	    	Parameter[2] = prevNonlocalEqStrain;									// nonlocal eq strain at the beginning of the time increment
	    	Parameter[3] = rNonlocalEqStrain;										// current nonlocal eq strain

	    	const Eigen::VectorXd ParameterList(Parameter);

	        // prepare starting Newton solver with respect to the "Unknown"
//	        NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> (GradientDamageEngineeringStressFatigue::*fdjacAn)
//	        	(const Eigen::VectorXd&,Eigen::VectorXd) const;
//
//	        Eigen::VectorXd (GradientDamageEngineeringStressFatigue::*residual)
//	        		(const Eigen::VectorXd&,Eigen::VectorXd) const;

	        // set Jacobi to analytical Jacobi
//	        fdjacAn = &GradientDamageEngineeringStressFatigue::DResidualAn;
//
//	        // set residual
//	        residual = &GradientDamageEngineeringStressFatigue::Residual;


	        boost::function<Eigen::VectorXd (const Eigen::VectorXd&,Eigen::VectorXd)> residualFunction;
	        residualFunction = boost::bind( &GradientDamageEngineeringStressFatigue::Residual, this, _1, _2 );

	        boost::function<NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>
	        	(const Eigen::VectorXd&,Eigen::VectorXd)> jacobiFunction;
	        jacobiFunction = boost::bind( &GradientDamageEngineeringStressFatigue::DResidualAn, this, _1, _2 );

	        // start Newton solver
	        try
	        {
#ifdef ENABLE_OPTIMIZE
	        	NuTo::NewtonRaphson myNonlinearSolver;

	        	// provide the solver with the pointers to the (!) residual, (!!) jacobi and (!!!) parameters necessary to evaluate the  both
	        	myNonlinearSolver.SetParameters(Parameter);
	        	myNonlinearSolver.SetResidualFunction(residualFunction);
	        	myNonlinearSolver.SetResidualDerivativeFunction(jacobiFunction);

	        	// solve
	        	myNonlinearSolver.Solve(Unknown);

//	        	// is the solution spurious?
//		        bool check;
//	        	check = myNonlinearSolver.GetCheckNewtonRaphson();
#else
    std::cout << "OPTIMIZE module not enabled - NewtonRaphson solver not available." << std::endl;
#endif // ENABLE_OPTIMIZE

	        }
	        catch (...)
	        {
	            rLogger << "[NuTo::GradientDamageEngineeringStressFatigue::ReturnMapping2D] No convergence in Newton." << "\n";
	        }

	        // provide omega and kappa
	        omega  = oldStaticData->mOmega + Unknown[0];
	        kappa += Unknown[1];

	    	// calculate nonlocal eq tangent
	        if (rNewTangentNonlocal!=0) {
	        	NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> matrixMultipl;

	        	matrixMultipl  = ((GradientDamageEngineeringStressFatigue::DResidualAn(ParameterList,Unknown)).fullPivLu().
	        			solve(GradientDamageEngineeringStressFatigue::DResidualDEpsNonlocalAn(ParameterList,Unknown))).eval();
	        	(*rNewTangentNonlocal) = matrixMultipl(0,0) * ElasticStiffness * rEngineeringStrain;
	        }
	    }

	    // calculate stress 2D
	    if (rNewStress2D!=0) {
	    	(*rNewStress2D) = ElasticStiffness * rEngineeringStrain;
	    }
	    // calculate stress 3D
	    if (rNewStress3D!=0) {
	    	(*rNewStress3D)[0] = (*rNewStress2D)[0];
	    	(*rNewStress3D)[1] = (*rNewStress2D)[1];
	    	(*rNewStress3D)[2] = 0.;
	    	(*rNewStress3D)[3] = (*rNewStress2D)[2];
	    	(*rNewStress3D)[4] = 0.;
	    	(*rNewStress3D)[5] = 0.;

	    	if (rElement->GetSection()->GetType() == Section::PLANE_STRAIN) {
	    		(*rNewStress3D)[2] = C12 * (rEngineeringStrain[0] + rEngineeringStrain[1]);
	    	}
	    }
	    // calculate local tangent
	    if (rNewTangent!=0) {
	    	(*rNewTangent) = (1. - omega) * ElasticStiffness;
	    }
	    // update statevs
	    if (rNewStaticData!=0) {
	    	rNewStaticData->mKappa = kappa;
	    	rNewStaticData->mOmega = omega;
	    	rNewStaticData->mPrevNonlocalEqStrain.SetValue(0,rNonlocalEqStrain);
	    	rNewStaticData->mPrevStrain = rEngineeringStrain;
	    	rNewStaticData->mPrevSigma = ElasticStiffness * rEngineeringStrain;
	    }

		return Error::SUCCESSFUL;
	}

//! @brief ... evaluate the constitutive relation in 3D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
//! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
NuTo::Error::eError NuTo::GradientDamageEngineeringStressFatigue::Evaluate3D(ElementBase* rElement, int rIp,
        const std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*>& rConstitutiveInput,
        std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>& rConstitutiveOutput)
{
	throw NuTo::MechanicsException("[NuTo::GradientDamageEngineeringStressFatigue::Evaluate3D] not implemented for 3D.");
//    // get section information determining which input on the constitutive level should be used
//    //const SectionBase* section(rElement->GetSection());
//
//    // check if parameters are valid
//    if (this->mParametersValid == false)
//    {
//        //throw an exception giving information related to the wrong parameter
//        CheckParameters();
//    }
//
//    // calculate engineering strain
//    EngineeringStrain3D strain3D;
//
//    if (rConstitutiveInput.find(NuTo::Constitutive::Input::DEFORMATION_GRADIENT_3D) == rConstitutiveInput.end())
//        throw MechanicsException("[NuTo::GradientDamageEngineeringStressFatigue::Evaluate2D] deformation gradient 3d needed to evaluate engineering strain3d.");
//    const DeformationGradient3D& deformationGradient(rConstitutiveInput.find(NuTo::Constitutive::Input::DEFORMATION_GRADIENT_3D)->second->GetDeformationGradient3D());
//    deformationGradient.GetEngineeringStrain(strain3D);
//
//    //Get previous ip_data
//    ConstitutiveStaticDataGradientDamage1DFatigue *oldStaticData = rElement->GetStaticData(rIp)->AsGradientDamage1DFatigue();
//
//    // nonlocal eq strain
//    if (rConstitutiveInput.find(NuTo::Constitutive::Input::NONLOCAL_EQ_STRAIN) == rConstitutiveInput.end())
//        throw MechanicsException("[NuTo::GradientDamageEngineeringStressFatigue::Evaluate3D] nonlocal eq strain needed to evaluate engineering stress 3d.");
//    const double nonlocalEqStrain = rConstitutiveInput.find(NuTo::Constitutive::Input::NONLOCAL_EQ_STRAIN)->second->GetNonlocalEqStrain().GetValue(0);
//
//    // last converged kappa
//    const double kappaLastConverged = oldStaticData->mKappa;
//
//    // kappa
//    const double kappa = std::max(nonlocalEqStrain, kappaLastConverged);
//
//    //calculate damage
//    double omega = CalculateDamage(kappa);
//
//    LocalEqStrain localEqStrain; // = CalculateLocalEqStrain3D(strain3D);
//    ConstitutiveTangentLocal<6, 1> localEqStrainTangent; // = CalculateLocalEqStrainTangent3D(strain3D);
//
//    CalculateLocalEqStrainAndDerivativeModifiedMises3D(strain3D, localEqStrain, localEqStrainTangent);
//
//    //set this to true, if update is in the map, perform the update after all other outputs have been calculated
//    bool performUpdateAtEnd(false);
//
//    for (std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>::iterator itOutput = rConstitutiveOutput.begin(); itOutput != rConstitutiveOutput.end(); itOutput++)
//    {
//        switch (itOutput->first)
//        {
//        case NuTo::Constitutive::Output::ENGINEERING_STRESS_3D:
//        {
//            EngineeringStress3D& engineeringStress3D = itOutput->second->GetEngineeringStress3D();
//
//            double C11, C12, C44;
//
//            this->CalculateCoefficients3D(C11, C12, C44);
//            // calculate Engineering stress
//            engineeringStress3D[0] = C11 * strain3D[0] + C12 * strain3D[1] + C12 * strain3D[2];
//            engineeringStress3D[1] = C11 * strain3D[1] + C12 * strain3D[0] + C12 * strain3D[2];
//            engineeringStress3D[2] = C11 * strain3D[2] + C12 * strain3D[0] + C12 * strain3D[1];
//            engineeringStress3D[3] = C44 * strain3D[3];
//            engineeringStress3D[4] = C44 * strain3D[4];
//            engineeringStress3D[5] = C44 * strain3D[5];
//
//            engineeringStress3D *= (1 - omega);
//        }
//            break;
//        case NuTo::Constitutive::Output::ENGINEERING_STRAIN_3D:
//        {
//            EngineeringStrain3D& engineeringStrain3D(itOutput->second->GetEngineeringStrain3D());
//            engineeringStrain3D[0] = strain3D[0];
//            engineeringStrain3D[1] = strain3D[1];
//            engineeringStrain3D[2] = strain3D[2];
//            engineeringStrain3D[3] = strain3D[3];
//            engineeringStrain3D[4] = strain3D[4];
//            engineeringStrain3D[5] = strain3D[5];
//        }
//            break;
//        case NuTo::Constitutive::Output::LOCAL_EQ_STRAIN:
//        {
//            itOutput->second->GetLocalEqStrain() = localEqStrain;
//        }
//            break;
//        case NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_3D:
//        {
//            ConstitutiveTangentLocal<6, 6>& tangent(itOutput->second->AsConstitutiveTangentLocal_6x6());
//            double C11, C12, C44;
//            // calculate coefficients of the material matrix
//
//            this->CalculateCoefficients3D(C11, C12, C44);
//
//            // store tangent at the output object
//            tangent(0, 0) = C11;
//            tangent(1, 0) = C12;
//            tangent(2, 0) = C12;
//            tangent(3, 0) = 0;
//            tangent(4, 0) = 0;
//            tangent(5, 0) = 0;
//
//            tangent(0, 1) = C12;
//            tangent(1, 1) = C11;
//            tangent(2, 1) = C12;
//            tangent(3, 1) = 0;
//            tangent(4, 1) = 0;
//            tangent(5, 1) = 0;
//
//            tangent(0, 2) = C12;
//            tangent(1, 2) = C12;
//            tangent(2, 2) = C11;
//            tangent(3, 2) = 0;
//            tangent(4, 2) = 0;
//            tangent(5, 2) = 0;
//
//            tangent(0, 3) = 0;
//            tangent(1, 3) = 0;
//            tangent(2, 3) = 0;
//            tangent(3, 3) = C44;
//            tangent(4, 3) = 0;
//            tangent(5, 3) = 0;
//
//            tangent(0, 4) = 0;
//            tangent(1, 4) = 0;
//            tangent(2, 4) = 0;
//            tangent(3, 4) = 0;
//            tangent(4, 4) = C44;
//            tangent(5, 4) = 0;
//
//            tangent(0, 5) = 0;
//            tangent(1, 5) = 0;
//            tangent(2, 5) = 0;
//            tangent(3, 5) = 0;
//            tangent(4, 5) = 0;
//            tangent(5, 5) = C44;
//
//            tangent *= (1 - omega);
//            tangent.SetSymmetry(true);
//        }
//            break;
//        case NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_NONLOCAL_EQ_STRAIN_3D:
//        {
//            ConstitutiveTangentLocal<6, 1>& tangent(itOutput->second->AsConstitutiveTangentLocal_6x1());
//            if (nonlocalEqStrain >= kappaLastConverged)
//            {
//                double damageDerivative = CalculateDerivativeDamage(nonlocalEqStrain);
//                // loading
//                double C11, C12, C44;
//                this->CalculateCoefficients3D(C11, C12, C44);
//                // calculate Engineering stress
//                tangent(0) = -damageDerivative * (C11 * strain3D[0] + C12 * strain3D[1] + C12 * strain3D[2] );
//                tangent(1) = -damageDerivative * (C11 * strain3D[1] + C12 * strain3D[0] + C12 * strain3D[2] );
//                tangent(2) = -damageDerivative * (C11 * strain3D[2] + C12 * strain3D[0] + C12 * strain3D[1] );
//
//                tangent(3) = -damageDerivative * (C44 * strain3D[3]);
//                tangent(4) = -damageDerivative * (C44 * strain3D[4]);
//                tangent(5) = -damageDerivative * (C44 * strain3D[5]);
//
//            } else
//            {
//                // unloading
//                tangent.setZero();
//            }
//        }
//            break;
//        case NuTo::Constitutive::Output::D_LOCAL_EQ_STRAIN_D_STRAIN_3D:
//        {
//            ConstitutiveTangentLocal<6, 1>& tangent(itOutput->second->AsConstitutiveTangentLocal_6x1());
//            tangent = localEqStrainTangent;
//        }
//            break;
////        case NuTo::Constitutive::Output::D_LOCAL_EQ_STRAIN_XI_D_STRAIN_3D:
////        {
////            ConstitutiveTangentLocal<6, 1>& tangent(itOutput->second->AsConstitutiveTangentLocal_6x1());
////
////            double xi = mNonlocalRadius;
////            double dXi = 0;
////
////            if (mNonlocalRadiusParameter != 0.)
////            {
////                throw MechanicsException("[NuTo::GradientDamageEngineeringStressFatigue::Evaluate3D] D_LOCAL_EQ_STRAIN_XI_D_STRAIN_3D not implemented for 3D for mNonlocalRadiusParameter != 0.");
////            }
////
////            double factor = 1. / xi - (localEqStrain[0] - nonlocalEqStrain) / (xi * xi) * dXi;
////
////            tangent = factor * localEqStrainTangent;
////
////        }
////            break;
////        case NuTo::Constitutive::Output::NONLOCAL_PARAMETER_XI:
////        {
////            ConstitutiveTangentLocal<1, 1>& variableNonlocalRadius(itOutput->second->AsConstitutiveTangentLocal_1x1());
////            variableNonlocalRadius[0] = mNonlocalRadius;
////            if (mNonlocalRadiusParameter != 0.)
////            {
////                throw MechanicsException("[NuTo::GradientDamageEngineeringStressFatigue::Evaluate3D] NONLOCAL_PARAMETER_XI not implemented for 3D for mNonlocalRadiusParameter != 0.");
////            }
////        }
////            break;
//        case NuTo::Constitutive::Output::DAMAGE:
//        {
//            itOutput->second->GetDamage().SetDamage(omega);
//        }
//            break;
//        case NuTo::Constitutive::Output::UPDATE_TMP_STATIC_DATA:
//        {
//            throw MechanicsException("[NuTo::GradientDamageEngineeringStressFatigue::Evaluate3D] tmp_static_data has to be updated without any other outputs, call it separately.");
//        }
//            break;
//        case NuTo::Constitutive::Output::UPDATE_STATIC_DATA:
//        {
//            performUpdateAtEnd = true;
//        }
//            break;
//        default:
//            throw MechanicsException(
//                    std::string("[NuTo::GradientDamageEngineeringStressFatigue::Evaluate3D] output object ") + NuTo::Constitutive::OutputToString(itOutput->first)
//                            + std::string(" could not be calculated, check the allocated material law and the section behavior."));
//        }
//    }
//
//    //update history variables
//    if (performUpdateAtEnd)
//    {
//        oldStaticData->mKappa = kappa;
//    }
//    return Error::SUCCESSFUL;
}

//! @brief ... create new static data object for an integration point
//! @return ... pointer to static data object
NuTo::ConstitutiveStaticDataBase* NuTo::GradientDamageEngineeringStressFatigue::AllocateStaticDataEngineeringStress_EngineeringStrain1D(const ElementBase* rElement) const
{
    return new ConstitutiveStaticDataGradientDamage1DFatigue();
}

//! @brief ... create new static data object for an integration point
//! @return ... pointer to static data object
NuTo::ConstitutiveStaticDataBase* NuTo::GradientDamageEngineeringStressFatigue::AllocateStaticDataEngineeringStress_EngineeringStrain2D(const ElementBase* rElement) const
{
    return new ConstitutiveStaticDataGradientDamage2DFatigue();
}

//! @brief ... create new static data object for an integration point
//! @return ... pointer to static data object
NuTo::ConstitutiveStaticDataBase* NuTo::GradientDamageEngineeringStressFatigue::AllocateStaticDataEngineeringStress_EngineeringStrain3D(const ElementBase* rElement) const
{
	throw MechanicsException(std::string("[NuTo::GradientDamageEngineeringStressFatigue::AllocateStaticDataEngineeringStress_EngineeringStrain3D] not implemented"));
}

// calculate coefficients of the material matrix
void NuTo::GradientDamageEngineeringStressFatigue::CalculateCoefficients3D(double& C11, double& C12, double& C44) const
{
    double factor = this->mE / ((1.0 + this->mNu) * (1.0 - 2.0 * this->mNu));
    C11 = factor * (1.0 - this->mNu);
    C12 = factor * this->mNu;
    C44 = this->mE / (2 * (1.0 + this->mNu));
}

// calculate coefficients of the material matrix
void NuTo::GradientDamageEngineeringStressFatigue::CalculateCoefficients2DPlaneStress(double& C11, double& C12, double& C33) const
{
    double factor = this->mE / (1.0 - (this->mNu * this->mNu));
    C11 = factor;
    C12 = factor * this->mNu;
    C33 = factor * 0.5 * (1.0 - this->mNu);
}
// parameters /////////////////////////////////////////////////////////////

//! @brief ... gets a parameter of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested parameter
//! @return ... value of the requested variable
double NuTo::GradientDamageEngineeringStressFatigue::GetParameterDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier) const
{
    switch(rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::COMPRESSIVE_STRENGTH:
        return mCompressiveStrength;
    case Constitutive::eConstitutiveParameter::DENSITY:
        return this->mRho;
    case Constitutive::eConstitutiveParameter::FRACTURE_ENERGY:
        return mFractureEnergy;
    case Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS:
        return this->mNonlocalRadius;
    case Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS_PARAMETER:
        return mNonlocalRadiusParameter;
    case Constitutive::eConstitutiveParameter::VISCOSITY_EXPONENT:
    	return mViscosityExponent;
    case Constitutive::eConstitutiveParameter::POISSONS_RATIO:
        return this->mNu;
    case Constitutive::eConstitutiveParameter::TENSILE_STRENGTH:
        return mTensileStrength;
    case Constitutive::eConstitutiveParameter::THERMAL_EXPANSION_COEFFICIENT:
        return this->mThermalExpansionCoefficient;
    case Constitutive::eConstitutiveParameter::YOUNGS_MODULUS:
        return this->mE;
    default:
    {
        throw MechanicsException("[NuTo::GradientDamageEngineeringStressFatigue::GetParameterDouble] Constitutive law does not have the requested variable");
    }
    }
}

//! @brief ... sets a parameter of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested parameter
//! @param rValue ... new value for requested variable
void NuTo::GradientDamageEngineeringStressFatigue::SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier, double rValue)
{
    switch(rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::COMPRESSIVE_STRENGTH:
    {
        this->CheckCompressiveStrength(rValue);
        this->mCompressiveStrength = rValue;
        this->SetParametersValid();
        break;
    }
    case Constitutive::eConstitutiveParameter::DENSITY:
    {
        this->CheckDensity(rValue);
        this->mRho = rValue;
        this->SetParametersValid();
        break;
    }
    case Constitutive::eConstitutiveParameter::FRACTURE_ENERGY:
    {
        this->CheckFractureEnergy(rValue);
        this->mFractureEnergy = rValue;
        this->SetParametersValid();
        break;
    }
    case Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS:
    {
        this->CheckNonlocalRadius(rValue);
        this->mNonlocalRadius = rValue;
        this->SetParametersValid();
        break;
    }
    case Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS_PARAMETER:
    {
        this->CheckNonlocalRadiusParameter(rValue);
        mNonlocalRadiusParameter = rValue;
        this->SetParametersValid();
        break;
    }
    case Constitutive::eConstitutiveParameter::VISCOSITY_EXPONENT:
    {
        this->CheckViscosityExponent(rValue);
        this->mViscosityExponent = rValue;
        this->SetParametersValid();
        break;
    }
    case Constitutive::eConstitutiveParameter::POISSONS_RATIO:
    {
        this->CheckPoissonsRatio(rValue);
        this->mNu = rValue;
        this->SetParametersValid();
        break;
    }
    case Constitutive::eConstitutiveParameter::TENSILE_STRENGTH:
    {
        this->CheckTensileStrength(rValue);
        this->mTensileStrength = rValue;
        this->SetParametersValid();
        break;
    }
    case Constitutive::eConstitutiveParameter::THERMAL_EXPANSION_COEFFICIENT:
    {
        this->CheckThermalExpansionCoefficient(rValue);
        this->mThermalExpansionCoefficient = rValue;
        this->SetParametersValid();
        break;
    }
    case Constitutive::eConstitutiveParameter::YOUNGS_MODULUS:
    {
        this->CheckYoungsModulus(rValue);
        this->mE = rValue;
        this->SetParametersValid();
        break;
    }
    default:
    {
        throw MechanicsException("[NuTo::GradientDamageEngineeringStressFatigue::SetParameterDouble] Constitutive law does not have the requested variable");
    }
    }
}


//! @brief ... gets a parameter of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested parameter
//! @return ... value of the requested variable
Eigen::VectorXd NuTo::GradientDamageEngineeringStressFatigue::GetParameterFullVectorDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier) const
{
    switch(rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::DAMAGE_LAW:
    {
        Eigen::VectorXd damageLaw(mDamageLawParameters.rows() + 1);

        damageLaw[0] = static_cast<double>(mDamageLawType);
        damageLaw.SetBlock(1, 0, mDamageLawParameters);
        return damageLaw;
        break;
    }
    default:
    {
        throw MechanicsException("[NuTo::GradientDamageEngineeringStressFatigue::GetParameterFullVectorDouble] Constitutive law does not have the requested variable");
    }
    }
}


//! @brief ... sets a parameter of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested parameter
//! @param rValue ... new value for requested variable
void NuTo::GradientDamageEngineeringStressFatigue::SetParameterFullVectorDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier, Eigen::VectorXd rValue)
{
    switch(rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::DAMAGE_LAW:
    {
        this->CheckDamageLaw(rValue);
        mDamageLawType = static_cast<int>(rValue[0]);
        int numDamageLawParameters = rValue.rows() - 1;
        if (numDamageLawParameters > 0)
            mDamageLawParameters = rValue.GetBlock(1, 0, numDamageLawParameters, 1);

        this->SetParametersValid();
        break;
    }
    default:
    {
        throw MechanicsException("[NuTo::GradientDamageEngineeringStressFatigue::SetParameterFullVectorDouble] Constitutive law does not have the requested variable");
    }
    }
}


///////////////////////////////////////////////////////////////////////////

//! @brief ... get type of constitutive relationship
//! @return ... type of constitutive relationship
//! @sa eConstitutiveType
NuTo::Constitutive::eConstitutiveType NuTo::GradientDamageEngineeringStressFatigue::GetType() const
{
    return NuTo::Constitutive::GRADIENT_DAMAGE_ENGINEERING_STRESS_FATIGUE;
}

//! @brief ... check compatibility between element type and type of constitutive relationship
//! @param rElementType ... element type
//! @return ... <B>true</B> if the element is compatible with the constitutive relationship, <B>false</B> otherwise.
bool NuTo::GradientDamageEngineeringStressFatigue::CheckElementCompatibility(NuTo::Element::eElementType rElementType) const
{
    switch (rElementType)
    {
    case NuTo::Element::ELEMENT1D:
        return true;
    case NuTo::Element::ELEMENT2D:
        return true;
    case NuTo::Element::ELEMENT3D:
//        return true;
    case NuTo::Element::BOUNDARYELEMENT1D:
        return true;
    case NuTo::Element::BOUNDARYELEMENT2D:
        return true;
//    case NuTo::Element::BOUNDARYELEMENT3D:
//        return true;
    default:
        return false;
    }
}

//! @brief ... check if density is non negative
//! @param rE ... Young's modulus
void NuTo::GradientDamageEngineeringStressFatigue::CheckDensity(double rRho) const
{
    if (rRho < 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::GradientDamageEngineeringStressFatigue::CheckDensity] The density must not be negative.");
    }
}

//! @brief ... check if Young's modulus is positive
//! @param rE ... Young's modulus
void NuTo::GradientDamageEngineeringStressFatigue::CheckYoungsModulus(double rE) const
{
    if (rE <= 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::GradientDamageEngineeringStressFatigue::CheckYoungsModulus] The Young's modulus must be a positive value.");
    }
}

//! @brief ... check if Poisson's ratio is valid \f$ (-1.0 < \nu < 0.5) \f$
//! @param rE ... Poisson's ratio
void NuTo::GradientDamageEngineeringStressFatigue::CheckPoissonsRatio(double rNu) const
{
    if (rNu <= -1.0)
    {
        throw NuTo::MechanicsException("[NuTo::GradientDamageEngineeringStressFatigue::CheckPoissonsRatio] Poisson's ratio must be greater or equal to -1.0.");
    }
    if (rNu >= 0.5)
    {
        throw NuTo::MechanicsException("[NuTo::GradientDamageEngineeringStressFatigue::CheckPoissonsRatio] Poisson's ratio must be smaller or equal to 0.5.");
    }
}

//! @brief ... check if the nonlocal radius is positive
//! @param rRadius ... nonlocal radius
void NuTo::GradientDamageEngineeringStressFatigue::CheckNonlocalRadius(double rRadius) const
{
    if (rRadius <= 0)
    {
        throw NuTo::MechanicsException("[NuTo::GradientDamageEngineeringStress::CheckNonlocalRadius] Nonlocal radius must be positive.");
    }
}

//! @brief ... check if the nonlocal radius parameter is greater than 1
//! @param rRadius ... nonlocal radius
void NuTo::GradientDamageEngineeringStressFatigue::CheckNonlocalRadiusParameter(double rRadiusParameter) const
{
    if (rRadiusParameter <= 1 and rRadiusParameter != 0.)
    {
        throw NuTo::MechanicsException("[NuTo::GradientDamageEngineeringStress::CheckNonlocalRadius] Nonlocal radius parameter must be greater than 1. or 0.");
    }
}

//! @brief ... check thermal expansion coefficient
//! @param rAlpha ... thermal expansion coefficient
void NuTo::GradientDamageEngineeringStressFatigue::CheckThermalExpansionCoefficient(double rAlpha) const
{
}

//! @brief ... check if tensile strength is positive
//! @param rTensileStrength ... nonlocal radius
void NuTo::GradientDamageEngineeringStressFatigue::CheckTensileStrength(double rTensileStrength) const
{
    if (rTensileStrength <= 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::GradientDamageEngineeringStressFatigue::CheckTensileStrength] The tensile strength must be a positive value.");
    }
}

//! @brief ... check if compressive strength is positive
//! @param rRadius ... compressive strength
void NuTo::GradientDamageEngineeringStressFatigue::CheckCompressiveStrength(double rCompressiveStrength) const
{
    if (rCompressiveStrength <= 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::GradientDamageEngineeringStressFatigue::CheckCompressiveStrength] The compressive strength must be a positive value.");
    }
}

//! @brief ... check if fracture energy is positive
//! @param rFractureEnergy ... fracture energy
void NuTo::GradientDamageEngineeringStressFatigue::CheckFractureEnergy(double rFractureEnergy) const
{
    if (rFractureEnergy <= 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::GradientDamageEngineeringStressFatigue::CheckFractureEnergy] The fracture energy must be a positive value.");
    }
}

//! @brief ... check whether the viscosity exponent is equal or larger than 1
//! @param rViscosity ... viscosity exponent
void NuTo::GradientDamageEngineeringStressFatigue::CheckViscosityExponent(double rViscosityExponent) const
{
    if (rViscosityExponent < 1.)
    {
        throw NuTo::MechanicsException("[NuTo::GradientDamageEngineeringStressFatigue::CheckViscosityExponent] The viscosity exponent must be equal or larger than 1.");
    }
}

//! @brief ... check damage law parameters
//! @param rDamageLawParameters ... damage law parameters
void NuTo::GradientDamageEngineeringStressFatigue::CheckDamageLaw(const Eigen::VectorXd& rDamageLaw) const
{
    int damageLawType = static_cast<int>(rDamageLaw[0]);
    int numDamageLawParameters = rDamageLaw.rows() - 1;

    switch (damageLawType)
    {
    case Constitutive::eDamageLawType::ISOTROPIC_NO_SOFTENING:
        break;
    case Constitutive::eDamageLawType::ISOTROPIC_LINEAR_SOFTENING:
        break;
    case Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING:
        break;
    case Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING_RES_LOAD:
        break;
    case Constitutive::eDamageLawType::ISOTROPIC_CUBIC_HERMITE:
        break;
    case Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING_SMOOTH:
    {
        if (numDamageLawParameters != 3)
            throw NuTo::MechanicsException(std::string("[NuTo::GradientDamageEngineeringStressFatigue::CheckDamageLaw] ") + std::string("ISOTROPIC_EXPONENTIAL_SOFTENING_SMOOTH needs exactly 3 parameters."));

    }
        break;
    default:
    {
        throw NuTo::MechanicsException(std::string("[NuTo::GradientDamageEngineeringStressFatigue::CheckDamageLaw] ") + std::string("The required damage law is not implemented. "));
        break;
    }

    }

}

//! @brief ... print information about the object
//! @param rVerboseLevel ... verbosity of the information
//! @param rLogger stream for the output
void NuTo::GradientDamageEngineeringStressFatigue::Info(unsigned short rVerboseLevel, Logger& rLogger) const
{
    this->ConstitutiveBase::Info(rVerboseLevel, rLogger);
    rLogger << "    Young's modulus      : " << this->mE << "\n";
    rLogger << "    Poisson's ratio      : " << this->mNu << "\n";
    rLogger << "    damage law params    : " << this->mDamageLawParameters << "\n";
    rLogger << "    thermal expansion coeff : " << this->mThermalExpansionCoefficient << "\n";
}

// check parameters
void NuTo::GradientDamageEngineeringStressFatigue::CheckParameters() const
{
    this->CheckDensity(mRho);
    this->CheckYoungsModulus(mE);
    this->CheckPoissonsRatio(mNu);
    this->CheckTensileStrength(mTensileStrength);
    this->CheckCompressiveStrength(mCompressiveStrength);
    this->CheckFractureEnergy(mFractureEnergy);
    this->CheckDamageLaw(GetParameterFullVectorDouble(Constitutive::eConstitutiveParameter::DAMAGE_LAW));
    this->CheckThermalExpansionCoefficient(mThermalExpansionCoefficient);
}

NuTo::LocalEqStrain NuTo::GradientDamageEngineeringStressFatigue::CalculateLocalEqStrain2D(const NuTo::EngineeringStrain2D& rStrain2D) const
{
    NuTo::LocalEqStrain localEqStrain;
    // calculate principal strains e1 and e2

    /**                          _____________________
     /    2              2
     exx   eyy      /  exy    /exx   eyy\
   e1 =   --- + --- +   /   ---- + |--- - ---|
     2     2    \/     4     \ 2     2 /

     A     +               B
     **/
    double A = (rStrain2D(0) + rStrain2D(1)) / 2.;
    double B = std::sqrt(std::pow(rStrain2D(0) - rStrain2D(1), 2.) / 4. + std::pow(rStrain2D(2), 2.) / 4.);

    // macaulay brackets of principal strains
    double e1 = std::max(A + B, 0.);
    double e2 = std::max(A - B, 0.);

    localEqStrain(0) = std::sqrt(e1 * e1 + e2 * e2);
    return localEqStrain;
}

NuTo::LocalEqStrain NuTo::GradientDamageEngineeringStressFatigue::CalculateLocalEqStrain3D(const NuTo::EngineeringStrain3D& rStrain3D) const
{
    throw NuTo::MechanicsException("[NuTo::GradientDamageEngineeringStress::CalculateLocalEqStrain3D] not implemented");
}

NuTo::ConstitutiveTangentLocal<3, 1> NuTo::GradientDamageEngineeringStressFatigue::CalculateLocalEqStrainTangent2D(const NuTo::EngineeringStrain2D& rStrain2D) const
{
    NuTo::ConstitutiveTangentLocal<3, 1> tangent;
    double A = (rStrain2D(0) + rStrain2D(1)) / 2.;
    double B = 0.5 * std::sqrt(std::pow(rStrain2D(0) - rStrain2D(1), 2.) + std::pow(rStrain2D(2), 2.));

    // macaulay brackets of principal strains
    double e1 = std::max(A + B, 0.);
    double e2 = std::max(A - B, 0.);

    if (e1 + e2 == 0)
    {
        tangent.setZero();
        return tangent;
    }
    if (e1 != 0 and e2 != 0)
    {
        double eMazar = std::sqrt(e1 * e1 + e2 * e2);
        tangent(0) = rStrain2D(0) / eMazar;
        tangent(1) = rStrain2D(1) / eMazar;
        tangent(2) = rStrain2D(2) / eMazar * 0.5;
        return tangent;
    }
    if (e1 != 0)
    {
        tangent(0) = (e1 - rStrain2D[1]) / (2 * B);
        tangent(1) = (e1 - rStrain2D[0]) / (2 * B);
        tangent(2) = 0.5 * rStrain2D[2] / (2 * B);
        return tangent;
    }
    if (e2 != 0)
    {
        tangent(0) = (rStrain2D[0] - e2) / (2 * B);
        tangent(1) = (rStrain2D[1] - e2) / (2 * B);
        tangent(2) = -0.5 * rStrain2D[2] / (2 * B);
        return tangent;
    }

//    std::cout << e1 << std::endl;
//    std::cout << e2 << std::endl;

    throw NuTo::MechanicsException("[NuTo::GradientDamageEngineeringStress::CalculateLocalEqStrainTangent2D] error");
}

NuTo::ConstitutiveTangentLocal<6, 1> NuTo::GradientDamageEngineeringStressFatigue::CalculateLocalEqStrainTangent3D(const NuTo::EngineeringStrain3D& rStrain3D) const
{
    throw NuTo::MechanicsException("[NuTo::GradientDamageEngineeringStressFatigue::CalculateLocalEqStrainTangent3D] not implemented");
}

//! @brief calculates the local eq strain and its derivatives for the modified mises model
//! @param rStrain2D ... 2d strain
//! @param rLocalEqStrain ... local eq strain
//! @param rLocalEqStrainTangent ... d local eq strain / d strain 2D
void NuTo::GradientDamageEngineeringStressFatigue::CalculateLocalEqStrainAndDerivativeModifiedMises2DPlaneStrain(const NuTo::EngineeringStrain2D& rStrain2D, NuTo::LocalEqStrain& rLocalEqStrain,
        NuTo::ConstitutiveTangentLocal<3, 1>& rLocalEqStrainTangent) const
{
    double k = mCompressiveStrength / mTensileStrength;

    double K1 = (k - 1.) / (2. * k * (1 - 2 * mNu));
    double K2 = 3 / (k * (1 + mNu) * (1 + mNu));

    double I1 = rStrain2D[0] + rStrain2D[1];
    double J2 = 1. / 3. * (rStrain2D[0] * rStrain2D[0] + rStrain2D[1] * rStrain2D[1] - rStrain2D[0] * rStrain2D[1]) + 0.25 * rStrain2D[2] * rStrain2D[2];

    double A = K1 * K1 * I1 * I1 + K2 * J2;
    if (A == 0) // prevent division by 0
    {
        rLocalEqStrain[0] = 0;

        rLocalEqStrainTangent[0] = K1;
        rLocalEqStrainTangent[1] = K1;
        rLocalEqStrainTangent[2] = 0;
        return;
    }

    rLocalEqStrain[0] = K1 * I1 + std::sqrt(A);

    double dJ2dexx = 1. / 3. * (2 * rStrain2D[0] - rStrain2D[1]);
    double dJ2deyy = 1. / 3. * (2 * rStrain2D[1] - rStrain2D[0]);
    double dJ2dgxy = 0.5 * rStrain2D[2];

    rLocalEqStrainTangent[0] = K1 + 1. / (2 * std::sqrt(A)) * (2 * K1 * K1 * I1 + K2 * dJ2dexx);
    rLocalEqStrainTangent[1] = K1 + 1. / (2 * std::sqrt(A)) * (2 * K1 * K1 * I1 + K2 * dJ2deyy);
    rLocalEqStrainTangent[2] = 1. / (2 * std::sqrt(A)) * (K2 * dJ2dgxy);

}

//! @brief calculates the local eq strain and its derivatives for the modified mises model
//! @param rStrain3D ... 3d strain
//! @param rLocalEqStrain ... local eq strain
//! @param rLocalEqStrainTangent ... d local eq strain / d strain 3D
void NuTo::GradientDamageEngineeringStressFatigue::CalculateLocalEqStrainAndDerivativeModifiedMises3D(const NuTo::EngineeringStrain3D& rStrain3D, NuTo::LocalEqStrain& rLocalEqStrain,
        NuTo::ConstitutiveTangentLocal<6, 1>& rLocalEqStrainTangent) const
{
    const double k = mCompressiveStrength / mTensileStrength;
    const double K1 = (k - 1.) / (2. * k * (1 - 2 * mNu));
    const double K2 = 3. / (k * (1 + mNu) * (1 + mNu));

    const double eps_xx = rStrain3D[0];
    const double eps_yy = rStrain3D[1];
    const double eps_zz = rStrain3D[2];
    const double eps_xy = 0.5*rStrain3D[3];
    const double eps_yz = 0.5*rStrain3D[4];
    const double eps_zx = 0.5*rStrain3D[5];

    double I1 = eps_xx + eps_yy + eps_zz;
    double J2 = 1. / 6. * (std::pow(eps_xx - eps_yy, 2) + std::pow(eps_yy - eps_zz, 2) + std::pow(eps_zz - eps_xx, 2))
            + std::pow(eps_xy, 2)+ std::pow(eps_yz, 2)+ std::pow(eps_zx, 2);

    double A = K1 * K1 * I1 * I1 + K2 * J2;

    if (A == 0) // prevent division by 0
    {
        rLocalEqStrain[0] = 0;

        rLocalEqStrainTangent[0] = (1 + mNu / (mNu - 1)) * K1;
        rLocalEqStrainTangent[1] = (1 + mNu / (mNu - 1)) * K1;
        rLocalEqStrainTangent[2] = (1 + mNu / (mNu - 1)) * K1;
        rLocalEqStrainTangent[3] = 0;
        rLocalEqStrainTangent[4] = 0;
        rLocalEqStrainTangent[5] = 0;

        return;
    }

    rLocalEqStrain[0] = K1 * I1 + std::sqrt(A);

    const double dJ2dexx = 1. / 3. * (2*eps_xx - eps_yy - eps_zz);
    const double dJ2deyy = 1. / 3. * (2*eps_yy - eps_xx - eps_zz);
    const double dJ2dezz = 1. / 3. * (2*eps_zz - eps_xx - eps_yy);
    const double dJ2dgxy = eps_xy;
    const double dJ2dgyz = eps_yz;
    const double dJ2dgzx = eps_zx;

    rLocalEqStrainTangent[0] = K1 + 1. / (2 * std::sqrt(A)) * (2 * K1 * K1 * I1 + K2 * dJ2dexx);
    rLocalEqStrainTangent[1] = K1 + 1. / (2 * std::sqrt(A)) * (2 * K1 * K1 * I1 + K2 * dJ2deyy);
    rLocalEqStrainTangent[2] = K1 + 1. / (2 * std::sqrt(A)) * (2 * K1 * K1 * I1 + K2 * dJ2dezz);
    rLocalEqStrainTangent[3] = 1. / (2 * std::sqrt(A)) * K2 * dJ2dgxy;
    rLocalEqStrainTangent[4] = 1. / (2 * std::sqrt(A)) * K2 * dJ2dgyz;
    rLocalEqStrainTangent[5] = 1. / (2 * std::sqrt(A)) * K2 * dJ2dgzx;

}

//! @brief calculates the local eq strain and its derivatives for the modified mises model in plane stress state
//! @param rStrain2D ... 2d strain
//! @param rLocalEqStrain ... local eq strain
//! @param rLocalEqStrainTangent ... d local eq strain / d strain 2D
void NuTo::GradientDamageEngineeringStressFatigue::CalculateLocalEqStrainAndDerivativeModifiedMises2DPlaneStress(const NuTo::EngineeringStrain2D& rStrain2D, NuTo::LocalEqStrain& rLocalEqStrain,
        NuTo::ConstitutiveTangentLocal<3, 1>& rLocalEqStrainTangent) const
{
    const double k = mCompressiveStrength / mTensileStrength;
    const double K1 = (k - 1.) / (2. * k * (1 - 2 * mNu));
    const double K2 = 3. / (k * (1 + mNu) * (1 + mNu));

    const double eps_xx = rStrain2D[0];
    const double eps_yy = rStrain2D[1];
    const double eps_zz = mNu / (mNu - 1) * (rStrain2D[0] + rStrain2D[1]);
    const double eps_xy = 0.5 * rStrain2D[2];

    double I1 = eps_xx + eps_yy + eps_zz;
    double J2 = 1. / 6. * (std::pow(eps_xx - eps_yy, 2) + std::pow(eps_yy - eps_zz, 2) + std::pow(eps_zz - eps_xx, 2)) + eps_xy * eps_xy;

    double A = K1 * K1 * I1 * I1 + K2 * J2;

    if (A == 0) // prevent division by 0
    {
        rLocalEqStrain[0] = 0;

        rLocalEqStrainTangent[0] = (1 + mNu / (mNu - 1)) * K1;
        rLocalEqStrainTangent[1] = (1 + mNu / (mNu - 1)) * K1;
        rLocalEqStrainTangent[2] = 0;

        return;
    }

    rLocalEqStrain[0] = K1 * I1 + std::sqrt(A);

    double dJ2dexx = 1. / 3. * (2. * eps_xx - eps_yy - 2. * eps_zz + 2 * mNu / (mNu - 1) * eps_zz);
    double dJ2deyy = 1. / 3. * (2. * eps_yy - eps_xx - 2. * eps_zz + 2 * mNu / (mNu - 1) * eps_zz);
    double dJ2dgxy = eps_xy;

    rLocalEqStrainTangent[0] = (1 + mNu / (mNu - 1)) * K1 + 1. / (2. * std::sqrt(A)) * (2. * (1 + mNu / (mNu - 1)) * K1 * K1 * I1 + K2 * dJ2dexx);
    rLocalEqStrainTangent[1] = (1 + mNu / (mNu - 1)) * K1 + 1. / (2. * std::sqrt(A)) * (2. * (1 + mNu / (mNu - 1)) * K1 * K1 * I1 + K2 * dJ2deyy);
    rLocalEqStrainTangent[2] = 1. / (2. * std::sqrt(A)) * (K2 * dJ2dgxy);

}
double NuTo::GradientDamageEngineeringStressFatigue::CalculateDamage(double rKappa) const
{
    double omega = 0;

    double e_0 = mTensileStrength / mE;
    double e_f = mFractureEnergy / mTensileStrength;
    double e_c = 2 * e_f; // or something

    switch (mDamageLawType)
    {
    case Constitutive::eDamageLawType::ISOTROPIC_NO_SOFTENING:
    {
        if (rKappa > e_0)
        {
            omega = 1 - e_0 / rKappa;
        }
        break;
    }
    case Constitutive::eDamageLawType::ISOTROPIC_LINEAR_SOFTENING:
    {
        if (rKappa > e_0)
        {
            omega = e_c / rKappa * (rKappa - e_0) / (e_c - e_0);
            omega = std::min(omega, MAX_OMEGA);
        }
        break;
    }
    case Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING:
    {
        if (rKappa > e_0)
        {
            omega = 1 - e_0 / rKappa * exp((e_0 - rKappa) / e_f);
        }
        break;
    }
    case Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING_RES_LOAD:
    {
        if (rKappa > e_0)
        {
            omega = 1 - e_0 / rKappa * (1 - MAX_OMEGA + MAX_OMEGA * exp((e_0 - rKappa) / e_f));
        }
        break;
    }
    case Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING_SMOOTH:
    {
        double alpha = mDamageLawParameters[0];
        double beta = mDamageLawParameters[1];
        double gamma = mDamageLawParameters[2];

        double e_Dev_e0 = rKappa / e_0;
        double e0_Dev_e = e_0 / rKappa;

        if (e_Dev_e0 <= alpha)
            break;

        if (e_Dev_e0 < beta)
        {
            // pre peak polynomial smoothing
            omega = 1. - e0_Dev_e * ((1 - alpha) * pow((e_Dev_e0 - beta) / (beta - alpha), 3.) + 1);
            break;
        }

        double term = 2. * e_f / (e_0 * (gamma - beta));
        double A = -1. / (1 + term);
        if (e_Dev_e0 < gamma)
        {
            // post peak polynomial smoothing
            omega = 1. - e0_Dev_e * (pow((e_Dev_e0 - beta) / (gamma - beta), 2.) * A + 1);
            break;
        }
        double e_s = e_0 * gamma + e_f * log(-A * term);

        // else
        omega = 1. - e0_Dev_e * exp((e_s - rKappa) / e_f);

        break;
    }
    case Constitutive::eDamageLawType::ISOTROPIC_CUBIC_HERMITE:
    {
        if (rKappa > e_0)
        {
            if (rKappa < e_c)
            {
                double kappa_scaled = (rKappa - e_0) / (e_c - e_0);
                omega = 1 - e_0 / rKappa * (2 * kappa_scaled * kappa_scaled * kappa_scaled - 3 * kappa_scaled * kappa_scaled + 1);
            } else
            {
                omega = 1.;
            }
        }
        break;
    }
    default:
        throw NuTo::MechanicsException(std::string("[NuTo::GradientDamageEngineeringStressFatigue::CalculateDamage] ") + std::string("The required damage law is not implemented. "));

        break;
    }

    return omega;
}

double NuTo::GradientDamageEngineeringStressFatigue::CalculateDerivativeDamage(double rKappa) const
{
    double DomegaDkappa = 0;

    double e_0 = mTensileStrength / mE;
    double e_f = mFractureEnergy / mTensileStrength;
    double e_c = 2 * e_f; // or something

    switch (mDamageLawType)
    {
    case Constitutive::eDamageLawType::ISOTROPIC_NO_SOFTENING:
    {
        if (rKappa > e_0)
        {
            DomegaDkappa = e_0 / (rKappa * rKappa);
        }
        break;
    }
    case Constitutive::eDamageLawType::ISOTROPIC_LINEAR_SOFTENING:
    {
        double termA = e_c / MAX_OMEGA / (e_c - e_0);
        double kappa_max = termA * e_0 / (termA - 1);

        if (rKappa > kappa_max)
            std::cout << CalculateDamage(kappa_max) << std::endl;

        if (rKappa > e_0 && rKappa < kappa_max)
        {
            DomegaDkappa = e_c * e_0 / (rKappa * rKappa * (e_c - e_0));
        }
        break;
    }
    case Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING:
    {
        if (rKappa > e_0)
        {
            DomegaDkappa = e_0 / rKappa * (1 / rKappa + 1 / e_f) * exp((e_0 - rKappa) / e_f);
        }
        break;
    }
    case Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING_RES_LOAD:
    {
        if (rKappa > e_0)
        {
            DomegaDkappa = e_0 / rKappa * ((1 / rKappa + 1 / e_f) * MAX_OMEGA * exp((e_0 - rKappa) / e_f) + (1 - MAX_OMEGA) / rKappa);
        }
        break;
    }
    case Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING_SMOOTH:
    {
        double alpha = mDamageLawParameters[0];
        double beta = mDamageLawParameters[1];
        double gamma = mDamageLawParameters[2];

        double e_Dev_e0 = rKappa / e_0;
        double e0_dev_e2 = e_0 / pow(rKappa, 2.);

        if (e_Dev_e0 <= alpha)
            break;

        if (e_Dev_e0 < beta)
        {
            // pre peak polynomial smoothing
            DomegaDkappa = e0_dev_e2 * (1. + (1. - alpha) / pow(beta - alpha, 3.) * pow(e_Dev_e0 - beta, 2.) * (-2. * e_Dev_e0 - beta));
            break;
        }

        double term = 2. * e_f / (e_0 * (gamma - beta));
        double A = -1. / (1 + term);
        if (e_Dev_e0 < gamma)
        {
            // post peak polynomial smoothing
            DomegaDkappa = e0_dev_e2 * (1. + (beta * beta - e_Dev_e0 * e_Dev_e0) / ((gamma - beta) * (gamma - beta - 2. * e_f / e_0)));
            break;
        }
        double e_s = e_0 * gamma + e_f * log(-A * term);

        // else
        DomegaDkappa = e0_dev_e2 * (rKappa / e_f + 1.) * exp((e_s - rKappa) / e_f);

        break;
    }
    case Constitutive::eDamageLawType::ISOTROPIC_CUBIC_HERMITE:
    {
        if (rKappa > e_0 && rKappa < e_c)
        {
            double kappa_scaled = (rKappa - e_0) / (e_c - e_0);
            DomegaDkappa = -6 * e_0 / rKappa / (e_c - e_0) * (kappa_scaled * kappa_scaled - kappa_scaled)
                    + e_0 / (rKappa * rKappa) * (2 * kappa_scaled * kappa_scaled * kappa_scaled - 3 * kappa_scaled * kappa_scaled + 1);

        }
        break;
    }
    default:
        throw NuTo::MechanicsException(std::string("[NuTo::GradientDamageEngineeringStressFatigue::CalculateDerivativeDamage] ") + std::string("The required damage law is not implemented. "));

        break;
    }

    return DomegaDkappa;
}

double NuTo::GradientDamageEngineeringStressFatigue::CalculateSecondDerivativeDamage(double rKappa) const
{
    double DDomegaDDkappa = 0;
    double DomegaDkappa = 0.;

    double e_0 = mTensileStrength / mE;
    double e_f = mFractureEnergy / mTensileStrength;
    double e_c = 2 * e_f; // or something

    switch (mDamageLawType)
    {
    case Constitutive::eDamageLawType::ISOTROPIC_NO_SOFTENING:
    {
        if (rKappa > e_0)
        {
            DDomegaDDkappa = -2. * e_0 / (rKappa * rKappa * rKappa);
        }
        break;
    }
    case Constitutive::eDamageLawType::ISOTROPIC_LINEAR_SOFTENING:
    {
        double termA = e_c / MAX_OMEGA / (e_c - e_0);
        double kappa_max = termA * e_0 / (termA - 1);

        if (rKappa > e_0 && rKappa < kappa_max)
        {
            DDomegaDDkappa = -2. * e_c * e_0 / (rKappa * rKappa * rKappa * (e_c - e_0));
        }
        break;
    }
    case Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING:
    {
        if (rKappa > e_0)
        {
            DomegaDkappa = e_0 / rKappa * (1 / rKappa + 1 / e_f) * exp((e_0 - rKappa) / e_f);
            DDomegaDDkappa = -(1 / rKappa + 1 / e_f) * DomegaDkappa;
            DDomegaDDkappa -= (e_0 / (rKappa * rKappa * rKappa)) * exp((e_0 - rKappa) / e_f);
        }
        break;
    }
    case Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING_RES_LOAD:
    {
        if (rKappa > e_0)
        {
            DomegaDkappa = MAX_OMEGA * e_0 / rKappa * (1 / rKappa + 1 / e_f) * exp((e_0 - rKappa) / e_f);
            DDomegaDDkappa = -(1 / rKappa + 1 / e_f) * DomegaDkappa;
            DDomegaDDkappa -= (MAX_OMEGA * e_0 / (rKappa * rKappa * rKappa)) * exp((e_0 - rKappa) / e_f);
            DDomegaDDkappa -= 2. * (1. - MAX_OMEGA) * e_0 /(rKappa * rKappa * rKappa);
        }
        break;
    }
    case Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING_SMOOTH:
    {
        throw NuTo::MechanicsException(std::string("[NuTo::GradientDamageEngineeringStressFatigue::CalculateSecondDerivativeDamage] ") + std::string("The second derivative is not implemented for the required damage law. "));

        break;
    }
    case Constitutive::eDamageLawType::ISOTROPIC_CUBIC_HERMITE:
    {
        if (rKappa > e_0 && rKappa < e_c)
        {
            double kappa_scaled = (rKappa - e_0) / (e_c - e_0);
            DDomegaDDkappa  = - (2. * kappa_scaled * kappa_scaled * kappa_scaled - 3. * kappa_scaled * kappa_scaled + 1.) / (rKappa * rKappa);
            DDomegaDDkappa += 6. * (kappa_scaled * kappa_scaled - kappa_scaled) / (rKappa * (e_c - e_0));
            DDomegaDDkappa -= 3. * (2. * kappa_scaled - 1.) / ((e_c - e_0) * (e_c - e_0));
            DDomegaDDkappa *= 2. * e_0 / rKappa;
        }
        break;
    }
    default:
        throw NuTo::MechanicsException(std::string("[NuTo::GradientDamageEngineeringStressFatigue::CalculateSecondDerivativeDamage] ") + std::string("The required damage law is not implemented. "));

        break;
    }

    return DomegaDkappa;
}

//! @brief ... returns true, if a material model has tmp static data (which has to be updated before stress or stiffness are calculated)
//! @return ... see brief explanation
bool NuTo::GradientDamageEngineeringStressFatigue::HaveTmpStaticData() const
{
    return false;
}
