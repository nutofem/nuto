#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain3D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress3D.h"
#include "nuto/mechanics/constitutive/mechanics/GradientDamagePlasticityEngineeringStress.h"
#include "nuto/mechanics/MechanicsException.h"

#include "nuto/math/SparseDirectSolverMUMPS.h"
#include "nuto/math/SparseMatrixCSRVector2General.h"

//#define PRINTRESULT
#include <eigen3/Eigen/Eigenvalues>
int main()
{
try
{
	bool testRoundedRankine3D(true);
	bool testDruckerPrager3D(true);


	if (testRoundedRankine3D)
    {
#ifdef PRINTRESULT
		std::cout << "test rounded rankine yield surface " << std::endl;
#endif //PRINTRESULT

		// test for the rounded rankine yield surface in 3D
		NuTo::GradientDamagePlasticityEngineeringStress myConstitutiveLaw;

		double fct(3.);
		Eigen::Matrix<double,10,6> stressCases;
		stressCases <<  4., -0.2, -0.2, 0., 0., 0.,        //D=0,q!=0 one positive
						4., 4., -0.2, 0., 0., 0.,      //D=0,q!=0 two positive
						4., 4., 5., 0., 0., 0.,    //D=0,q!=0 three positive
						4., -0.1, -0.2, 0.0, 0.0, 0.0,     //D!=0, one positive
						4., 0.2, -0.2, 0.01, 0.0, 0.0,     //D!=0, two positive
						4., 0.2,  0.2, 0.02, 0.0, 0.0,     //D!=0, three positive
						4., 4., 4., 0., 0., 0.,    //D=0,q=0 three positive
						4., 0.001, 0.001, 0., 0., 0.,      //D=0,q!0 three positive check the next three for the same gradient
						4., -0.001, -0.001, 0., 0., 0.,    //D=0,q=0 one positive
						4., -0.001, 0.001, 0., 0., 0.;     //D!=0   three positive
		Eigen::Matrix<double,6,1> stressVector;
		Eigen::Matrix<double,3,3> stressMatrix;
		Eigen::Matrix<double,3,3> stressMatrixRot;
		Eigen::Matrix<double,6,1> dF_dSigma_1, dF_dSigma_2, dF_dSigma_cdf;
		Eigen::Matrix<double,6,6> d2F_d2Sigma_1, d2F_d2Sigma_cdf;
		Eigen::Matrix<double,3,3> rotation1,rotation2,rotation3, rotation;
		Eigen::Matrix<double,1,3> angleCases;

		Eigen::Matrix<double,10,6> resultGradient;
		Eigen::Matrix<double,10,6> resultGradientRef;
		resultGradientRef << 0.592008 , 0.345492 ,  0.0625    ,  0.904508     ,-0.293893    , -0.38471,
				             0.639584 , 0.707107 ,  0.0675227 , -7.85414e-17  ,-1.13232e-19 , -0.415627,
				             0.542461 , 0.529813 ,  0.649618  , -5.10763e-17  , 2.38359e-17 ,  0.0778541,
				             0.592008 , 0.345492 ,  0.0625    ,  0.904508     ,-0.293893    , -0.38471,
				             0.604724 , 0.380117 ,  0.0638424 ,  0.859673     ,-0.279325    , -0.392973,
				             0.606577 , 0.382009 ,  0.108647  ,  0.860053     ,-0.279448    , -0.361768,
				             0.57735  , 0.57735  ,  0.57735   , -5.7347e-17   , 2.07797e-17 , -6.40988e-17,
				             0.59211  , 0.345655 ,  0.0627344 ,  0.904282     ,-0.293819    , -0.384614,
				             0.592008 , 0.345492 ,  0.0625    ,  0.904508     ,-0.293893    , -0.38471,
				             0.592032 , 0.345491 ,  0.0627261 ,  0.904508     ,-0.293893    , -0.384563;

		//angleCases << 0,0,0,
		//              M_PI*0.1,M_PI*0.2,M_PI*0.3;
		angleCases << M_PI*0.5,M_PI*.1,M_PI*0.3;
		double psi,theta,phi;
		double F_1, F_2;

		//check yield surface Rankine
		double delta=1e-8;
		for (int countAngle=0; countAngle<angleCases.rows(); countAngle++)
		{
			psi = angleCases(countAngle,0);
			theta = angleCases(countAngle,1);
			phi = angleCases(countAngle,2);

			rotation1 << cos(psi), sin(psi), 0, -sin(psi), cos(psi), 0, 0,0,1;
			rotation2 << 1,0,0,0, cos(theta), sin(theta), 0.,-sin(theta), cos(theta);
			rotation3 << cos(phi), sin(phi), 0, -sin(phi), cos(phi), 0, 0,0,1;

			rotation.noalias() = (rotation1*rotation2)*rotation3;

			for (int countStress=0; countStress<stressCases.rows(); countStress++)
			{
				stressMatrix << stressCases(countStress,0) , stressCases(countStress,3) , stressCases(countStress,4)
							  , stressCases(countStress,3) , stressCases(countStress,1) , stressCases(countStress,5)
							  , stressCases(countStress,4) , stressCases(countStress,5) , stressCases(countStress,2);
				stressMatrixRot = rotation * stressMatrix * rotation.transpose();

				//stressVector = stressCases.row(countStress).transpose();
				stressVector << stressMatrixRot(0,0) , stressMatrixRot(1,1) , stressMatrixRot(2,2),
								stressMatrixRot(0,1) , stressMatrixRot(1,2) , stressMatrixRot(0,2);

				F_1 = myConstitutiveLaw.YieldSurfaceRankine3DRounded(stressVector,fct,&dF_dSigma_1,&d2F_d2Sigma_1);
				for (int count=0; count<6; count++)
				{
					stressVector << stressMatrixRot(0,0) , stressMatrixRot(1,1) , stressMatrixRot(2,2),
									stressMatrixRot(0,1) , stressMatrixRot(1,2) , stressMatrixRot(0,2);
					stressVector(count,0)+=delta;
					F_2 = myConstitutiveLaw.YieldSurfaceRankine3DRounded(stressVector,fct,&dF_dSigma_2,0);

					dF_dSigma_cdf(count,0) = (F_2-F_1)/delta;
					d2F_d2Sigma_cdf.row(count) = (dF_dSigma_2-dF_dSigma_1)/delta;
				}
				double normDiffGrad((dF_dSigma_1-dF_dSigma_cdf).norm());
				double normDiffHess((d2F_d2Sigma_cdf-d2F_d2Sigma_1).norm());

#ifdef PRINTRESULT
				std::cout << "stress\n" << stressMatrixRot << "\n diff grad " << normDiffGrad << ", diff hess " << normDiffHess << std::endl;;
				std::cout << "norm Diff (grad, hess) " << normDiffGrad << "  " << normDiffHess << std::endl<< std::endl;

				std::cout << "dF_dSigma_1 \n  " << dF_dSigma_1.transpose() << std::endl;
#endif// PRINTRESULT
				resultGradient.row(countStress) = dF_dSigma_1.transpose();
				if (normDiffGrad > 1e-5 || (normDiffHess > 1e-5 && countStress<7))
				{
#ifdef PRINTRESULT
					std::cout << "dF_dSigma_cdf \n  " << dF_dSigma_cdf.transpose() << std::endl;

					std::cout << "d2F_d2Sigma_1 \n  " << d2F_d2Sigma_1 << std::endl;
					std::cout << "d2F_d2Sigma_cdf \n  " << d2F_d2Sigma_cdf << std::endl;
					throw NuTo::MechanicsException("Second derivative for rounded Rankine yield surface in GradientDamagePlasticityEngineeringStress is wrong.");
#endif// PRINTRESULT
				}
#ifdef PRINTRESULT
				std::cout << "======================================================================="<< std::endl;
#endif// PRINTRESULT
			}
#ifdef PRINTRESULT
			std::cout << "======================================================================="<< std::endl;
#endif// PRINTRESULT
		}
		if ((resultGradient-resultGradientRef).norm()>1e-5)
		{
#ifdef PRINTRESULT
			std::cout << "resultGradient" << resultGradient << std::endl << "resultGradientRef" << resultGradientRef << std::endl;

			std::cout << "norm Diff:" <<  (resultGradient-resultGradientRef).norm() << std::endl;
#endif// PRINTRESULT
			throw NuTo::MechanicsException("Result for rounded Rankine yield surface in GradientDamagePlasticityEngineeringStress is different than reference value.");
		}
    }

	if (testDruckerPrager3D)
    {
#ifdef PRINTRESULT
		std::cout << "test drucker prager yield surface " << std::endl;
#endif //PRINTRESULT

		// test for the rounded rankine yield surface in 3D
		NuTo::GradientDamagePlasticityEngineeringStress myConstitutiveLaw;

		double f_c1(30); //uniaxial compressive strength
		double f_c2(f_c1*1.12); //biaxial compressive strength
		double beta = sqrt(3)*(f_c2-f_c1) / (2.*f_c2-f_c1);
		double hp  = f_c2*f_c1 / (sqrt(3)*(2.*f_c2-f_c1));

		Eigen::Matrix<double,1,6> stressCases;
		stressCases <<  4., -0.2, -0.2, 0., 0., 0.;

		Eigen::Matrix<double,6,1> stressVector;
		Eigen::Matrix<double,3,3> stressMatrix;
		Eigen::Matrix<double,3,3> stressMatrixRot;
		Eigen::Matrix<double,6,1> dF_dSigma_1, dF_dSigma_2, dF_dSigma_cdf;
		Eigen::Matrix<double,6,6> d2F_d2Sigma_1, d2F_d2Sigma_cdf;
		Eigen::Matrix<double,3,3> rotation1,rotation2,rotation3, rotation;
		Eigen::Matrix<double,1,3> angleCases;

		Eigen::Matrix<double,1,6> resultGradient;
		Eigen::Matrix<double,1,6> resultGradientRef;
		resultGradientRef << 0.279892, 0.0664019, -0.178676,  0.783327, -0.254518, -0.333169;

		//angleCases << 0,0,0,
		//              M_PI*0.1,M_PI*0.2,M_PI*0.3;
		angleCases << M_PI*0.5,M_PI*.1,M_PI*0.3;
		double psi,theta,phi;
		double F_1, F_2;

		//check yield surface Drucker Prager
		double delta=1e-8;
		for (int countAngle=0; countAngle<angleCases.rows(); countAngle++)
		{
			psi = angleCases(countAngle,0);
			theta = angleCases(countAngle,1);
			phi = angleCases(countAngle,2);

			rotation1 << cos(psi), sin(psi), 0, -sin(psi), cos(psi), 0, 0,0,1;
			rotation2 << 1,0,0,0, cos(theta), sin(theta), 0.,-sin(theta), cos(theta);
			rotation3 << cos(phi), sin(phi), 0, -sin(phi), cos(phi), 0, 0,0,1;

			rotation.noalias() = (rotation1*rotation2)*rotation3;

			for (int countStress=0; countStress<stressCases.rows(); countStress++)
			{
				stressMatrix << stressCases(countStress,0) , stressCases(countStress,3) , stressCases(countStress,4)
							  , stressCases(countStress,3) , stressCases(countStress,1) , stressCases(countStress,5)
							  , stressCases(countStress,4) , stressCases(countStress,5) , stressCases(countStress,2);
				stressMatrixRot = rotation * stressMatrix * rotation.transpose();

				//stressVector = stressCases.row(countStress).transpose();
				stressVector << stressMatrixRot(0,0) , stressMatrixRot(1,1) , stressMatrixRot(2,2),
								stressMatrixRot(0,1) , stressMatrixRot(1,2) , stressMatrixRot(0,2);

				bool errorDerivatives(false);
				F_1 = myConstitutiveLaw.YieldSurfaceDruckerPrager3D(stressVector,beta,hp,&dF_dSigma_1,&d2F_d2Sigma_1,errorDerivatives);

				for (int count=0; count<6; count++)
				{
					stressVector << stressMatrixRot(0,0) , stressMatrixRot(1,1) , stressMatrixRot(2,2),
									stressMatrixRot(0,1) , stressMatrixRot(1,2) , stressMatrixRot(0,2);
					stressVector(count,0)+=delta;
					F_2 = myConstitutiveLaw.YieldSurfaceDruckerPrager3D(stressVector,beta,hp,&dF_dSigma_2,0, errorDerivatives);

					dF_dSigma_cdf(count,0) = (F_2-F_1)/delta;
					d2F_d2Sigma_cdf.row(count) = (dF_dSigma_2-dF_dSigma_1)/delta;
				}
				double normDiffGrad((dF_dSigma_1-dF_dSigma_cdf).norm());
				double normDiffHess((d2F_d2Sigma_cdf-d2F_d2Sigma_1).norm());

#ifdef PRINTRESULT
				std::cout << "stress\n" << stressMatrixRot << "\n diff grad " << normDiffGrad << ", diff hess " << normDiffHess << std::endl;;
				std::cout << "norm Diff (grad, hess) " << normDiffGrad << "  " << normDiffHess << std::endl<< std::endl;

				std::cout << "dF_dSigma_1 \n  " << dF_dSigma_1.transpose() << std::endl;
#endif// PRINTRESULT
				resultGradient.row(countStress) = dF_dSigma_1.transpose();
				if (normDiffGrad > 1e-5 || (normDiffHess > 1e-5 ))
				{
#ifdef PRINTRESULT
					std::cout << "dF_dSigma_cdf \n  " << dF_dSigma_cdf.transpose() << std::endl;

					std::cout << "d2F_d2Sigma_1 \n  " << d2F_d2Sigma_1 << std::endl;
					std::cout << "d2F_d2Sigma_cdf \n  " << d2F_d2Sigma_cdf << std::endl;
					throw NuTo::MechanicsException("Second derivative for Drucker Prager yield surface in GradientDamagePlasticityEngineeringStress is wrong.");
#endif// PRINTRESULT
				}
#ifdef PRINTRESULT
				std::cout << "======================================================================="<< std::endl;
#endif// PRINTRESULT
			}
#ifdef PRINTRESULT
			std::cout << "======================================================================="<< std::endl;
#endif// PRINTRESULT
		}
		if ((resultGradient-resultGradientRef).norm()>1e-5)
		{
#ifdef PRINTRESULT
			std::cout << "resultGradient" << resultGradient << std::endl << "resultGradientRef" << resultGradientRef << std::endl;

			std::cout << "norm Diff:" <<  (resultGradient-resultGradientRef).norm() << std::endl;
#endif// PRINTRESULT
			throw NuTo::MechanicsException("Result for rounded Rankine yield surface in GradientDamagePlasticityEngineeringStress is different than reference value.");
		}
    }

}
catch (NuTo::Exception& e)
{
    std::cout << e.ErrorMessage() << std::endl;
    return -1;
}
return 0;
}
