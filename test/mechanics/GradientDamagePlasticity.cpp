#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain1D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain3D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress3D.h"
#include "nuto/mechanics/constitutive/mechanics/GradientDamagePlasticityEngineeringStress.h"
#include "nuto/mechanics/MechanicsException.h"

#include "nuto/math/SparseDirectSolverMUMPS.h"
#include "nuto/math/SparseMatrixCSRVector2General.h"
#include "nuto/base/Logger.h"

#define PRINTRESULT
#include <eigen3/Eigen/Eigenvalues>
int main()
{
try
{
    bool testRoundedRankineYieldSurface3D(true);
    bool testDruckerPragerYieldSurface3D(true);
    bool testRoundedRankineYieldSurface1D(true);
    bool testDruckerPragerYieldSurface1D(true);
    bool testReturnMapping3D(true);
    bool testReturnMapping1D(false); //I'm still working on that one

    NuTo::Logger logger;

    NuTo::GradientDamagePlasticityEngineeringStress myConstitutiveLaw;
    myConstitutiveLaw.SetDensity(1.);
    myConstitutiveLaw.SetYoungsModulus(40000.);
    myConstitutiveLaw.SetPoissonsRatio(0.2);
    myConstitutiveLaw.SetNonlocalRadius(2.);
    myConstitutiveLaw.SetTensileStrength(3.);
    myConstitutiveLaw.SetCompressiveStrength(30.);
    myConstitutiveLaw.SetBiaxialCompressiveStrength(1.12*30.);
    myConstitutiveLaw.SetFractureEnergy(0.1);
    myConstitutiveLaw.SetThermalExpansionCoefficient(0.);

    double f_c1(myConstitutiveLaw.GetCompressiveStrength()); //uniaxial compressive strength
    double f_c2(myConstitutiveLaw.GetBiaxialCompressiveStrength()); //biaxial compressive strength
    double beta = sqrt(3)*(f_c2-f_c1) / (2*f_c2-f_c1);
    double hp  = f_c2*f_c1 / (sqrt(3)*(2*f_c2-f_c1));
    double fct(myConstitutiveLaw.GetTensileStrength());

    // *********************************
    // rounded rankine yield surface 3D
    // *********************************
    if (testRoundedRankineYieldSurface3D)
    {
#ifdef PRINTRESULT
        std::cout << "test rounded rankine yield surface 3D" << std::endl;
#endif //PRINTRESULT


        Eigen::Matrix<double,6,1> stressVector;
        Eigen::Matrix<double,3,3> stressMatrix;
        Eigen::Matrix<double,3,3> stressMatrixRot;

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
                        4., -0.001, 0.001, 0., 0., 0.;    //D!=0   two positive

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

        //angleCases << 0.0,0.0,0.0;
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

                F_1 = myConstitutiveLaw.YieldSurfaceRoundedRankine3D(stressVector,fct,&dF_dSigma_1,&d2F_d2Sigma_1);
                for (int count=0; count<6; count++)
                {
                    stressVector << stressMatrixRot(0,0) , stressMatrixRot(1,1) , stressMatrixRot(2,2),
                                    stressMatrixRot(0,1) , stressMatrixRot(1,2) , stressMatrixRot(0,2);
                    stressVector(count,0)+=delta;
                    F_2 = myConstitutiveLaw.YieldSurfaceRoundedRankine3D(stressVector,fct,&dF_dSigma_2,0);

                    dF_dSigma_cdf(count,0) = (F_2-F_1)/delta;
                    d2F_d2Sigma_cdf.col(count) = (dF_dSigma_2-dF_dSigma_1)/delta;
                }
#ifdef PRINTRESULT
                double normDiffGrad((dF_dSigma_1-dF_dSigma_cdf).norm());
                double normDiffHess((d2F_d2Sigma_cdf-d2F_d2Sigma_1).norm());

                std::cout << countStress << ".stress\n" << stressMatrixRot << "\n diff grad " << normDiffGrad << ", diff hess " << normDiffHess << std::endl;;
                std::cout << "norm Diff (grad, hess) " << normDiffGrad << "  " << normDiffHess << std::endl<< std::endl;

                std::cout << "dF_dSigma_1 \n  " << dF_dSigma_1.transpose() << std::endl<< std::endl;
#endif// PRINTRESULT
                resultGradient.row(countStress) = dF_dSigma_1.transpose();
                //if (normDiffGrad > 1e-5 || (normDiffHess > 1e-5 && countStress<7))
                {
#ifdef PRINTRESULT
                    std::cout << "dF_dSigma_cdf \n  " << dF_dSigma_cdf.transpose() << std::endl<< std::endl;

                    std::cout << "d2F_d2Sigma_1 \n  " << d2F_d2Sigma_1 << std::endl<< std::endl;
                    std::cout << "d2F_d2Sigma_cdf \n  " << d2F_d2Sigma_cdf << std::endl<< std::endl;
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
            std::cout << "resultGradient\n" << resultGradient << std::endl << "resultGradientRef\n" << resultGradientRef << std::endl;

            std::cout << "norm Diff:" <<  (resultGradient-resultGradientRef).norm() << std::endl;
#endif// PRINTRESULT
            throw NuTo::MechanicsException("Result for rounded Rankine yield surface in GradientDamagePlasticityEngineeringStress is different than reference value.");
        }
    }

    // ********************************
    // Drucker Prager yield surface 3D
    // ********************************
    if (testDruckerPragerYieldSurface3D)
    {
#ifdef PRINTRESULT
        std::cout << "test Drucker Prager yield surface 3D" << std::endl;
#endif //PRINTRESULT

        Eigen::Matrix<double,1,6> stressCases;
        Eigen::Matrix<double,6,1> stressVector;
        Eigen::Matrix<double,3,3> stressMatrix;
        Eigen::Matrix<double,3,3> stressMatrixRot;
        stressCases <<  4., -0.2, -0.2, 0.1, 0.2, 0.3;        //
        Eigen::Matrix<double,6,1> dF_dSigma_1, dF_dSigma_2, dF_dSigma_cdf;
        Eigen::Matrix<double,6,6> d2F_d2Sigma_1, d2F_d2Sigma_cdf;
        Eigen::Matrix<double,3,3> rotation1,rotation2,rotation3, rotation;
        Eigen::Matrix<double,1,3> angleCases;

        Eigen::Matrix<double,1,6> resultGradient;
        Eigen::Matrix<double,1,6> resultGradientRef;
        resultGradientRef << 0.261482, 0.0856598, -0.179524,   0.74077, -0.395079, -0.313089;

        //angleCases << 0,0,0;
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

                bool(error);
                F_1 = myConstitutiveLaw.YieldSurfaceDruckerPrager3D(stressVector,beta,hp,&dF_dSigma_1,&d2F_d2Sigma_1, error);

                for (int count=0; count<6; count++)
                {
                    stressVector << stressMatrixRot(0,0) , stressMatrixRot(1,1) , stressMatrixRot(2,2),
                                    stressMatrixRot(0,1) , stressMatrixRot(1,2) , stressMatrixRot(0,2);
                    stressVector(count,0)+=delta;
                    F_2 = myConstitutiveLaw.YieldSurfaceDruckerPrager3D(stressVector,beta,hp,&dF_dSigma_2,0,error);

                    dF_dSigma_cdf(count,0) = (F_2-F_1)/delta;
                    d2F_d2Sigma_cdf.col(count) = (dF_dSigma_2-dF_dSigma_1)/delta;
                }
#ifdef PRINTRESULT
                double normDiffGrad((dF_dSigma_1-dF_dSigma_cdf).norm());
                double normDiffHess((d2F_d2Sigma_cdf-d2F_d2Sigma_1).norm());

                std::cout << countStress << ".stress\n" << stressMatrixRot << "\n diff grad " << normDiffGrad << ", diff hess " << normDiffHess << std::endl;;
                std::cout << "norm Diff (grad, hess) " << normDiffGrad << "  " << normDiffHess << std::endl<< std::endl;

                std::cout << "dF_dSigma_1 \n  " << dF_dSigma_1.transpose() << std::endl<< std::endl;
#endif// PRINTRESULT
                resultGradient.row(countStress) = dF_dSigma_1.transpose();
#ifdef PRINTRESULT
                if (normDiffGrad > 1e-5)
                {
                    std::cout << "dF_dSigma_cdf \n  " << dF_dSigma_cdf.transpose() << std::endl<< std::endl;

                    std::cout << "d2F_d2Sigma_1 \n  " << d2F_d2Sigma_1 << std::endl<< std::endl;
                    std::cout << "d2F_d2Sigma_cdf \n  " << d2F_d2Sigma_cdf << std::endl<< std::endl;
                }
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
            std::cout << "resultGradient\n" << resultGradient << std::endl << "resultGradientRef\n" << resultGradientRef << std::endl;

            std::cout << "norm Diff:" <<  (resultGradient-resultGradientRef).norm() << std::endl;
#endif// PRINTRESULT
            throw NuTo::MechanicsException("Result for rounded Rankine yield surface in GradientDamagePlasticityEngineeringStress is different than reference value.");
        }
    }

    // *********************************************************
    // Rounded Rankine yield surface 1D (uniaxial stress state)
    // *********************************************************
    if (testRoundedRankineYieldSurface1D)
    {
#ifdef PRINTRESULT
        std::cout << "test rounded Rankine yield surface 1D" << std::endl;
#endif //PRINTRESULT

        Eigen::Matrix<double,1,6> stressCases3D;
        stressCases3D <<  2*fct, 0.000001, 0.000001, 0., 0., 0.;        //
        Eigen::Matrix<double,1,1> stressVector1D;
        Eigen::Matrix<double,6,1> stressVector3D;

        Eigen::Matrix<double,1,1> dF_dSigma1_1D;
        Eigen::Matrix<double,5,1> dF_dSigma2_1D;
        Eigen::Matrix<double,6,1> dF_dSigma_1D;
        Eigen::Matrix<double,6,1> dF_dSigma_3D;

        Eigen::Matrix<double,6,6> d2F_d2Sigma_3D;
        Eigen::Matrix<double,1,1> d2F_d2Sigma1_1D;
        Eigen::Matrix<double,5,1> d2F_dSigma2dSigma1_1D;
        Eigen::Matrix<double,6,1> d2F_dSigma_1D;

        double deltaF,normDiffGrad,normDiffHess;
        //check yield surface Drucker Prager
		for (int countStress=0; countStress<stressCases3D.rows(); countStress++)
		{
			stressVector3D = stressCases3D.row(countStress).transpose();
			stressVector1D = stressVector3D.row(0);

			double f3D = myConstitutiveLaw.YieldSurfaceRoundedRankine3D(stressVector3D,fct,
					&dF_dSigma_3D,&d2F_d2Sigma_3D);
			double f1D = myConstitutiveLaw.YieldSurfaceRoundedRankine1D(stressVector1D,fct,
					&dF_dSigma1_1D, &dF_dSigma2_1D ,&d2F_d2Sigma1_1D, &d2F_dSigma2dSigma1_1D);

			deltaF = fabs(f3D-f1D);

			dF_dSigma_1D.block<1,1>(0,0) = dF_dSigma1_1D;
			dF_dSigma_1D.block<5,1>(1,0) = dF_dSigma2_1D;
			normDiffGrad = (dF_dSigma_1D-dF_dSigma_3D).norm();

			d2F_dSigma_1D.block<1,1>(0,0) = d2F_d2Sigma1_1D;
			d2F_dSigma_1D.block<5,1>(1,0) = d2F_dSigma2dSigma1_1D;
			normDiffHess = (d2F_dSigma_1D-d2F_d2Sigma_3D.col(0)).norm();

#ifdef PRINTRESULT
			std::cout << countStress << ".stress\n" << stressVector3D.transpose() << "\nf " << f1D << "\n diff yield function " << deltaF << " diff grad " << normDiffGrad << ", diff hess " << normDiffHess << std::endl;;

			std::cout << "dF_dSigma_1D \n  " << dF_dSigma_1D.transpose() << std::endl<< std::endl;
			std::cout << "d2F_dSigma_1D \n  " << d2F_dSigma_1D.transpose() << std::endl<< std::endl;
#endif// PRINTRESULT
			if (deltaF>1e-5 || normDiffGrad > 1e-5 || normDiffHess>1e-5)
			{
#ifdef PRINTRESULT
				std::cout << "dF_dSigma_3D \n  " << dF_dSigma_3D.transpose() << std::endl<< std::endl;
				std::cout << "d2F_d2Sigma_3D.col(0) \n  " << d2F_d2Sigma_3D.col(0).transpose() << std::endl<< std::endl;
#endif// PRINTRESULT
			}
#ifdef PRINTRESULT
            std::cout << "======================================================================="<< std::endl;
#endif// PRINTRESULT
        }
		if (deltaF>1e-5 || normDiffGrad > 1e-5 || normDiffHess>1e-5)
        {
            throw NuTo::MechanicsException("Result for rounded Rankine yield surface in 1D in GradientDamagePlasticityEngineeringStress in 1D and 3D.");
        }
    }

    // ********************************************************
    // Drucker Prager yield surface 1D (uniaxial stress state)
    // *******************************************************
    if (testDruckerPragerYieldSurface1D)
    {
#ifdef PRINTRESULT
        std::cout << "test Drucker Prager yield surface 1D" << std::endl;
#endif //PRINTRESULT

        Eigen::Matrix<double,1,6> stressCases3D;
        stressCases3D <<  -2.*f_c1, 0.000001, 0.000001, 0., 0., 0.;        //
        Eigen::Matrix<double,1,1> stressVector1D;
        Eigen::Matrix<double,6,1> stressVector3D;

        Eigen::Matrix<double,1,1> dF_dSigma1_1D;
        Eigen::Matrix<double,5,1> dF_dSigma2_1D;
        Eigen::Matrix<double,6,1> dF_dSigma_1D;
        Eigen::Matrix<double,6,1> dF_dSigma_3D;

        Eigen::Matrix<double,6,6> d2F_d2Sigma_3D;
        Eigen::Matrix<double,1,1> d2F_d2Sigma1_1D;
        Eigen::Matrix<double,5,1> d2F_dSigma2dSigma1_1D;
        Eigen::Matrix<double,6,1> d2F_dSigma_1D;

        double deltaF,normDiffGrad,normDiffHess;
        //check yield surface Drucker Prager
		for (int countStress=0; countStress<stressCases3D.rows(); countStress++)
		{
			stressVector3D = stressCases3D.row(countStress).transpose();
			stressVector1D = stressVector3D.row(0);

			bool(error);
			double f3D = myConstitutiveLaw.YieldSurfaceDruckerPrager3D(stressVector3D,beta,hp,
					&dF_dSigma_3D,&d2F_d2Sigma_3D, error);
			double f1D = myConstitutiveLaw.YieldSurfaceDruckerPrager1D(stressVector1D,beta,hp,
					&dF_dSigma1_1D, &dF_dSigma2_1D ,&d2F_d2Sigma1_1D, &d2F_dSigma2dSigma1_1D, error);

			deltaF = fabs(f3D-f1D);

			dF_dSigma_1D.block<1,1>(0,0) = dF_dSigma1_1D;
			dF_dSigma_1D.block<5,1>(1,0) = dF_dSigma2_1D;
			normDiffGrad = (dF_dSigma_1D-dF_dSigma_3D).norm();

			d2F_dSigma_1D.block<1,1>(0,0) = d2F_d2Sigma1_1D;
			d2F_dSigma_1D.block<5,1>(1,0) = d2F_dSigma2dSigma1_1D;
			normDiffHess = (d2F_dSigma_1D-d2F_d2Sigma_3D.col(0)).norm();

#ifdef PRINTRESULT
			std::cout << countStress << ".stress\n" << stressVector3D.transpose() << "\nf " << f1D << "\n diff yield function " << deltaF << " diff grad " << normDiffGrad << ", diff hess " << normDiffHess << std::endl;;

			std::cout << "dF_dSigma_1D \n  " << dF_dSigma_1D.transpose() << std::endl<< std::endl;
			std::cout << "d2F_dSigma_1D \n  " << d2F_dSigma_1D.transpose() << std::endl<< std::endl;
#endif// PRINTRESULT
			if (deltaF>1e-5 || normDiffGrad > 1e-5 || normDiffHess>1e-5)
			{
#ifdef PRINTRESULT
				std::cout << "dF_dSigma_3D \n  " << dF_dSigma_3D.transpose() << std::endl<< std::endl;
				std::cout << "d2F_d2Sigma_3D.col(0) \n  " << d2F_d2Sigma_3D.col(0).transpose() << std::endl<< std::endl;
#endif// PRINTRESULT
			}
#ifdef PRINTRESULT
            std::cout << "======================================================================="<< std::endl;
#endif// PRINTRESULT
        }
		if (deltaF>1e-5 || normDiffGrad > 1e-5 || normDiffHess>1e-5)
        {
            throw NuTo::MechanicsException("Result for Drucker-Prager yield surface in 1D in GradientDamagePlasticityEngineeringStress in 1D and 3D.");
        }
    }

    // ********************************
    // Return mapping 3D
    // ********************************
    if (testReturnMapping3D)
    {
#ifdef PRINTRESULT
        std::cout << "test return mapping in 3D" << std::endl;
#endif //PRINTRESULT

        Eigen::Matrix<double,6,18> strainCases; //total strains at t and t+1 and plastic strains at t
        strainCases <<   0    ,0,0,0,0,0 ,       0.0002, -0.000001, 0.000001, 0.      , 0.       , 0.          , 0      ,       0,       0,       0,       0,       0, //only Rankine uniaxial
                         0    ,0,0,0,0,0 ,       -0.002, 0.       , 0.      , 0.      , 0.       , 0.          , 0      ,       0,       0,       0,       0,       0, //only Drucker Prager uniaxial
                         0    ,0,0,0,0,0 ,       0.0002, 0.0001   , 0.00005 , 0.00005 , 0.00005  , 0.00005     , 0.00005, 0.00004, 0.00006, 0.00001, 0.00002, 0.00003, //only Rankine triaxial
                         0    ,0,0,0,0,0 ,       -0.002, 0.0001   , 0.00005 , 0.00005 , 0.00005  , 0.00005     , 0.00005, 0.00004, 0.00006, 0.00001, 0.00002, 0.00003, //only Drucker Prager triaxial
                         0.001,0,0,0,0,0 ,       -0.002, 0.002    , 0.0018  , 0.00005 , 0.00005  , 0.00005     , 0.00005, 0.00004, 0.00006, 0.00001, 0.00002, 0.00003,   //mixed triaxial
                         0.001,0,0,0,0,0 ,       -0.000, 0.000    , 0.000   , 0.00000 , 0.00000  , 0.00000     , 0.00105, 0.00004, 0.00006, 0.00001, 0.00002, 0.00003;   //unloading

        Eigen::Matrix<double,6,1> strainVector;
        Eigen::Matrix<double,3,3> strainMatrix;
        Eigen::Matrix<double,3,3> strainMatrixRot;
        Eigen::Matrix<double,6,1> sigma_1, sigma_2;
        Eigen::Matrix<double,6,6> d2sigma_d2epsilon_1, d2sigma_d2epsilon_cdf;
        Eigen::Matrix<double,3,3> rotation1,rotation2,rotation3, rotation;
        Eigen::Matrix<double,1,3> angleCases;

        Eigen::Matrix<double,6,6> resultStress;
        Eigen::Matrix<double,6,6> resultStressRef;
        resultStressRef <<
        		 1.86056,   1.18068,  0.421063,   1.24369,   -0.4041, -0.522927,
        		 -58.3571,  -48.2591,  -36.6669,  -18.5256,   6.01935,   7.87943,
        		    1.586,   1.82156,  0.367397,  0.914886, -0.593383, -0.571522,
        		 -59.6457,  -48.5145,  -37.8072,  -18.8036,   5.75406,    7.7786,
        		 -16.2175,  -9.17583,    1.0038,  -12.6593,   4.00921,   6.15978,
        		  -33.546,  -25.2426,  -16.5447,  -15.3908,   5.63212,   6.14616;

        angleCases << M_PI*0.5,M_PI*.1,M_PI*0.3;
        double psi,theta,phi;

        double delta=1e-9;
        for (int countAngle=0; countAngle<angleCases.rows(); countAngle++)
        {
            psi = angleCases(countAngle,0);
            theta = angleCases(countAngle,1);
            phi = angleCases(countAngle,2);

            rotation1 << cos(psi), sin(psi), 0, -sin(psi), cos(psi), 0, 0,0,1;
            rotation2 << 1,0,0,0, cos(theta), sin(theta), 0.,-sin(theta), cos(theta);
            rotation3 << cos(phi), sin(phi), 0, -sin(phi), cos(phi), 0, 0,0,1;

            rotation.noalias() = (rotation1*rotation2)*rotation3;

            for (int countStrain=0; countStrain<strainCases.rows(); countStrain++)
            {
            	//std::cout << "strainCase " << countStrain << std::endl;
                //prev total strain
                strainMatrix << strainCases(countStrain,0)     , 0.5*strainCases(countStrain,3) , 0.5*strainCases(countStrain,4)
                              , 0.5*strainCases(countStrain,3) , 0.5*strainCases(countStrain,1) , 0.5*strainCases(countStrain,5)
                              , 0.5*strainCases(countStrain,4) , 0.5*strainCases(countStrain,5) , strainCases(countStrain,2);
                strainMatrixRot = rotation * strainMatrix * rotation.transpose();

                double prevStrain3D[6];
                prevStrain3D[0] = strainMatrixRot(0,0);
                prevStrain3D[1] = strainMatrixRot(1,1);
                prevStrain3D[2] = strainMatrixRot(2,2);
                prevStrain3D[3] = 2.*strainMatrixRot(0,1);
                prevStrain3D[4] = 2.*strainMatrixRot(1,2);
                prevStrain3D[5] = 2.*strainMatrixRot(0,2);
                NuTo::EngineeringStrain3D prevEngineeringStrain3D;
                prevEngineeringStrain3D.SetData(prevStrain3D);

                //total strain at current time
                strainMatrix << strainCases(countStrain,0+6)     , 0.5*strainCases(countStrain,3+6) , 0.5*strainCases(countStrain,4+6)
                              , 0.5*strainCases(countStrain,3+6) , 0.5*strainCases(countStrain,1+6) , 0.5*strainCases(countStrain,5+6)
                              , 0.5*strainCases(countStrain,4+6) , 0.5*strainCases(countStrain,5+6) , strainCases(countStrain,2+6);
                strainMatrixRot = rotation * strainMatrix * rotation.transpose();

                double strain3D[6];
                strain3D[0] = strainMatrixRot(0,0);
                strain3D[1] = strainMatrixRot(1,1);
                strain3D[2] = strainMatrixRot(2,2);
                strain3D[3] = 2.*strainMatrixRot(0,1);
                strain3D[4] = 2.*strainMatrixRot(1,2);
                strain3D[5] = 2.*strainMatrixRot(0,2);
                NuTo::EngineeringStrain3D engineeringStrain3D;
                engineeringStrain3D.SetData(strain3D);

                //previous plastic strain
                strainMatrix << strainCases(countStrain,0+12)     , 0.5*strainCases(countStrain,3+12) , 0.5*strainCases(countStrain,4+12)
                              , 0.5*strainCases(countStrain,3+12) , 0.5*strainCases(countStrain,1+12) , 0.5*strainCases(countStrain,5+12)
                              , 0.5*strainCases(countStrain,4+12) , 0.5*strainCases(countStrain,5+12) , strainCases(countStrain,2+12);
                strainMatrixRot = rotation * strainMatrix * rotation.transpose();

                double prevPlasticStrain3D[6];
                prevPlasticStrain3D[0] = strainMatrixRot(0,0);
                prevPlasticStrain3D[1] = strainMatrixRot(1,1);
                prevPlasticStrain3D[2] = strainMatrixRot(2,2);
                prevPlasticStrain3D[3] = 2.*strainMatrixRot(0,1);
                prevPlasticStrain3D[4] = 2.*strainMatrixRot(1,2);
                prevPlasticStrain3D[5] = 2.*strainMatrixRot(0,2);

                double deltaEqPlasticStrain;
                Eigen::Matrix<double,6,6> dSigmadEpsilon,dSigmadEpsilonCDF;
                Eigen::Matrix<double,6,6> dEpsilonPdEpsilon,dEpsilonPdEpsilonCDF;
                Eigen::Matrix<double,6,1> newEpsilonP1,newEpsilonP2;
                Eigen::Matrix<double,6,1> newStress1,newStress2;

                myConstitutiveLaw.ReturnMapping3D(engineeringStrain3D, prevPlasticStrain3D, prevEngineeringStrain3D,
                        newStress1, newEpsilonP1,    deltaEqPlasticStrain, &dSigmadEpsilon, &dEpsilonPdEpsilon, logger);

                for (int count=0; count<6; count++)
                {
                	//std::cout << "checkComponent " << count << std::endl;
                    strainMatrix << strainCases(countStrain,0+6)     , 0.5*strainCases(countStrain,3+6) , 0.5*strainCases(countStrain,4+6)
                                  , 0.5*strainCases(countStrain,3+6) , 0.5*strainCases(countStrain,1+6) , 0.5*strainCases(countStrain,5+6)
                                  , 0.5*strainCases(countStrain,4+6) , 0.5*strainCases(countStrain,5+6) , strainCases(countStrain,2+6);

                    strainMatrixRot = rotation * strainMatrix * rotation.transpose();

                    strain3D[0] = strainMatrixRot(0,0);
                    strain3D[1] = strainMatrixRot(1,1);
                    strain3D[2] = strainMatrixRot(2,2);
                    strain3D[3] = 2.*strainMatrixRot(0,1);
                    strain3D[4] = 2.*strainMatrixRot(1,2);
                    strain3D[5] = 2.*strainMatrixRot(0,2);

                    strain3D[count]+=delta;

                    NuTo::EngineeringStrain3D engineeringStrain3D;
                    engineeringStrain3D.SetData(strain3D);

                    myConstitutiveLaw.ReturnMapping3D(engineeringStrain3D, prevPlasticStrain3D, prevEngineeringStrain3D,
                            newStress2, newEpsilonP2,    deltaEqPlasticStrain, 0, 0, logger);

                    dSigmadEpsilonCDF.col(count) = (newStress2-newStress1)/delta;
                    dEpsilonPdEpsilonCDF.col(count) = (newEpsilonP2-newEpsilonP1)/delta;
                }
                double normDiffdsde((dSigmadEpsilonCDF-dSigmadEpsilon).norm());
                double normDiffdepde((dEpsilonPdEpsilonCDF-dEpsilonPdEpsilon).norm());

#ifdef PRINTRESULT
                std::cout << countStrain<< ".strain \n\n" << strainMatrixRot << "\n\ndiff dsigma depsilon " << normDiffdsde << ", diff depsilonP depsilon " << normDiffdepde << std::endl << std::endl;;

                std::cout << "stress \n  " << newStress1.transpose() << std::endl<< std::endl;
#endif// PRINTRESULT
                resultStress.row(countStrain) = newStress1.transpose();
                if (normDiffdsde > 2e-1 || normDiffdepde > 1e-5 )
                {
#ifdef PRINTRESULT
                    std::cout << "dsigma depsilon exact\n  " << dSigmadEpsilon.transpose() << std::endl<< std::endl;
                    std::cout << "dsigma depsilon cdf  \n  " << dSigmadEpsilonCDF.transpose() << std::endl<< std::endl;

                    std::cout << "depsilonP depsilon exact\n  " << dEpsilonPdEpsilon.transpose() << std::endl<< std::endl;
                    std::cout << "depsilonP depsilon cdf  \n  " << dEpsilonPdEpsilonCDF.transpose() << std::endl<< std::endl;
#endif// PRINTRESULT
                    throw NuTo::MechanicsException("Derivatives for return mapping 3D in GradientDamagePlasticityEngineeringStress are different compared to central differences.");
                }
#ifdef PRINTRESULT
                std::cout << "======================================================================="<< std::endl;
#endif// PRINTRESULT
            }
#ifdef PRINTRESULT
            std::cout << "======================================================================="<< std::endl;
#endif// PRINTRESULT
        }
        if ((resultStress-resultStressRef).norm()>1e-2)
        {
#ifdef PRINTRESULT
            std::cout << "resultStress" << std::endl << resultStress << std::endl << "resultStressRef" << std::endl << resultStressRef << std::endl;

            std::cout << "norm Diff:" <<  (resultStress-resultStressRef).norm() << std::endl;
#endif// PRINTRESULT
            throw NuTo::MechanicsException("Result for return mapping 3D in GradientDamagePlasticityEngineeringStress is different than reference value.");
        }
    }

    // ********************************
    // Return mapping 1D
    // ********************************
    if (testReturnMapping1D)
    {
#ifdef PRINTRESULT
        std::cout << "test return mapping in 1D test" << std::endl;
#endif //PRINTRESULT

        Eigen::Matrix<double,1,3> strainCases; //total strains at t and t+1 and plastic strains at t (component in axial direction)
        strainCases <<   0. , 0.0004, 0.;

        Eigen::Matrix<double,1,1> sigma_1, sigma_2;
        Eigen::Matrix<double,1,1> d2sigma_d2epsilon_1, d2sigma_d2epsilon_cdf;

        Eigen::Matrix<double,1,1> resultStress;
        Eigen::Matrix<double,1,1> resultStressRef;
        resultStressRef <<   3.;


        double delta=1e-9;

		for (int countStrain=0; countStrain<strainCases.rows(); countStrain++)
		{
			//prev total strain
			NuTo::EngineeringStrain1D prevEngineeringStrain1D;
			prevEngineeringStrain1D[0] = strainCases(countStrain,0);

			//total strain at current time
			NuTo::EngineeringStrain1D engineeringStrain1D;
			engineeringStrain1D[0] = strainCases(countStrain,1);

			//previous plastic strain
			NuTo::EngineeringStrain1D prevPlasticStrain;
			prevPlasticStrain[0] = strainCases(countStrain,2);

			Eigen::Matrix<double,1,1> dSigmadEpsilon,dSigmadEpsilonCDF;
			Eigen::Matrix<double,2,1> dKappadEpsilon,dKappadEpsilonCDF;
			Eigen::Matrix<double,1,1> newEpsilonP1,newEpsilonP2;
			Eigen::Matrix<int,2,1> yieldConditionFlag;
			Eigen::Matrix<double,1,1> newStress1,newStress2;
			Eigen::Matrix<double,2,1> deltaEqPlasticStrain;

			myConstitutiveLaw.ReturnMapping1DNew(engineeringStrain1D, prevPlasticStrain, prevEngineeringStrain1D,
					newStress1, newEpsilonP1, yieldConditionFlag, deltaEqPlasticStrain, &dSigmadEpsilon, &dKappadEpsilon, logger);

			engineeringStrain1D[0] = strainCases(countStrain,2)+delta;

			myConstitutiveLaw.ReturnMapping1DNew(engineeringStrain1D, prevPlasticStrain, prevEngineeringStrain1D,
					newStress2, newEpsilonP2, yieldConditionFlag, deltaEqPlasticStrain, 0, 0, logger);

			dSigmadEpsilonCDF.col(0) = (newStress2-newStress1)/delta;
			//dEpsilonPdEpsilonCDF(0,0) = (newEpsilonP2(0,0)-newEpsilonP1(0,0))/delta;

#ifdef PRINTRESULT
			double normDiffdsde((dSigmadEpsilonCDF-dSigmadEpsilon).norm());
			//double normDiffdepde((dEpsilonPdEpsilonCDF-dEpsilonPdEpsilon).norm());

			std::cout << countStrain<< ".strains \n" << strainCases.row(countStrain) << "\n\ndiff dsigma depsilon "
					//<< normDiffdsde << ", diff depsilonP depsilon " << normDiffdepde
					  << std::endl << std::endl;
			std::cout << "stress \n  " << newStress1 << std::endl<< std::endl;
			std::cout << "new plastic strain \n  " << newEpsilonP1.transpose() << std::endl<< std::endl;
#endif// PRINTRESULT
			resultStress.row(countStrain) = newStress1.transpose();
			std::cout << "normDiffdsde \n  " << normDiffdsde << std::endl<< std::endl;
			//std::cout << "normDiffdepde \n  " << normDiffdepde << std::endl<< std::endl;
			if (normDiffdsde > 1e-1 )
			{
#ifdef PRINTRESULT
				std::cout << "dsigma depsilon exact\n  " << dSigmadEpsilon.transpose() << std::endl<< std::endl;
				std::cout << "dsigma depsilon cdf  \n  " << dSigmadEpsilonCDF.transpose() << std::endl<< std::endl;

				//std::cout << "depsilonP depsilon exact\n  " << dEpsilonPdEpsilon.transpose() << std::endl<< std::endl;
				//std::cout << "depsilonP depsilon cdf  \n  " << dEpsilonPdEpsilonCDF.transpose() << std::endl<< std::endl;
#endif// PRINTRESULT
				throw NuTo::MechanicsException("Derivatives for return mapping 3D in GradientDamagePlasticityEngineeringStress are different compared to central differences.");
			}
#ifdef PRINTRESULT
			std::cout << "======================================================================="<< std::endl;
#endif// PRINTRESULT
		}
        if ((resultStress-resultStressRef).norm()>1e-2)
        {
#ifdef PRINTRESULT
            std::cout << "resultStress" << std::endl << resultStress << std::endl << "resultStressRef" << std::endl << resultStressRef << std::endl;

            std::cout << "norm Diff:" <<  (resultStress-resultStressRef).norm() << std::endl;
#endif// PRINTRESULT
            throw NuTo::MechanicsException("Result for return mapping 3D in GradientDamagePlasticityEngineeringStress is different than reference value.");
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
