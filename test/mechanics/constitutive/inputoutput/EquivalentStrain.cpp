#include "BoostUnitTest.h"

#include "mechanics/constitutive/inputoutput/EquivalentStrain.h"

using namespace NuTo;

//! param planeState Defaults to plane stress; for 1D/3D irrelevant
template <int TDim>
void CheckLocalEqStrainDerivativesMises(ePlaneState planeState = ePlaneState::PLANE_STRESS)
{
    constexpr int VDim = ConstitutiveIOBase::GetVoigtDim(TDim);

    double nu = 0.25;
    double k  = 9.81;

    std::vector<Eigen::Matrix<double, VDim, 1>> strainCases(VDim * 2, Eigen::Matrix<double, VDim, 1>::Zero());

    // define test cases with exactly one value =! 0
    for (int i = 0; i < VDim; ++i)
    {
        strainCases[i][i]        = M_PI;
        strainCases[i + VDim][i] = -M_PI;
    }

    // plus a random case
    strainCases.push_back(Eigen::Matrix<double, VDim, 1>::Random());

    double delta = 1.e-8;

    // check derivatives for plane stress and plane strain
    for (auto strainCase : strainCases)
    {
        EngineeringStrain<TDim> strain;
        strain.AsVector() = strainCase;

        EquivalentStrainModifiedMises<TDim> modMises(strain, k, nu, planeState);

        double localEqStrain0 = modMises.Get();

        ConstitutiveVector<VDim> tangent = modMises.GetDerivative();
        ConstitutiveVector<VDim> tangent_CDF;

        // calculate derivative numerically
        for (int i = 0; i < VDim; ++i)
        {
            strain[i] += delta;
            EquivalentStrainModifiedMises<TDim> modMises1(strain, k, nu, planeState);
            tangent_CDF[i] = (modMises1.Get() - localEqStrain0) / delta;
            strain[i] -= delta;
        }

        BOOST_CHECK_SMALL((tangent - tangent_CDF).cwiseAbs().maxCoeff(), 1.e-6);
    }
}


BOOST_AUTO_TEST_CASE(CheckEqStrainDerivatives)
{
    CheckLocalEqStrainDerivativesMises<1>();

    CheckLocalEqStrainDerivativesMises<2>(ePlaneState::PLANE_STRAIN);

    CheckLocalEqStrainDerivativesMises<2>(ePlaneState::PLANE_STRESS);

    CheckLocalEqStrainDerivativesMises<3>();
}


BOOST_AUTO_TEST_CASE(EqStrain)
{
    const double k  = 10;
    const double nu = 0.;
    EngineeringStrain<3> strain;
    strain[0] = 2.;
    BOOST_CHECK_CLOSE(EquivalentStrainModifiedMises<3>(strain, k, nu).Get(), 2., 1.e-6);
    strain[0] = 0.;
    strain[2] = - 2 * k;
    BOOST_CHECK_CLOSE(EquivalentStrainModifiedMises<3>(strain, k, nu).Get(), 2., 1.e-6);
}
