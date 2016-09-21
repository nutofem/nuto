#include <iostream>
#include "nuto/mechanics/constitutive/inputoutput/EngineeringStrain.h"
#include "nuto/mechanics/constitutive/inputoutput/EquivalentStrain.h"
#include "nuto/base/Timer.h"
#include <eigen3/Eigen/Dense>
#include <vector>

constexpr int numRuns = 10;

void CheckInvariants(NuTo::EngineeringStrain<3> r)
{

    Eigen::Matrix3d m;
    m <<     r[0], .5* r[5], .5* r[4],
         .5* r[5],     r[1], .5* r[3],
         .5* r[4], .5* r[3],     r[2];

    if (std::abs(r.InvariantI1() - m.trace()) > 1.e-10)
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "first strain invariant incorrect.");

    if (std::abs(r.InvariantI2() - .5 * (m.trace() * m.trace() - (m * m).trace())) > 1.e-10)
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "second strain invariant incorrect.");

    if (std::abs(r.InvariantI3() - m.determinant()) > 1.e-10)
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "third strain invariant incorrect.");

    // note: the invariants are not equal but with opposite sign
    if (std::abs(r.InvariantJ2() +  r.Deviatoric().InvariantI2()) > 1.e-10)
    {
        if (std::abs(r.InvariantI2() - r.Deviatoric().InvariantI2()) > 1.e-10)
            throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "second deviatoric strain defined with negative algebraic sign.");

        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "second deviatoric strain invariant incorrect.");
    }

    // note: the invariants are not equal but with opposite sign
    if (std::abs(r.InvariantJ2() -  (1./3. * r.InvariantI1() * r.InvariantI1() - r.InvariantI2())) > 1.e-10)
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "second deviatoric strain invariant incorrect.");
}

//! param planeState Defaults to plane stress; for 1D/3D irrelevant
template <int TDim>
void CheckLocalEqStrainDerivativesMises(NuTo::ePlaneState planeState = NuTo::ePlaneState::PLANE_STRESS)
{
    constexpr int VDim = NuTo::ConstitutiveIOBase::GetVoigtDim(TDim);


    double nu = 0.25;
    double k = 9.81;

    std::vector<Eigen::Matrix<double, VDim, 1>> strainCases(VDim*2, Eigen::Matrix<double, VDim, 1>::Zero());

    // define test cases with exactly one value =! 0
    for (int i = 0; i < VDim; ++i)
    {
        strainCases[i     ][i] =   M_PI;
        strainCases[i+VDim][i] = - M_PI;
    }

    // plus some random cases
    for (int i = 0; i < numRuns; ++i)
        strainCases.push_back(Eigen::Matrix<double, VDim, 1>::Random());

    double delta = 1.e-8;

    // check derivatives for plane stress and plane strain
    for (auto strainCase : strainCases)
    {
        NuTo::EngineeringStrain<TDim> strain;
        strain.AsVector() = strainCase;

        NuTo::EquivalentStrainModifiedMises<TDim> modMises(strain, k, nu, planeState);

        double localEqStrain0 = modMises.Get();

        NuTo::ConstitutiveVector<VDim> tangent = modMises.GetDerivative();
        NuTo::ConstitutiveVector<VDim> tangent_CDF;

        // calculate derivative numerically
        for (int i = 0; i < VDim; ++i)
        {
            strain[i] += delta;
            NuTo::EquivalentStrainModifiedMises<TDim> modMises1(strain, k, nu, planeState);
            tangent_CDF[i] = (modMises1.Get() - localEqStrain0) / delta;
            strain[i] -= delta;
        }

        if ((tangent - tangent_CDF).cwiseAbs().maxCoeff() > 1.e-6)
        {
            Eigen::IOFormat cleanFormat(7, 0, " ", "\n", "|", " |");
            std::cout << "strain:       " << strain.transpose().format(cleanFormat) << "\n";
            std::cout << "eq strain:    " << localEqStrain0 << "\n";
            std::cout << "tangent algo: " << tangent.transpose().format(cleanFormat) <<  "\n";
            std::cout << "tangent cdf : " << tangent_CDF.transpose().format(cleanFormat) << "\n";
            throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "wrong derivatives!");
        }
    }
}


int main()
{



    try
    {
        NuTo::Timer timer("Check invariants");

        for (int i = 0; i < numRuns; ++i)
        {
            NuTo::EngineeringStrain<1> e1D;
            NuTo::EngineeringStrain<2> e2D;
            NuTo::EngineeringStrain<3> e3D;

            e1D.AsVector() = Eigen::Matrix<double, 1, 1>::Random();
            e2D.AsVector() = Eigen::Matrix<double, 3, 1>::Random();
            e3D.AsVector() = Eigen::Matrix<double, 6, 1>::Random();

            CheckInvariants(e1D.As3D(0.3));
            CheckInvariants(e2D.As3D(0.3, NuTo::ePlaneState::PLANE_STRAIN));
            CheckInvariants(e2D.As3D(0.3, NuTo::ePlaneState::PLANE_STRESS));
            CheckInvariants(e3D);
        }
        timer.Reset("CheckLocalEqStrainDerivativesMises<1>");
        CheckLocalEqStrainDerivativesMises<1>();

        timer.Reset("CheckLocalEqStrainDerivativesMises<2> PLANE_STRAIN");
        CheckLocalEqStrainDerivativesMises<2>(NuTo::ePlaneState::PLANE_STRAIN);

        timer.Reset("CheckLocalEqStrainDerivativesMises<2> PLANE_STRESS");
        CheckLocalEqStrainDerivativesMises<2>(NuTo::ePlaneState::PLANE_STRESS);

        timer.Reset("CheckLocalEqStrainDerivativesMises<3>");
        CheckLocalEqStrainDerivativesMises<3>();

        timer.Reset("Finishing.");

    } catch (NuTo::MechanicsException& e)
    {
        std::cout << "Errors occurred. \n" << e.ErrorMessage() << std::endl;
        return EXIT_FAILURE;
    } catch (...)
    {
        std::cout << "Errors occurred. \n";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
