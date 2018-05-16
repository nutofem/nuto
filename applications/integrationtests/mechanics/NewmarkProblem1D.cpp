#include "BoostUnitTest.h"


#include "nuto/mechanics/cell/CellIpData.h"
#include "nuto/mechanics/dofs/DofMatrix.h"
#include "nuto/mechanics/dofs/DofType.h"
#include "nuto/mechanics/dofs/DofVector.h"
#include "nuto/mechanics/constitutive/EngineeringStress.h"
#include "nuto/mechanics/constitutive/EngineeringTangent.h"

#include <cmath>

namespace NuTo
{
namespace Integrands
{

//! @remark The area is fixed to 1m²
class MassDamperSpring1D
{
public:
    MassDamperSpring1D(DofType dofType, double e, double d, double rho)
        : mDofType(dofType)
        , mE{e}
        , mD{d}
        , mRho{rho}
    {
    }

    DofVector<double> Gradient(const CellIpData& cellIpData, double delta_t)
    {
        DofVector<double> gradient;

        Eigen::VectorXd dis = cellIpData.NodeValueVector(mDofType);
        Eigen::VectorXd vel = cellIpData.NodeValueVector(mDofType, 1);
        Eigen::VectorXd acc = cellIpData.NodeValueVector(mDofType, 2);

        gradient[mDofType] = Hessian0(cellIpData, delta_t)(mDofType, mDofType) * dis +
                             Hessian1(cellIpData, delta_t)(mDofType, mDofType) * vel +
                             Hessian2(cellIpData, delta_t)(mDofType, mDofType) * acc;

        return gradient;
    }

    DofMatrix<double> Hessian0(const CellIpData& cellIpData, double delta_t)
    {
        DofMatrix<double> hessian0;

        BMatrixStrain B = cellIpData.B(mDofType, Nabla::Strain());
        hessian0(mDofType, mDofType) = mE * B.transpose() * B;

        return hessian0;
    }

    DofMatrix<double> Hessian1(const CellIpData& cellIpData, double delta_t)
    {
        DofMatrix<double> hessian0;

        BMatrixStrain B = cellIpData.B(mDofType, Nabla::Strain());
        hessian0(mDofType, mDofType) = mD * B.transpose() * B;

        return hessian0;
    }

    DofMatrix<double> Hessian2(const CellIpData& cellIpData, double delta_t)
    {

        NMatrix N = cellIpData.N(mDofType);
        DofMatrix<double> hessian2;

        hessian2(mDofType, mDofType) = N.transpose() * N * mRho;
        return hessian2;
    }

protected:
    DofType mDofType;
    double mE = 0.;
    double mD = 0.;
    double mRho = 0.;
};
} /* Integrand */
} /* NuTo */

//! @brief Theoretical solution for a specific point in time
//! @param t : time
//! @param x_0 : initial displacement
//! @param e : stiffness
//! @param d : damping coefficient
//! @param m : mass
//! @remark source:
//! https://www.mrt.tu-berlin.de/fileadmin/fg20/Dokumente/Ueb-Skripte/MesstechnUebungen_I/MT_I_Zweimassenschwinger.pdf
double TheoreticalSolution(double t, double x_0, double e, double d, double m)
{
    double omega_0 = std::sqrt(e / m);
    double omega_d = omega_0 * std::sqrt(1 - d * d);

    return x_0 * std::exp(-omega_0 * d * t) * std::cos(omega_d * t);
}

BOOST_AUTO_TEST_CASE(MassSpringDamper)
{
    std::cout << TheoreticalSolution(0, 1., 15., 0.1, 1) << std::endl;
    std::cout << TheoreticalSolution(1, 1., 15., 0.1, 1) << std::endl;
    std::cout << TheoreticalSolution(1.5, 1., 15., 0.1, 1) << std::endl;
}
