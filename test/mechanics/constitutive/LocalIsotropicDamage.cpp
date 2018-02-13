#include "BoostUnitTest.h"
#include "math/EigenCompanion.h"
#include "mechanics/constitutive/LocalIsotropicDamage.h"
#include "mechanics/constitutive/ModifiedMisesStrainNorm.h"
#include "mechanics/constitutive/LinearElastic.h"
#include "mechanics/constitutive/damageLaws/DamageLawExponential.h"

#include "math/NewtonRaphson.h"

using namespace NuTo;

// Material parameters
const double E = 20000;
const double nu = 0.2;
NuTo::Laws::LinearElastic<1> elasticLaw(E, nu);

const double ft = 4;
const double fc = 10 * ft;
const double gf = 0.01;
const double kappa0 = ft / E;
const double beta = ft / gf;
const double alpha = 1;

template <int TDim>
auto TestLaw()
{
    // Define the law using policy based design principles, hopefully applied correctly
    using Damage = Constitutive::DamageLawExponential;
    using StrainNorm = Constitutive::ModifiedMisesStrainNorm<TDim>;
    using Evolution = Laws::EvolutionImplicit<TDim, StrainNorm>;
    using ElasticLaw = Laws::LinearElastic<TDim>;
    using Law = Laws::LocalIsotropicDamage<TDim, Damage, Evolution, ElasticLaw>;

    Damage dmg(kappa0, beta, alpha);
    StrainNorm strainNorm(nu, fc / ft);
    Evolution evolutionEq(StrainNorm(nu, fc / ft), /*numCells=*/1, /*numIps=*/1);
    ElasticLaw linearElastic(E, nu);

    return Law(dmg, evolutionEq, linearElastic);
}

BOOST_AUTO_TEST_CASE(OneDimensional)
{
    auto localDamageLaw = TestLaw<1>();

    // Test the law by following the decreasing part of the load displacement curve after the peak load
    //
    /* Something like:
     *
     * stress
     *    ^
     *    |
     * ft |...
     *    |   /\
     *    |  /  `.
     *    | /     `•.
     *    |/          ` •  . _
     *    0----|---------------> strain
     *      kappa0
     */

    auto analyticStress = [&](double strain) {
        if (strain < kappa0)
            return E * strain;

        return ft * std::exp(ft / gf * (kappa0 - strain));
    };

    double stress = ft * 0.99;
    auto R = [&](double strain) { return localDamageLaw.Stress(EigenCompanion::ToEigen(strain), 0, 0, 0)[0] - stress; };
    auto DR = [&](double strain) { return localDamageLaw.Tangent(EigenCompanion::ToEigen(strain), 0, 0, 0)[0] - 1; };
    auto Norm = [&](double stress) { return std::abs(stress); };
    auto Info = [](int i, double x, double R) { BOOST_TEST_MESSAGE("" << i << ": x = " << x << " R = " << R); };
    auto problem = NuTo::NewtonRaphson::DefineProblem(R, DR, Norm, 1.e-14, Info);

    double strain = kappa0;
    localDamageLaw.Update(EigenCompanion::ToEigen(strain), 0, 0, 0);
    for (; stress > 0; stress -= ft * 0.1)
    {
        strain = NuTo::NewtonRaphson::Solve(problem, strain, NuTo::NewtonRaphson::DoubleSolver());
        localDamageLaw.Update(EigenCompanion::ToEigen(strain), 0, 0, 0);
        BOOST_CHECK_CLOSE(stress, analyticStress(strain), 1.e-10);
    }
}

void CheckTangent(std::initializer_list<double> values, double kappa)
{
    EngineeringStrain<3> strain = NuTo::EigenCompanion::ToEigen(values);
    auto law = TestLaw<3>();
    law.mEvolution.mKappas[0] = kappa;

    Eigen::Matrix<double, 6, 6> tangent = law.Tangent(strain, 0, 0, 0);
    Eigen::Matrix<double, 6, 6> tangent_cdf = law.Tangent(strain, 0, 0, 0) * 0.;

    const double delta = 1.e-8;

    for (int i = 0; i < 6; ++i)
    {
        strain[i] -= delta / 2.;
        auto s0 = law.Stress(strain, 0, 0, 0);
        strain[i] += delta;
        auto s1 = law.Stress(strain, 0, 0, 0);
        strain[i] -= delta / 2.;
        tangent_cdf.col(i) = (s1 - s0) / delta;
    }

    BOOST_CHECK_LT((tangent_cdf - tangent).cwiseAbs().maxCoeff(), 1.e-6);
}

BOOST_AUTO_TEST_CASE(Tangents)
{
    double kappa = kappa0 / 3.;
    CheckTangent({0., 0., 0., 0., 0., 0.}, kappa);
    CheckTangent({1.e-5, 0., 0., 0., 0., 0.}, kappa);
    CheckTangent({-1.e-5, 0., 0., 0., 0., 0.}, kappa);
    CheckTangent({1.e-5, 1.e-5, 0., 0., 0., 0.}, kappa);
    CheckTangent({2.e-5, 1.e-5, 0., 0., 0., 0.}, kappa);
    CheckTangent({2.e-5, -1.e-5, 0., 0., 0., 0.}, kappa);
    CheckTangent({0, 0, 2.e-5, 0., 0., 0.}, kappa);
    CheckTangent({1.e-5, 1.e-5, 2.e-5, 0., 0., 0.}, kappa);
    CheckTangent({1.e-5, -2.e-5, 2.e-5, 0., 0., 0.}, kappa);

    // some test in damaged loading
    kappa = 2 * kappa0;
    double eps = 1.e-6; // small load increment = damaged loading
    CheckTangent({kappa + eps, 0., 0., 0., 0., 0.}, kappa);
    CheckTangent({kappa, eps, 0., 0., 0., 0.}, kappa);
    CheckTangent({kappa, 0., eps, 0., 0., 0.}, kappa);

    CheckTangent({kappa + eps, +eps, 0., 0., 0., 0.}, kappa);
    CheckTangent({kappa, eps, eps, 0., 0., 0.}, kappa);
    CheckTangent({kappa + eps, 0., eps, 0., 0., 0.}, kappa);


    // decrement = elastic unloading
    CheckTangent({kappa - eps, 0., 0., 0., 0., 0.}, kappa);
    CheckTangent({kappa, -eps, 0., 0., 0., 0.}, kappa);
    CheckTangent({kappa, 0., -eps, 0., 0., 0.}, kappa);

    CheckTangent({kappa - eps, -eps, 0., 0., 0., 0.}, kappa);
    CheckTangent({kappa, -eps, -eps, 0., 0., 0.}, kappa);
    CheckTangent({kappa - eps, 0., -eps, 0., 0., 0.}, kappa);
}

BOOST_AUTO_TEST_CASE(EvolutionEdgeCase)
{
    auto law = TestLaw<3>();
    EngineeringStrain<3> strain = EigenCompanion::ToEigen({1, 2, 3, 4, 5, 6}) * 1.e-5;
    double kappa = law.mEvolution.mStrainNorm.Value(strain);
    law.mEvolution.mKappas[0] = kappa;

    BOOST_CHECK_GT(law.mEvolution.DkappaDstrain(strain, 0, 0, 0)(0, 0), 0.);
}
