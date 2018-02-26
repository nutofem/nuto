#include <benchmark/benchmark.h>
#include "mechanics/nodes/NodeSimple.h"
#include "mechanics/integrationtypes/IntegrationTypeTensorProduct.h"
#include "mechanics/elements/ElementFem.h"
#include "mechanics/cell/Cell.h"
#include "mechanics/interpolation/InterpolationQuadSerendipity.h"
#include "mechanics/integrands/MomentumBalance.h"
#include "mechanics/constitutive/LinearElastic.h"

/*
 * Measures/Compares time for the calculation of a linear elastic gradient with quadratic quad elements
 *   - current NuTo implementation
 *   - (hopefully?) as fast as possible hardcode implementation
 */

template <int TDim>
class FixNode
{
public:
    FixNode(Eigen::Matrix<double, TDim, 1> rValues)
        : mValues(rValues)
    {
    }
    const Eigen::Matrix<double, TDim, 1>& Get() const
    {
        return mValues;
    }

private:
    Eigen::Matrix<double, TDim, 1> mValues;
};


class IntegrationTypeQuad
{
public:
    constexpr IntegrationTypeQuad(){};

    constexpr int GetNumIntegrationPoints() const
    {
        return 4;
    }
    Eigen::Vector2d GetLocalIntegrationPointCoordinates(int rIndex) const
    {
        constexpr double a = 0.57;
        return std::array<Eigen::Vector2d, 4>({Eigen::Vector2d({-a, -a}), Eigen::Vector2d({a, -a}),
                                               Eigen::Vector2d({-a, a}), Eigen::Vector2d({a, a})})[rIndex];
    }
};

class HardCodeElement8N
{

public:
    HardCodeElement8N(const std::vector<FixNode<2>*>& rNodes)
        : mNodes(rNodes)
    {
        constexpr IntegrationTypeQuad it;
        for (int i = 0; i < it.GetNumIntegrationPoints(); ++i)
        {
            const Eigen::Vector2d ip = it.GetLocalIntegrationPointCoordinates(i);
            mDerivativeShapeCache[i] = NuTo::ShapeFunctions2D::DerivativeShapeFunctionsQuadOrder2(ip);
        }
    }

    Eigen::Matrix<double, 16, 1> BuildInternalGradient() const
    {
        constexpr IntegrationTypeQuad it;
        Eigen::Matrix<double, 16, 1> result = Eigen::Matrix<double, 16, 1>::Zero();

        const auto disp = GetDisp();
        const auto coords = GetCoordinatesModified();

        for (int i = 0; i < it.GetNumIntegrationPoints(); ++i)
        {
            const auto& derivativeShapeFunctions = mDerivativeShapeCache[i];
            const auto J = GetJacobian(derivativeShapeFunctions, coords);
            const auto B = GetB(derivativeShapeFunctions, J);

            result += B.transpose() * mLaw.Stress(B * disp) * J.determinant();
        }
        return result;
    }

    Eigen::Matrix<double, 16, 1> GetDisp() const
    {
        Eigen::Matrix<double, 16, 1> disp;
        for (int i = 0; i < 8; ++i)
            disp.segment<2>(2 * i) = mNodes[i]->Get();
        return disp;
    }

    Eigen::Matrix<double, 2, 8> GetCoordinatesModified() const
    {
        Eigen::Matrix<double, 2, 8> coordinates;
        for (int i = 0; i < 8; ++i)
            coordinates.block<2, 1>(0, i) = mNodes[i]->Get();
        return coordinates;
    }

    Eigen::Matrix2d GetJacobian(const Eigen::Matrix<double, 8, 2>& rDerivativeShapeFunctions,
                                const Eigen::Matrix<double, 2, 8>& rCoords) const
    {
        return rCoords * rDerivativeShapeFunctions;
    }

    Eigen::Matrix<double, 3, 16> GetB(const Eigen::Matrix<double, 8, 2>& rDerivativeShapeFunctions,
                                      const Eigen::Matrix<double, 2, 2>& J) const
    {
        const Eigen::Matrix<double, 8, 2> derivativeShapeFunctionsJ = rDerivativeShapeFunctions * J.inverse();
        Eigen::Matrix<double, 3, 16> B = Eigen::Matrix<double, 3, 16>::Zero();

        for (int iNode = 0, iColumn = 0; iNode < 8; ++iNode, iColumn += 2)
        {
            const double dNdX = derivativeShapeFunctionsJ(iNode, 0);
            const double dNdY = derivativeShapeFunctionsJ(iNode, 1);

            B(0, iColumn) = dNdX;
            B(1, iColumn + 1) = dNdY;
            B(2, iColumn) = dNdY;
            B(2, iColumn + 1) = dNdX;
        }
        return B;
    }

private:
    std::vector<FixNode<2>*> mNodes;
    std::array<Eigen::Matrix<double, 8, 2>, 4> mDerivativeShapeCache;
    const NuTo::Laws::LinearElastic<2> mLaw = NuTo::Laws::LinearElastic<2>(20000, 0.3, NuTo::ePlaneState::PLANE_STRAIN);
};

static void Hardcode(benchmark::State& state)
{
    std::vector<FixNode<2>*> nodes;

    FixNode<2> n0(Eigen::Vector2d({0, 0}));
    FixNode<2> n1(Eigen::Vector2d({1, 0}));
    FixNode<2> n2(Eigen::Vector2d({1, 1}));
    FixNode<2> n3(Eigen::Vector2d({0, 1}));
    FixNode<2> n4(Eigen::Vector2d({0.5, 0}));
    FixNode<2> n5(Eigen::Vector2d({1, 0.5}));
    FixNode<2> n6(Eigen::Vector2d({0.5, 1}));
    FixNode<2> n7(Eigen::Vector2d({0, 0.5}));
    nodes.push_back(&n0);
    nodes.push_back(&n1);
    nodes.push_back(&n2);
    nodes.push_back(&n3);
    nodes.push_back(&n4);
    nodes.push_back(&n5);
    nodes.push_back(&n6);
    nodes.push_back(&n7);

    HardCodeElement8N e(nodes);

    for (auto _ : state)
        e.BuildInternalGradient();
}
BENCHMARK(Hardcode);

static void NuToPde(benchmark::State& state)
{
    NuTo::NodeSimple n0(Eigen::Vector2d({0, 0}));
    NuTo::NodeSimple n1(Eigen::Vector2d({1, 0}));
    NuTo::NodeSimple n2(Eigen::Vector2d({1, 1}));
    NuTo::NodeSimple n3(Eigen::Vector2d({0, 1}));
    NuTo::NodeSimple n4(Eigen::Vector2d({0.5, 0}));
    NuTo::NodeSimple n5(Eigen::Vector2d({1, 0.5}));
    NuTo::NodeSimple n6(Eigen::Vector2d({0.5, 1}));
    NuTo::NodeSimple n7(Eigen::Vector2d({0, 0.5}));
    std::vector<NuTo::NodeSimple*> coordNodes({&n0, &n1, &n2, &n3, &n4, &n5, &n6, &n7});
    NuTo::InterpolationQuadSerendipity coordInterpolation;
    NuTo::ElementFem coordElement(coordNodes, coordInterpolation);

    NuTo::NodeSimple nd0(Eigen::Vector2d({0, 0}));
    NuTo::NodeSimple nd1(Eigen::Vector2d({0, 0}));
    NuTo::NodeSimple nd2(Eigen::Vector2d({0, 0}));
    NuTo::NodeSimple nd3(Eigen::Vector2d({0, 0}));
    NuTo::NodeSimple nd4(Eigen::Vector2d({0, 0}));
    NuTo::NodeSimple nd5(Eigen::Vector2d({0, 0}));
    NuTo::NodeSimple nd6(Eigen::Vector2d({0, 0}));
    NuTo::NodeSimple nd7(Eigen::Vector2d({0, 0}));
    std::vector<NuTo::NodeSimple*> displNodes({&nd0, &nd1, &nd2, &nd3, &nd4, &nd5, &nd6, &nd7});
    NuTo::InterpolationQuadSerendipity displInterpolation;
    NuTo::ElementFem displElement(displNodes, displInterpolation);

    NuTo::DofType displDof("Displacements", 2);
    NuTo::ElementCollectionFem element(coordElement);
    element.AddDofElement(displDof, displElement);

    NuTo::Laws::LinearElastic<2> law(20000, 0.3, NuTo::ePlaneState::PLANE_STRAIN);
    NuTo::Integrands::MomentumBalance<2> integrand(displDof, law);
    const NuTo::IntegrationTypeTensorProduct<2> it(2, NuTo::eIntegrationMethod::GAUSS);

    NuTo::Cell cell(element, it, 0);

    auto Gradient = [&](const NuTo::CellIpData& cellIpData) { return integrand.Gradient(cellIpData, 0); };

    for (auto _ : state)
        cell.Integrate(Gradient);
}
BENCHMARK(NuToPde);
BENCHMARK_MAIN();
