#include "Benchmark.h"
#include <eigen3/Eigen/Core>
#include "mechanics/elements/ElementOutputBlockVectorDouble.h"
#include "mechanics/dofSubMatrixStorage/BlockFullVector.h"
#include "mechanics/MechanicsEnums.h"
#include "mechanics/elements/ElementShapeFunctions.h"
#include "mechanics/constitutive/laws/EngineeringStressHelper.h"
#include "mechanics/nodes/NodeDof.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/dofSubMatrixStorage/BlockFullMatrix.h"
#include "mechanics/integrationtypes/IntegrationTypeTensorProduct.h"
#include "mechanics/elements/ElementBase.h"
#include "mechanics/cell/Cell.h"
#include "mechanics/interpolation/InterpolationQuadSerendipity.h"
#include "mechanics/cell/IntegrandLinearElastic.h"
#include "mechanics/sections/SectionPlane.h"

namespace Benchmark
{

class HardCodeElement8NDynamic
{

public:
    HardCodeElement8NDynamic(const std::vector<NuTo::NodeBase*>& rNodes)
        : mNodes(rNodes)
    {
    }

    Eigen::VectorXd BuildInternalGradient() const
    {
        const NuTo::IntegrationTypeTensorProduct<2> it(2, NuTo::eIntegrationMethod::GAUSS);

        Eigen::VectorXd result = Eigen::VectorXd::Zero(16, 1);

        const auto disp = GetDisp();
        const auto coords = GetCoordinatesModified();

        for (int i = 0; i < it.GetNumIntegrationPoints(); ++i)
        {
            const Eigen::Vector2d ip = it.GetLocalIntegrationPointCoordinates(i);

            const auto derivativeShapeFunctions = NuTo::ShapeFunctions2D::DerivativeShapeFunctionsQuadOrder2(ip);
            const auto J = GetJacobian(derivativeShapeFunctions, coords);

            const auto B = GetB(derivativeShapeFunctions, J);

            result += B.transpose() * mLaw.Stress(B * disp) * J.determinant();
        }
        return result;
    }

    Eigen::VectorXd GetDisp() const
    {
        Eigen::VectorXd disp(16);
        for (int i = 0; i < 8; ++i)
        {
            const auto& nodeValues = mNodes[i]->Get(NuTo::Node::eDof::DISPLACEMENTS);
            disp(2 * i) = nodeValues[0];
            disp(2 * i + 1) = nodeValues[1];
        }
        return disp;
    }

    Eigen::MatrixXd GetCoordinatesModified() const
    {

        Eigen::Matrix<double, 2, Eigen::Dynamic> coordinates(2, 8);
        for (int i = 0; i < 8; ++i)
        {
            const auto& nodeCoordinate = mNodes[i]->Get(NuTo::Node::eDof::COORDINATES);
            coordinates(0, i) = nodeCoordinate[0];
            coordinates(1, i) = nodeCoordinate[1];
        }
        return coordinates;
    }

    Eigen::Matrix2d GetJacobian(const Eigen::MatrixXd& rDerivativeShapeFunctions, const Eigen::MatrixXd& rCoords) const
    {
        return rCoords * rDerivativeShapeFunctions;
    }

    Eigen::MatrixXd GetB(const Eigen::MatrixXd& rDerivativeShapeFunctions, const Eigen::Matrix2d& J) const
    {
        const Eigen::MatrixXd derivativeShapeFunctionsJ = rDerivativeShapeFunctions * J.inverse();
        Eigen::MatrixXd B = Eigen::MatrixXd::Zero(3, 16);

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
    std::vector<NuTo::NodeBase*> mNodes;
    const NuTo::LinearElasticLaw2D mLaw = NuTo::LinearElasticLaw2D(20000, 0.3);
};

template <int TDim>
class Node
{
public:
    Node(Eigen::Matrix<double, TDim, 1> rValues)
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

class HardCodeElement8NDynamicStaticAlloc
{

    template <int TMaxRows>
    using SemiFixVector = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor, TMaxRows, 1>;
    template <int TMaxRows, int TMaxCols>
    using SemiFixMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor, TMaxRows, TMaxCols>;

public:
    HardCodeElement8NDynamicStaticAlloc(const std::vector<NuTo::NodeBase*>& rNodes)
        : mNodes(rNodes)
    {
    }

    SemiFixVector<NuTo::maxDim * NuTo::maxNumNodes> BuildInternalGradient() const
    {
        const NuTo::IntegrationTypeTensorProduct<2> it(2, NuTo::eIntegrationMethod::GAUSS);
        SemiFixVector<NuTo::maxDim* NuTo::maxNumNodes> result =
                SemiFixVector<NuTo::maxDim * NuTo::maxNumNodes>::Zero(16, 1);
        const auto disp = GetDisp();
        const auto coords = GetCoordinatesModified();


        for (int i = 0; i < it.GetNumIntegrationPoints(); ++i)
        {
            const Eigen::Vector2d ip = it.GetLocalIntegrationPointCoordinates(i);

            const auto derivativeShapeFunctions = NuTo::ShapeFunctions2D::DerivativeShapeFunctionsQuadOrder2(ip);
            const auto J = GetJacobian(derivativeShapeFunctions, coords);

            const auto B = GetB(derivativeShapeFunctions, J);

            result += B.transpose() * mLaw.Stress(B * disp) * J.determinant();
        }
        return result;
    }

    SemiFixVector<NuTo::maxDim * NuTo::maxNumNodes> GetDisp() const
    {
        SemiFixVector<NuTo::maxDim * NuTo::maxNumNodes> disp;
        for (int i = 0; i < 8; ++i)
        {
            const auto& nodeValues = mNodes[i]->Get(NuTo::Node::eDof::DISPLACEMENTS);
            disp(2 * i) = nodeValues[0];
            disp(2 * i + 1) = nodeValues[1];
        }
        return disp;
    }


    SemiFixMatrix<NuTo::maxDim, NuTo::maxNumNodes> GetCoordinatesModified() const
    {
        SemiFixMatrix<NuTo::maxDim, NuTo::maxNumNodes> coordinates;
        for (int i = 0; i < 8; ++i)
        {
            const auto& nodeCoordinate = mNodes[i]->Get(NuTo::Node::eDof::COORDINATES);
            coordinates(0, i) = nodeCoordinate[0];
            coordinates(1, i) = nodeCoordinate[1];
        }
        return coordinates;
    }

    Eigen::Matrix2d GetJacobian(const SemiFixMatrix<NuTo::maxNumNodes, NuTo::maxDim>& rDerivativeShapeFunctions,
                                const SemiFixMatrix<NuTo::maxDim, NuTo::maxNumNodes>& rCoords) const
    {
        return rCoords * rDerivativeShapeFunctions;
    }

    SemiFixMatrix<6, NuTo::maxDim * NuTo::maxNumNodes>
    GetB(const SemiFixMatrix<NuTo::maxNumNodes, NuTo::maxDim>& rDerivativeShapeFunctions,
         const Eigen::Matrix2d& J) const
    {
        const SemiFixMatrix<6, NuTo::maxDim* NuTo::maxNumNodes> derivativeShapeFunctionsJ =
                rDerivativeShapeFunctions * J.inverse();
        SemiFixMatrix<6, NuTo::maxDim* NuTo::maxNumNodes> B =
                SemiFixMatrix<6, NuTo::maxDim * NuTo::maxNumNodes>::Zero(3, 16);

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
    std::vector<NuTo::NodeBase*> mNodes;
    const NuTo::LinearElasticLaw2D mLaw = NuTo::LinearElasticLaw2D(20000, 0.3);
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
    HardCodeElement8N(const std::vector<Node<2>*>& rNodes)
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
    std::vector<Node<2>*> mNodes;
    std::array<Eigen::Matrix<double, 8, 2>, 4> mDerivativeShapeCache;
    const NuTo::LinearElasticLaw2D mLaw = NuTo::LinearElasticLaw2D(20000, 0.3);
};

template <int TDim>
class SemiHardCodeElement8N
{

public:
    SemiHardCodeElement8N(const std::vector<Node<TDim>*>& rNodes)
        : mNodes(rNodes)
    {
        constexpr IntegrationTypeQuad it;
        for (int i = 0; i < it.GetNumIntegrationPoints(); ++i)
        {
            const auto ip = it.GetLocalIntegrationPointCoordinates(i);
            mDerivativeShapeCache[i] = NuTo::ShapeFunctions2D::DerivativeShapeFunctionsQuadOrder2(ip);
        }
    }

    Eigen::VectorXd BuildInternalGradient() const
    {
        constexpr IntegrationTypeQuad it;
        Eigen::VectorXd result = Eigen::VectorXd::Zero(mNodes.size() * TDim);
        const Eigen::VectorXd disp = GetDisp();
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

    Eigen::VectorXd GetDisp() const
    {
        Eigen::VectorXd disp(mNodes.size() * TDim);
        for (unsigned i = 0; i < mNodes.size(); ++i)
            disp.segment<2>(2 * i) = mNodes[i]->Get();
        return disp;
    }

    Eigen::Matrix<double, TDim, Eigen::Dynamic> GetCoordinatesModified() const
    {
        Eigen::Matrix<double, TDim, Eigen::Dynamic> coordinates(TDim, mNodes.size());
        for (unsigned i = 0; i < mNodes.size(); ++i)
            coordinates.template block<TDim, 1>(0, i) = mNodes[i]->Get();
        return coordinates;
    }

    Eigen::Matrix<double, TDim, TDim>
    GetJacobian(const Eigen::Matrix<double, Eigen::Dynamic, TDim>& rDerivativeShapeFunctions,
                const Eigen::Matrix<double, TDim, Eigen::Dynamic>& rCoords) const
    {
        return rCoords * rDerivativeShapeFunctions;
    }

    Eigen::Matrix<double, 3, Eigen::Dynamic>
    GetB(const Eigen::Matrix<double, Eigen::Dynamic, TDim>& rDerivativeShapeFunctions,
         const Eigen::Matrix<double, TDim, TDim>& J) const
    {
        const Eigen::Matrix<double, Eigen::Dynamic, TDim> derivativeShapeFunctionsJ =
                rDerivativeShapeFunctions * J.inverse();
        Eigen::Matrix<double, 3, Eigen::Dynamic> B =
                Eigen::Matrix<double, 3, Eigen::Dynamic>::Zero(3, mNodes.size() * TDim);

        for (unsigned iNode = 0, iColumn = 0; iNode < mNodes.size(); ++iNode, iColumn += 2)
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
    std::vector<Node<2>*> mNodes;
    std::array<Eigen::Matrix<double, Eigen::Dynamic, TDim>, 4> mDerivativeShapeCache;
    const NuTo::LinearElasticLaw2D mLaw = NuTo::LinearElasticLaw2D(20000, 0.3);
};
}

BENCHMARK(BuildGradient, HardcodeFixEverything, runner)
{
    std::vector<Benchmark::Node<2>*> nodes;

    Benchmark::Node<2> n0(Eigen::Vector2d({0, 0}));
    Benchmark::Node<2> n1(Eigen::Vector2d({1, 0}));
    Benchmark::Node<2> n2(Eigen::Vector2d({1, 1}));
    Benchmark::Node<2> n3(Eigen::Vector2d({0, 1}));
    Benchmark::Node<2> n4(Eigen::Vector2d({0.5, 0}));
    Benchmark::Node<2> n5(Eigen::Vector2d({1, 0.5}));
    Benchmark::Node<2> n6(Eigen::Vector2d({0.5, 1}));
    Benchmark::Node<2> n7(Eigen::Vector2d({0, 0.5}));
    nodes.push_back(&n0);
    nodes.push_back(&n1);
    nodes.push_back(&n2);
    nodes.push_back(&n3);
    nodes.push_back(&n4);
    nodes.push_back(&n5);
    nodes.push_back(&n6);
    nodes.push_back(&n7);

    Benchmark::HardCodeElement8N e(nodes);

    while (runner.KeepRunningIterations(1e6))
    {
        e.BuildInternalGradient();
    }
}

BENCHMARK(BuildGradient, HardcodeFixNodeDimension, runner)
{
    std::vector<Benchmark::Node<2>*> nodes;

    Benchmark::Node<2> n0(Eigen::Vector2d({0, 0}));
    Benchmark::Node<2> n1(Eigen::Vector2d({1, 0}));
    Benchmark::Node<2> n2(Eigen::Vector2d({1, 1}));
    Benchmark::Node<2> n3(Eigen::Vector2d({0, 1}));
    Benchmark::Node<2> n4(Eigen::Vector2d({0.5, 0}));
    Benchmark::Node<2> n5(Eigen::Vector2d({1, 0.5}));
    Benchmark::Node<2> n6(Eigen::Vector2d({0.5, 1}));
    Benchmark::Node<2> n7(Eigen::Vector2d({0, 0.5}));
    nodes.push_back(&n0);
    nodes.push_back(&n1);
    nodes.push_back(&n2);
    nodes.push_back(&n3);
    nodes.push_back(&n4);
    nodes.push_back(&n5);
    nodes.push_back(&n6);
    nodes.push_back(&n7);

    Benchmark::SemiHardCodeElement8N<2> e(nodes);

    while (runner.KeepRunningIterations(1e6))
    {
        e.BuildInternalGradient();
    }
}

BENCHMARK(BuildGradient, HardcodeDynamicStaticAlloc, runner)
{
    std::vector<NuTo::NodeBase*> nodes;

    NuTo::NodeDofInfo info;
    info.mDimension = 2;
    info.mNumTimeDerivatives = 0;
    info.mIsDof = true;

    std::map<NuTo::Node::eDof, NuTo::NodeDofInfo> infos;
    infos[NuTo::Node::eDof::COORDINATES] = info;
    infos[NuTo::Node::eDof::DISPLACEMENTS] = info;

    NuTo::NodeDof n0(infos);
    NuTo::NodeDof n1(infos);
    NuTo::NodeDof n2(infos);
    NuTo::NodeDof n3(infos);
    NuTo::NodeDof n4(infos);
    NuTo::NodeDof n5(infos);
    NuTo::NodeDof n6(infos);
    NuTo::NodeDof n7(infos);
    n0.Set(NuTo::Node::eDof::COORDINATES, 0, Eigen::Vector2d({0, 0}));
    n1.Set(NuTo::Node::eDof::COORDINATES, 0, Eigen::Vector2d({1, 0}));
    n2.Set(NuTo::Node::eDof::COORDINATES, 0, Eigen::Vector2d({1, 1}));
    n3.Set(NuTo::Node::eDof::COORDINATES, 0, Eigen::Vector2d({0, 1}));
    n4.Set(NuTo::Node::eDof::COORDINATES, 0, Eigen::Vector2d({0.5, 0}));
    n5.Set(NuTo::Node::eDof::COORDINATES, 0, Eigen::Vector2d({1, 0.5}));
    n6.Set(NuTo::Node::eDof::COORDINATES, 0, Eigen::Vector2d({0.5, 1}));
    n7.Set(NuTo::Node::eDof::COORDINATES, 0, Eigen::Vector2d({0, 0.5}));
    nodes.push_back(&n0);
    nodes.push_back(&n1);
    nodes.push_back(&n2);
    nodes.push_back(&n3);
    nodes.push_back(&n4);
    nodes.push_back(&n5);
    nodes.push_back(&n6);
    nodes.push_back(&n7);

    Benchmark::HardCodeElement8NDynamicStaticAlloc e(nodes);

    while (runner.KeepRunningIterations(1e6))
    {
        e.BuildInternalGradient();
    }
}

BENCHMARK(BuildGradient, HardcodeDynamic, runner)
{
    std::vector<NuTo::NodeBase*> nodes;

    NuTo::NodeDofInfo info;
    info.mDimension = 2;
    info.mNumTimeDerivatives = 0;
    info.mIsDof = true;

    std::map<NuTo::Node::eDof, NuTo::NodeDofInfo> infos;
    infos[NuTo::Node::eDof::COORDINATES] = info;
    infos[NuTo::Node::eDof::DISPLACEMENTS] = info;

    NuTo::NodeDof n0(infos);
    NuTo::NodeDof n1(infos);
    NuTo::NodeDof n2(infos);
    NuTo::NodeDof n3(infos);
    NuTo::NodeDof n4(infos);
    NuTo::NodeDof n5(infos);
    NuTo::NodeDof n6(infos);
    NuTo::NodeDof n7(infos);
    n0.Set(NuTo::Node::eDof::COORDINATES, 0, Eigen::Vector2d({0, 0}));
    n1.Set(NuTo::Node::eDof::COORDINATES, 0, Eigen::Vector2d({1, 0}));
    n2.Set(NuTo::Node::eDof::COORDINATES, 0, Eigen::Vector2d({1, 1}));
    n3.Set(NuTo::Node::eDof::COORDINATES, 0, Eigen::Vector2d({0, 1}));
    n4.Set(NuTo::Node::eDof::COORDINATES, 0, Eigen::Vector2d({0.5, 0}));
    n5.Set(NuTo::Node::eDof::COORDINATES, 0, Eigen::Vector2d({1, 0.5}));
    n6.Set(NuTo::Node::eDof::COORDINATES, 0, Eigen::Vector2d({0.5, 1}));
    n7.Set(NuTo::Node::eDof::COORDINATES, 0, Eigen::Vector2d({0, 0.5}));
    nodes.push_back(&n0);
    nodes.push_back(&n1);
    nodes.push_back(&n2);
    nodes.push_back(&n3);
    nodes.push_back(&n4);
    nodes.push_back(&n5);
    nodes.push_back(&n6);
    nodes.push_back(&n7);

    Benchmark::HardCodeElement8NDynamic e(nodes);

    while (runner.KeepRunningIterations(1e6))
    {
        e.BuildInternalGradient();
    }
}

BENCHMARK(BuildGradient, NuTo, runner)
{
    NuTo::Structure s(2);

    s.NodeCreate(0, Eigen::Vector2d({0, 0}));
    s.NodeCreate(1, Eigen::Vector2d({1, 0}));
    s.NodeCreate(2, Eigen::Vector2d({1, 1}));
    s.NodeCreate(3, Eigen::Vector2d({0, 1}));
    s.NodeCreate(4, Eigen::Vector2d({0.5, 0}));
    s.NodeCreate(5, Eigen::Vector2d({1, .5}));
    s.NodeCreate(6, Eigen::Vector2d({0.5, 1}));
    s.NodeCreate(7, Eigen::Vector2d({0, .5}));


    int myInterpolationType = s.InterpolationTypeCreate("Quad2D");
    s.InterpolationTypeAdd(myInterpolationType, NuTo::Node::eDof::COORDINATES,
                           NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    s.InterpolationTypeAdd(myInterpolationType, NuTo::Node::eDof::DISPLACEMENTS,
                           NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    s.InterpolationTypeSetIntegrationType(myInterpolationType, NuTo::eIntegrationType::IntegrationType2D4NGauss4Ip);

    int elementID = s.ElementCreate(myInterpolationType, {0, 1, 2, 3, 4, 5, 6, 7});

    s.ConstitutiveLawCreate(0, NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    s.ConstitutiveLawSetParameterDouble(0, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, 12);
    s.ConstitutiveLawSetParameterDouble(0, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, 0.3);
    s.ElementTotalSetConstitutiveLaw(0);

    s.ElementSetConstitutiveLaw(elementID, 0);

    auto section = NuTo::SectionPlane::Create(3., true);
    s.ElementTotalSetSection(section);

    s.ElementTotalConvertToInterpolationType();
    std::map<NuTo::ElementEnum::eOutput, std::shared_ptr<NuTo::ElementOutputBase>> elementOutputMap;
    elementOutputMap[NuTo::ElementEnum::eOutput::INTERNAL_GRADIENT] =
            std::make_shared<NuTo::ElementOutputBlockVectorDouble>(s.GetDofStatus());

    NuTo::ElementBase* element = s.ElementGetElementPtr(elementID);

    while (runner.KeepRunningIterations(1e6))
    {
        element->Evaluate(elementOutputMap);
    }
}


BENCHMARK(BuildGradient, NuToPDE, runner)
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
    NuTo::InterpolationQuadSerendipity coordInterpolation(2);
    NuTo::ElementSimple coordElement(coordNodes, coordInterpolation);

    NuTo::NodeSimple nd0(Eigen::Vector2d({0, 0}));
    NuTo::NodeSimple nd1(Eigen::Vector2d({0, 0}));
    NuTo::NodeSimple nd2(Eigen::Vector2d({0, 0}));
    NuTo::NodeSimple nd3(Eigen::Vector2d({0, 0}));
    NuTo::NodeSimple nd4(Eigen::Vector2d({0, 0}));
    NuTo::NodeSimple nd5(Eigen::Vector2d({0, 0}));
    NuTo::NodeSimple nd6(Eigen::Vector2d({0, 0}));
    NuTo::NodeSimple nd7(Eigen::Vector2d({0, 0}));
    std::vector<NuTo::NodeSimple*> displNodes({&nd0, &nd1, &nd2, &nd3, &nd4, &nd5, &nd6, &nd7});
    NuTo::InterpolationQuadSerendipity displInterpolation(2);
    NuTo::ElementSimple displElement(displNodes, displInterpolation);

    NuTo::DofType displDof("Displacements", 2, 0);
    NuTo::DofContainer<NuTo::ElementSimple*> elements;
    elements[displDof] = &displElement;

    NuTo::LinearElasticLaw2D law(20000, 0.3);
    NuTo::IntegrandLinearElastic<2> integrand(displDof, law);
    const NuTo::IntegrationTypeTensorProduct<2> it(2, NuTo::eIntegrationMethod::GAUSS);

    NuTo::Cell<2> cell(coordElement, elements, it, integrand);

    while (runner.KeepRunningIterations(1e6))
    {
        cell.Gradient();
    }
}
