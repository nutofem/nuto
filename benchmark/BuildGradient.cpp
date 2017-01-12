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
#include "mechanics/integrationtypes/IntegrationType2D4NGauss4Ip.h"
#include "mechanics/elements/ElementBase.h"


namespace Benchmark
{

class HardCodeElement8NDynamic
{

public:

    HardCodeElement8NDynamic(const std::vector<NuTo::NodeBase*>& rNodes) : mNodes(rNodes) {}

    Eigen::VectorXd BuildInternalGradient() const
    {
        NuTo::IntegrationType2D4NGauss4Ip it;
        Eigen::Vector2d ip;
        double ipCoords[2];

        Eigen::VectorXd result = Eigen::VectorXd::Zero(16,1);

        for (int i = 0; i < it.GetNumIntegrationPoints(); ++i)
        {
            it.GetLocalIntegrationPointCoordinates2D(i, ipCoords);
            ip[0] = ipCoords[0];
            ip[1] = ipCoords[1];

            auto derivativeShapeFunctions = NuTo::ShapeFunctions2D::DerivativeShapeFunctionsQuadOrder2(ip);
            auto J = GetJacobian(derivativeShapeFunctions);

            auto B = GetB(derivativeShapeFunctions, J);

            result += B.transpose() * (C() * (B * GetDisp())) * J.determinant();
        }
        return result;
    }

    Eigen::VectorXd GetDisp() const
    {
        Eigen::VectorXd disp(16);
        for (int i = 0; i < 8; ++i)
        {
            const auto& nodeValues = mNodes[i]->Get(NuTo::Node::eDof::DISPLACEMENTS);
            disp(2*i) = nodeValues[0];
            disp(2*i+1) = nodeValues[1];
        }
        return disp;
    }

    Eigen::Matrix2d GetJacobian(const Eigen::MatrixXd& rDerivativeShapeFunctions) const
    {
        Eigen::Matrix<double, 2, Eigen::Dynamic> coordinates(2,8);
        for (int i = 0; i < 8; ++i)
        {
            const auto& nodeCoordinate = mNodes[i]->Get(NuTo::Node::eDof::COORDINATES);
            coordinates(0,i) = nodeCoordinate[0];
            coordinates(1,i) = nodeCoordinate[1];
        }
        return coordinates * rDerivativeShapeFunctions;
    }

    Eigen::MatrixXd GetB(const Eigen::MatrixXd& rDerivativeShapeFunctions, const Eigen::Matrix2d& J) const
    {
        Eigen::MatrixXd derivativeShapeFunctionsJ = rDerivativeShapeFunctions * J.inverse();
        Eigen::MatrixXd B = Eigen::MatrixXd::Zero(3,16);

        for (int iNode = 0, iColumn = 0; iNode < 8; ++iNode, iColumn += 2)
        {
            double dNdX = derivativeShapeFunctionsJ(iNode, 0);
            double dNdY = derivativeShapeFunctionsJ(iNode, 1);

            B(0, iColumn) = dNdX;
            B(1, iColumn + 1) = dNdY;
            B(2, iColumn) = dNdY;
            B(2, iColumn + 1) = dNdX;
        }
        return B;
    }

    Eigen::Matrix3d C() const
    {
        Eigen::Matrix3d c = Eigen::Matrix3d::Zero();

        double C11, C12, C33;
        std::tie(C11, C12, C33) = NuTo::EngineeringStressHelper::CalculateCoefficients2DPlaneStress(20000, 0.3);

        c(0, 0) = C11;
        c(1, 0) = C12;
        c(0, 1) = C12;
        c(1, 1) = C11;
        c(2, 2) = C33;
        return c;
    }


private:
    std::vector<NuTo::NodeBase*> mNodes;
};


class HardCodeElement8N
{

public:

    HardCodeElement8N(const std::vector<NuTo::NodeBase*>& rNodes) : mNodes(rNodes) {}

    Eigen::Matrix<double, 16, 1> BuildInternalGradient() const
    {
        NuTo::IntegrationType2D4NGauss4Ip it;
        Eigen::Vector2d ip;
        double ipCoords[2];

        Eigen::Matrix<double, 16, 1> result = Eigen::Matrix<double, 16, 1>::Zero();

        for (int i = 0; i < it.GetNumIntegrationPoints(); ++i)
        {
            it.GetLocalIntegrationPointCoordinates2D(i, ipCoords);
            ip[0] = ipCoords[0];
            ip[1] = ipCoords[1];

            auto derivativeShapeFunctions = NuTo::ShapeFunctions2D::DerivativeShapeFunctionsQuadOrder2(ip);
            auto J = GetJacobian(derivativeShapeFunctions);

            auto B = GetB(derivativeShapeFunctions, J);

            result += B.transpose() * (C() * (B * GetDisp())) * J.determinant();
        }
        return result;
    }

    Eigen::Matrix<double, 16, 1> GetDisp() const
    {
        Eigen::Matrix<double, 16, 1> disp;
        for (int i = 0; i < 8; ++i)
        {
            const auto& nodeValues = mNodes[i]->Get(NuTo::Node::eDof::DISPLACEMENTS);
            disp(2*i) = nodeValues[0];
            disp(2*i+1) = nodeValues[1];
        }
        return disp;
    }

    Eigen::Matrix2d GetJacobian(const Eigen::Matrix<double, 8, 2>& rDerivativeShapeFunctions) const
    {
        Eigen::Matrix<double, 2, 8> coordinates;
        for (int i = 0; i < 8; ++i)
        {
            const auto& nodeCoordinate = mNodes[i]->Get(NuTo::Node::eDof::COORDINATES);
            coordinates(0,i) = nodeCoordinate[0];
            coordinates(1,i) = nodeCoordinate[1];
        }
        return coordinates * rDerivativeShapeFunctions;
    }

    Eigen::Matrix<double, 3, 16> GetB(const Eigen::Matrix<double, 8, 2>& rDerivativeShapeFunctions, const Eigen::Matrix<double, 2, 2>& J) const
    {
        Eigen::Matrix<double, 8, 2> derivativeShapeFunctionsJ = rDerivativeShapeFunctions * J.inverse();
        Eigen::Matrix<double, 3, 16> B = Eigen::Matrix<double, 3, 16>::Zero();

        for (int iNode = 0, iColumn = 0; iNode < 8; ++iNode, iColumn += 2)
        {
            double dNdX = derivativeShapeFunctionsJ(iNode, 0);
            double dNdY = derivativeShapeFunctionsJ(iNode, 1);

            B(0, iColumn) = dNdX;
            B(1, iColumn + 1) = dNdY;
            B(2, iColumn) = dNdY;
            B(2, iColumn + 1) = dNdX;
        }
        return B;
    }

    Eigen::Matrix3d C() const
    {
        Eigen::Matrix3d c = Eigen::Matrix3d::Zero();

        double C11, C12, C33;
        std::tie(C11, C12, C33) = NuTo::EngineeringStressHelper::CalculateCoefficients2DPlaneStress(20000, 0.3);

        c(0, 0) = C11;
        c(1, 0) = C12;
        c(0, 1) = C12;
        c(1, 1) = C11;
        c(2, 2) = C33;
        return c;
    }


private:
    std::vector<NuTo::NodeBase*> mNodes;
};
}


BENCHMARK(BuildGradient, WarmUp, runner)
{
    while (runner.KeepRunningIterations(1e6))
    {
    }
}

BENCHMARK(BuildGradient, HardcodeFix, runner)
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
    n0.Set(NuTo::Node::eDof::COORDINATES, 0, Eigen::Vector2d({0  ,0}));
    n1.Set(NuTo::Node::eDof::COORDINATES, 0, Eigen::Vector2d({1  ,0}));
    n2.Set(NuTo::Node::eDof::COORDINATES, 0, Eigen::Vector2d({1  ,1}));
    n3.Set(NuTo::Node::eDof::COORDINATES, 0, Eigen::Vector2d({0  ,1}));
    n4.Set(NuTo::Node::eDof::COORDINATES, 0, Eigen::Vector2d({0.5,0}));
    n5.Set(NuTo::Node::eDof::COORDINATES, 0, Eigen::Vector2d({1  ,0.5}));
    n6.Set(NuTo::Node::eDof::COORDINATES, 0, Eigen::Vector2d({0.5,1}));
    n7.Set(NuTo::Node::eDof::COORDINATES, 0, Eigen::Vector2d({0  ,0.5}));
    nodes.push_back(&n0);
    nodes.push_back(&n1);
    nodes.push_back(&n2);
    nodes.push_back(&n3);
    nodes.push_back(&n4);
    nodes.push_back(&n5);
    nodes.push_back(&n6);
    nodes.push_back(&n7);

    Benchmark::HardCodeElement8N e(nodes);

    while(runner.KeepRunningIterations(1e6))
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
    n0.Set(NuTo::Node::eDof::COORDINATES, 0, Eigen::Vector2d({0  ,0}));
    n1.Set(NuTo::Node::eDof::COORDINATES, 0, Eigen::Vector2d({1  ,0}));
    n2.Set(NuTo::Node::eDof::COORDINATES, 0, Eigen::Vector2d({1  ,1}));
    n3.Set(NuTo::Node::eDof::COORDINATES, 0, Eigen::Vector2d({0  ,1}));
    n4.Set(NuTo::Node::eDof::COORDINATES, 0, Eigen::Vector2d({0.5,0}));
    n5.Set(NuTo::Node::eDof::COORDINATES, 0, Eigen::Vector2d({1  ,0.5}));
    n6.Set(NuTo::Node::eDof::COORDINATES, 0, Eigen::Vector2d({0.5,1}));
    n7.Set(NuTo::Node::eDof::COORDINATES, 0, Eigen::Vector2d({0  ,0.5}));
    nodes.push_back(&n0);
    nodes.push_back(&n1);
    nodes.push_back(&n2);
    nodes.push_back(&n3);
    nodes.push_back(&n4);
    nodes.push_back(&n5);
    nodes.push_back(&n6);
    nodes.push_back(&n7);

    Benchmark::HardCodeElement8NDynamic e(nodes);

    while(runner.KeepRunningIterations(1e6))
    {
        e.BuildInternalGradient();
    }
}

BENCHMARK(BuildGradient, NuTo, runner)
{
    NuTo::Structure s(2);

    s.NodeCreate(0, Eigen::Vector2d({0,0}));
    s.NodeCreate(1, Eigen::Vector2d({1,0}));
    s.NodeCreate(2, Eigen::Vector2d({1,1}));
    s.NodeCreate(3, Eigen::Vector2d({0,1}));

    int myInterpolationType = s.InterpolationTypeCreate("Quad2D");
    s.InterpolationTypeAdd(myInterpolationType, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    s.InterpolationTypeAdd(myInterpolationType, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    s.InterpolationTypeSetIntegrationType(myInterpolationType, NuTo::eIntegrationType::IntegrationType2D4NGauss4Ip);

    int elementID = s.ElementCreate(myInterpolationType, {0,1,2,3});

    s.ConstitutiveLawCreate(0, NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    s.ConstitutiveLawSetParameterDouble(0,NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, 12);
    s.ConstitutiveLawSetParameterDouble(0,NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO,0.3);
    s.ElementTotalSetConstitutiveLaw(0);

    s.ElementSetConstitutiveLaw(elementID, 0);

    int section = s.SectionCreate("Plane_Stress");
    s.SectionSetThickness(section, 3);
    s.ElementTotalSetSection(section);

    s.ElementTotalConvertToInterpolationType();
    std::map<NuTo::Element::eOutput,std::shared_ptr<NuTo::ElementOutputBase>> elementOutputMap;
    elementOutputMap[NuTo::Element::eOutput::INTERNAL_GRADIENT] = std::make_shared<NuTo::ElementOutputBlockVectorDouble>(s.GetDofStatus());

    NuTo::ElementBase* element = s.ElementGetElementPtr(elementID);

    while(runner.KeepRunningIterations(1e6))
    {
        element->Evaluate(elementOutputMap);
    }
}