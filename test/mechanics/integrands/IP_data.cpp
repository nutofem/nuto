#include "BoostUnitTest.h"

#include "base/Group.h"
#include "mechanics/cell/Cell.h"
#include "mechanics/cell/CellInterface.h"
#include "mechanics/cell/SimpleAssember.h"
#include "mechanics/constitutive/laws/LinearElastic.h"
#include "mechanics/constitutive/laws/MechanicsInterface.h"
#include "mechanics/constraintsPde/Constraints.h"
#include "mechanics/constraintsPde/ConstraintCompanion.h"
#include "mechanics/dofs/DofNumbering.h"
#include "mechanics/integrands/MomentumBalance.h"
#include "mechanics/integrationtypes/IntegrationTypeTensorProduct.h"
#include "mechanics/interpolation/InterpolationTrussLinear.h"
#include "mechanics/mesh/MeshFem.h"
#include "mechanics/mesh/MeshFemDofConvert.h"
#include "mechanics/nodes/NodeSimple.h"

#include <map>
#include <functional>
#include <vector>

using namespace NuTo;
using namespace NuTo::Groups;


class CreepLaw : public Laws::MechanicsInterface<1>
{

    double mE;
    double mNu;
    std::map<int, std::vector<float>> mIPdata;

public:
    using typename Laws::MechanicsInterface<1>::MechanicsTangent;


    CreepLaw(double E, double Nu)
        : mE(E)
        , mNu(Nu)
    {
    }

    void InitializeIPData(const Groups::Group<CellInterface>& cells, DofType dof)
    {
        for (CellInterface& cell : cells)
        {

            cell.Apply([this](const NuTo::CellData& cellData, const NuTo::CellIpData& cellIpData) {
                if (mIPdata.find(cellData.GetCellId()) == mIPdata.end())
                    mIPdata.emplace(cellData.GetCellId(), std::vector<float>(cellData.GetNumIntegrationPoints()));
            });
        }
    }

    EngineeringStressPDE<1> Stress(EngineeringStrainPDE<1> strain, double delta_t, int cellNum,
                                   int ipNum) const override
    {
        return mE * strain;
    }

    MechanicsTangent Tangent(EngineeringStrainPDE<1>, double delta_t, int cellNum, int ipNum) const override
    {
        return MechanicsTangent::Constant(mE);
    }
};


MeshFem Mesh1D()
{
    MeshFem mesh;
    const InterpolationSimple& interpolation = mesh.CreateInterpolation(InterpolationTrussLinear(1));
    NodeSimple* nr = nullptr;
    for (unsigned int i = 0; i < 21; ++i)
    {
        NodeSimple& nl = mesh.Nodes.Add({i * 0.5});
        if (i > 0)
            mesh.Elements.Add({{{nl, *nr}, interpolation}});
        nr = &nl;
    }


    return mesh;
}

using namespace std::placeholders;

BOOST_AUTO_TEST_CASE(IP_data)
{
    // Create mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MeshFem mesh = Mesh1D();
    DofType displ("displacements", 1);
    const auto& interpolation = mesh.CreateInterpolation(InterpolationTrussLinear(1));
    AddDofInterpolation(&mesh, displ, interpolation);


    // Create constraints %%%%%%%%%%%%%%%%%%%%%%%
    ConstraintPde::Constraints constraints;
    auto& nodeLeft = mesh.NodeAtCoordinate(Eigen::VectorXd::Zero(1), displ);
    constraints.Add(displ, ConstraintPde::Component(nodeLeft, {eDirection::X}));


    // DOF numbering %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DofNumbering::DofInfo dofInfo = DofNumbering::Build(mesh.NodesTotal(displ), displ, constraints);


    // Create law %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    constexpr double E = 20000;
    constexpr double nu = 0.2;
    Laws::LinearElastic<1> linearElasticLaw(E, nu);
    CreepLaw creepLaw(E, nu);


    // Create integrand %%%%%%%%%%%%%%%%%%%%%%%%%
    Integrands::MomentumBalance<1> momentumBalance(displ, creepLaw);
    auto MomentumGradientF = std::bind(&Integrands::MomentumBalance<1>::Gradient, momentumBalance, _1, _2, 0.);
    auto MomentumHessian0F = std::bind(&Integrands::MomentumBalance<1>::Hessian0, momentumBalance, _1, _2, 0.);


    // Create integration type %%%%%%%%%%%%%%%%%%
    IntegrationTypeTensorProduct<1> integrationType(2, eIntegrationMethod::GAUSS);


    // Create cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    boost::ptr_vector<CellInterface> cellContainer;
    Group<CellInterface> momentumBalanceCells;
    for (ElementCollection& element : mesh.Elements)
    {
        cellContainer.push_back(new Cell(element, integrationType));
        momentumBalanceCells.Add(cellContainer.back());
    }

    // Initialize IP data %%%%%%%%%%%%%%%%%%%%%%%
    creepLaw.InitializeIPData(momentumBalanceCells, displ);


    // Assemble system %%%%%%%%%%%%%%%%%%%%%%%%%%
    SimpleAssembler assembler(dofInfo.numIndependentDofs, dofInfo.numDependentDofs);
    GlobalDofVector gradient = assembler.BuildVector(momentumBalanceCells, {displ}, MomentumGradientF);
    GlobalDofMatrixSparse hessian = assembler.BuildMatrix(momentumBalanceCells, {displ}, MomentumHessian0F);


    // Build external Force %%%%%%%%%%%%%%%%%%%%%
    GlobalDofVector extF;
    extF.J[displ].setZero(dofInfo.numIndependentDofs[displ]);
    extF.K[displ].setZero(dofInfo.numDependentDofs[displ]);
    NodeSimple& nodeRight = mesh.NodeAtCoordinate(Eigen::VectorXd::Ones(1) * 10, displ);
    extF.J[displ][nodeRight.GetDofNumber(0)] = 2000.;


    // Solve system %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Eigen::MatrixXd hessianDense(hessian.JJ(displ, displ));
    Eigen::VectorXd newDisplacements = hessianDense.ldlt().solve(gradient.J[displ] - extF.J[displ]);
}
