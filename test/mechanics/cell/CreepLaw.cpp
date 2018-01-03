#include "BoostUnitTest.h"

#include "base/Group.h"
#include "mechanics/cell/Cell.h"
#include "mechanics/cell/CellInterface.h"
#include "mechanics/cell/SimpleAssember.h"
#include "mechanics/constitutive/LinearElastic.h"
#include "mechanics/constitutive/MechanicsInterface.h"
#include "mechanics/constraints/Constraints.h"
#include "mechanics/constraints/ConstraintCompanion.h"
#include "mechanics/dofs/DofNumbering.h"
#include "mechanics/integrands/MomentumBalance.h"
#include "mechanics/integrationtypes/IntegrationTypeTensorProduct.h"
#include "mechanics/interpolation/InterpolationTrussLinear.h"
#include "mechanics/mesh/MeshFem.h"
#include "mechanics/mesh/MeshFemDofConvert.h"
#include "mechanics/mesh/UnitMeshFem.h"
#include "mechanics/nodes/NodeSimple.h"

#include <cassert>
#include <functional>
#include <map>
#include <vector>

#include <boost/ptr_container/ptr_vector.hpp>

using namespace Eigen;
using namespace NuTo;

constexpr double SpecimenLength = 1.;
constexpr unsigned int numElements = 4;

//! @brief Standard history data management object.
//! @tparam T: Data type of the history data
template <typename T>
class HistoryDataContiguousMemory
{
    unsigned int mIpsPerCell = 0;
    std::vector<T> mHistoryData;

public:
    //! @brief Returns the history data for a specific integration point
    //! @param cellNum: Number of the cell containing the integration point
    //! @param ipNum: Number of the integration point
    //! @return History data of the integration point
    const T& GetIpHistoryData(const unsigned int cellNum, const unsigned int ipNum) const
    {
        assert(cellNum < mHistoryData.size() / mIpsPerCell && "Have you initialized the history data?");
        assert(ipNum < mIpsPerCell);
        return mHistoryData[cellNum * mIpsPerCell + ipNum];
    }

    //! @brief Returns the history data for a specific integration point
    //! @param cellNum: Number of the cell containing the integration point
    //! @param ipNum: Number of the integration point
    //! @return History data of the integration point
    T& GetIpHistoryData(const unsigned int cellNum, const unsigned int ipNum)
    {
        assert(cellNum < mHistoryData.size() / mIpsPerCell && "Have you initialized the history data?");
        assert(ipNum < mIpsPerCell);
        return mHistoryData[cellNum * mIpsPerCell + ipNum];
    }

    //! @brief Initializes the history data
    //! @param numCells: Number of cells that access the history data
    //! @param ipsPerCell: Number of integration points per cell
    void InitializeHistoryData(const unsigned int numCells, const unsigned int ipsPerCell)
    {
        if (mIpsPerCell > 0)
            throw Exception(__PRETTY_FUNCTION__, "History data is already initialized!");
        assert(ipsPerCell > 0);
        assert(numCells > 0);
        mIpsPerCell = ipsPerCell;
        mHistoryData.resize(ipsPerCell * numCells);
    }
};


// %%% Custom law with history data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//! @brief Structure that holds the history data of the creep law
struct CreepHistoryData
{
    EngineeringStrain<1> prevStrain;
    EngineeringStress<1> prevStress;
    VectorXd gamma = Eigen::VectorXd::Zero(2);

    CreepHistoryData()
    {
        prevStrain[0] = 0.;
        prevStress[0] = 0.;
    }
};

//! @brief Creep law using the exponential algorithm
class CreepLaw : public Laws::MechanicsInterface<1>, public HistoryDataContiguousMemory<CreepHistoryData>
{

    double mE; //!< Youngs modulus
    VectorXd mE_KC; //!< Kelvin chain stiffness
    VectorXd mT_KC; //!< Kelvin chain retardation time


public:
    using typename Laws::MechanicsInterface<1>::MechanicsTangent;


    //! @brief Ctor
    //! @param E: Youngs modulus
    //! @param E_KC: Kelvin chain stiffness
    //! @param T_KC: Kelvin chain retardation time
    CreepLaw(double E, VectorXd E_KC, VectorXd T_KC)
        : mE(E)
        , mE_KC(E_KC)
        , mT_KC(T_KC)
    {
        assert(mE_KC.rows() == mT_KC.rows());
    }

    //! @brief Calculates the stress at an integration point
    //! @param strain: Strain at integration point
    //! @param delta_t: Time increment
    //! @param cellNum: Number of currently evaluated cell
    //! @param ipNum: Number of currently evaluated integration point
    //! @return Stress at integration point
    EngineeringStress<1> Stress(EngineeringStrain<1> strain, double delta_t, int cellNum, int ipNum) const override
    {
        // Get history data
        const auto& hisData = GetIpHistoryData(cellNum, ipNum);

        // Calc strain increment
        EngineeringStrain<1> deltaStrain = strain - hisData.prevStrain;

        // Calc creep strain increment
        EngineeringStrain<1> deltaCreep{DeltaCreep(hisData, delta_t)};


        // Calc Stress
        EngineeringStress<1> deltaStress = Tangent(strain, delta_t, cellNum, ipNum) * (deltaStrain - deltaCreep);
        return hisData.prevStress + deltaStress;
    }

    //! @brief Calculates the mechanical tangent(stiffness) at an integration point
    //! @param strain: Strain at integration point
    //! @param delta_t: Time increment
    //! @param cellNum: Number of currently evaluated cell
    //! @param ipNum: Number of currently evaluated integration point
    //! @return Mechanical tangent(stiffness) at an integration point
    MechanicsTangent Tangent(EngineeringStrain<1>, double delta_t, int, int) const override
    {
        // Calc Kelvin Chain compliance
        double chainCompliance = 1. / mE;
        for (unsigned int i = 0; i < mE_KC.rows(); ++i)
            chainCompliance += (1. - Lambda(delta_t, i)) / mE_KC[i];

        // Calc Kelvin Chain stiffness
        return MechanicsTangent::Constant(1. / chainCompliance);
    }

    //! @brief Updates the history data
    //! @param cellData: Cell related data
    //! @param cellIpData: IP related data
    //! @param dofType: Dof type (needed to calculate strains)
    //! @param delta_t: Time increment
    void UpdateHistoryData(const NuTo::CellData& cellData, const NuTo::CellIpData& cellIpData, DofType dofType,
                           double delta_t)
    {
        // Get history data
        auto& hisData = GetIpHistoryData(cellData.GetCellId(), cellIpData.GetIpId());

        // Calculate necessary values for update
        NuTo::BMatrixStrain B = cellIpData.GetBMatrixStrain(dofType);
        NuTo::NodeValues u = cellData.GetNodeValues(dofType);
        EngineeringStrain<1> deltaCreep{DeltaCreep(hisData, delta_t)};
        NuTo::EngineeringStrain<1> strain = B * u;
        NuTo::EngineeringStress<1> stress = Stress(strain, delta_t, cellData.GetCellId(), cellIpData.GetIpId());
        MechanicsTangent E = Tangent(strain, delta_t, cellData.GetCellId(), cellIpData.GetIpId());
        NuTo::EngineeringStrain<1> deltaStrain = strain - hisData.prevStrain;

        // The actual update
        hisData.prevStrain = strain;
        hisData.prevStress = stress;
        for (unsigned int i = 0; i < mE_KC.rows(); ++i)
            hisData.gamma[i] = Lambda(delta_t, i) * E / mE_KC[i] * (deltaStrain - deltaCreep) +
                               Beta(delta_t, i) * hisData.gamma[i];
    }

private:
    //! @brief Calculates the algorithm specific parameter beta
    //! @param delta_t: Time increment
    //! @param index: Index of the Kelvin Unit
    //! @return Algorithm specific parameter beta
    double Beta(double delta_t, unsigned int index) const
    {
        assert(index < mT_KC.rows());
        return std::exp(-delta_t / mT_KC[index]);
    }

    //! @brief Calculates the algorithm specific parameter lambda
    //! @param delta_t: Time increment
    //! @param index: Index of the Kelvin Unit
    //! @return Algorithm specific parameter lambda
    double Lambda(double delta_t, unsigned int index) const
    {
        assert(index < mT_KC.rows());
        return mT_KC[index] / delta_t * (1 - Beta(delta_t, index));
    }

    //! @brief Calculates the creep strain increment
    //! @param hisData: History data object
    //! @param delta_t: Time increment
    //! @return Creep strain increment
    EngineeringStrain<1> DeltaCreep(const CreepHistoryData& hisData, double delta_t) const
    {
        EngineeringStrain<1> deltaCreep;
        deltaCreep[0] = 0.;
        for (unsigned int i = 0; i < mE_KC.rows(); ++i)
            deltaCreep[0] += (1. - Beta(delta_t, i)) * hisData.gamma[i];
        return deltaCreep;
    }
};


// %%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

using namespace std::placeholders;

BOOST_AUTO_TEST_CASE(History_Data)
{
    // Create mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MeshFem mesh = UnitMeshFem::CreateLines(numElements);
    //    MeshFem mesh = UnitMeshFem::Transform(UnitMeshFem::CreateLines(numElements),
    //                                          [](Eigen::VectorXd vec) { return SpecimenLength * vec; });

    DofType displ("displacements", 1);
    const auto& interpolation = mesh.CreateInterpolation(InterpolationTrussLinear());
    AddDofInterpolation(&mesh, displ, interpolation);


    // Create constraints %%%%%%%%%%%%%%%%%%%%%%%
    Constraint::Constraints constraints;
    auto& nodeLeft = mesh.NodeAtCoordinate(Eigen::VectorXd::Zero(1), displ);
    constraints.Add(displ, Constraint::Component(nodeLeft, {eDirection::X}));


    // DOF numbering %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DofNumbering::DofInfo dofInfo = DofNumbering::Build(mesh.NodesTotal(displ), displ, constraints);


    // Create law %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    constexpr double E = 40000;
    const Vector2d E_KC(60000., 120000.);
    const Vector2d D_KC(1., 2.);
    CreepLaw creepLaw(E, E_KC, D_KC);


    // Create integrand %%%%%%%%%%%%%%%%%%%%%%%%%
    Integrands::MomentumBalance<1> momentumBalance(displ, creepLaw);
    double delta_t = 0.1;
    auto MomentumGradientF =
            std::bind(&Integrands::MomentumBalance<1>::Gradient, momentumBalance, _1, _2, std::ref(delta_t));
    auto MomentumHessian0F =
            std::bind(&Integrands::MomentumBalance<1>::Hessian0, momentumBalance, _1, _2, std::ref(delta_t));
    auto MomentumUpdateHistoryDataF =
            std::bind(&CreepLaw::UpdateHistoryData, std::ref(creepLaw), _1, _2, displ, std::ref(delta_t));

    // Create integration type %%%%%%%%%%%%%%%%%%
    IntegrationTypeTensorProduct<1> integrationType(2, eIntegrationMethod::GAUSS);


    // Create cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    boost::ptr_vector<CellInterface> cellContainer;
    Group<CellInterface> momentumBalanceCells;
    int cellId = 0;
    for (ElementCollection& element : mesh.Elements)
    {
        cellContainer.push_back(new Cell(element, integrationType, cellId++));
        momentumBalanceCells.Add(cellContainer.back());
    }

    // Initialize IP data %%%%%%%%%%%%%%%%%%%%%%%
    creepLaw.InitializeHistoryData(momentumBalanceCells.Size(), integrationType.GetNumIntegrationPoints());


    // Get gradient %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SimpleAssembler assembler(dofInfo.numIndependentDofs, dofInfo.numDependentDofs);
    VectorXd displacements = VectorXd::Zero(dofInfo.numIndependentDofs[displ]);

    // Build external Force %%%%%%%%%%%%%%%%%%%%%
    constexpr double rhsForce = 2000.;
    GlobalDofVector extF;
    extF.J[displ].setZero(dofInfo.numIndependentDofs[displ]);
    extF.K[displ].setZero(dofInfo.numDependentDofs[displ]);
    NodeSimple& nodeRight = mesh.NodeAtCoordinate(Eigen::VectorXd::Ones(1) * SpecimenLength, displ);
    extF.J[displ][nodeRight.GetDofNumber(0)] = rhsForce;

    // Post processing stuff %%%%%%%%%%%%%%%%%%%%
    std::vector<double> rhsDispNumerical;

    // Time loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    constexpr double timeFinal = 10.;
    constexpr unsigned int maxIter = 10;
    double time = 0.;
    delta_t = 0.01;

    while (time < timeFinal)
    {
        unsigned int numIter = 0;
        time += delta_t;

        // Calculate residual %%%%%%%%%%%%%%%%%%%
        GlobalDofVector gradient = assembler.BuildVector(momentumBalanceCells, {displ}, MomentumGradientF);
        Eigen::VectorXd residual = gradient.J[displ] - extF.J[displ];

        // Iterate for equilibrium %%%%%%%%%%%%%%
        while (numIter < maxIter && residual.lpNorm<Infinity>() > 1e-9)
        {
            numIter++;
            // Build and solve system %%%%%%%%%%%
            GlobalDofMatrixSparse hessian = assembler.BuildMatrix(momentumBalanceCells, {displ}, MomentumHessian0F);
            Eigen::MatrixXd hessianDense(hessian.JJ(displ, displ));
            Eigen::VectorXd deltaDisplacements = hessianDense.ldlt().solve(residual);
            displacements -= deltaDisplacements;

            // Merge dof values %%%%%%%%%%%%%%%%%
            int numUnconstrainedDofs = dofInfo.numIndependentDofs[displ];
            for (NodeSimple& node : mesh.NodesTotal(displ))
            {
                int dofNumber = node.GetDofNumber(0);
                if (dofNumber < numUnconstrainedDofs)
                    node.SetValue(0, displacements[dofNumber]);
            }

            // Calculate new residual %%%%%%%%%%%
            gradient = assembler.BuildVector(momentumBalanceCells, {displ}, MomentumGradientF);
            residual = gradient.J[displ] - extF.J[displ];
        }
        if (numIter >= maxIter)
        {
            std::cout << residual.lpNorm<Infinity>() << std::endl;
            std::cout << time << std::endl;
            throw Exception(__PRETTY_FUNCTION__, "No convergence");
        }
        // Update history data %%%%%%%%%%%%%%%%%%
        for (auto& cell : momentumBalanceCells)
            cell.Apply(MomentumUpdateHistoryDataF);

        // Store rhs displacement %%%%%%%%%%%%%%%
        if (std::abs(time - std::round(time)) < delta_t / 2.) // <--- store only if time is an integer
            rhsDispNumerical.push_back(displacements[displacements.rows() - 1]);
    }

    // Theoretical solution %%%%%%%%%%%%%%%%%%%%%
    delta_t = 1.;
    std::vector<double> rhsDispTheoretical;
    for (float time = delta_t; time <= timeFinal; time += delta_t)
    {
        double totalStrain = 0.0;
        totalStrain += rhsForce / E;
        for (unsigned int i = 0; i < E_KC.rows(); ++i)
            totalStrain += rhsForce / E_KC[i] * (1. - std::exp(-time / D_KC[i]));
        rhsDispTheoretical.push_back(totalStrain * SpecimenLength);
    }

    // Compare results %%%%%%%%%%%%%%%%%%%%%%%%%%

    assert(rhsDispNumerical.size() == rhsDispTheoretical.size());
    double maxDeviation = 0.;
    for (unsigned int i = 0; i < rhsDispNumerical.size(); ++i)
    {
        double deviation = std::abs((rhsDispNumerical[i] - rhsDispTheoretical[i]) / rhsDispTheoretical[i]);
        if (deviation > maxDeviation)
            maxDeviation = deviation;
    }
    if (maxDeviation > 0.01)
        throw Exception(__PRETTY_FUNCTION__, "Difference between theoretical and numerical solution is too high.");
}
