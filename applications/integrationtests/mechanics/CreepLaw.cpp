#include "BoostUnitTest.h"

#include "nuto/base/Group.h"
#include "nuto/math/EigenSparseSolve.h"

#include "nuto/mechanics/cell/Cell.h"
#include "nuto/mechanics/cell/CellInterface.h"
#include "nuto/mechanics/cell/SimpleAssembler.h"
#include "nuto/mechanics/constitutive/LinearElastic.h"
#include "nuto/mechanics/constitutive/MechanicsInterface.h"
#include "nuto/mechanics/constraints/Constraints.h"
#include "nuto/mechanics/constraints/ConstraintCompanion.h"
#include "nuto/mechanics/dofs/DofNumbering.h"
#include "nuto/mechanics/integrands/MomentumBalance.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeTensorProduct.h"
#include "nuto/mechanics/interpolation/InterpolationTrussLinear.h"
#include "nuto/mechanics/mesh/MeshFem.h"
#include "nuto/mechanics/mesh/GeometryMeshFem.h"
#include "nuto/mechanics/mesh/MeshFemDofConvert.h"
#include "nuto/mechanics/mesh/UnitMeshFem.h"
#include "nuto/mechanics/nodes/DofNode.h"

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
//! @tparam T Data type of the history data
template <typename T>
class HistoryDataContiguousMemory
{
    unsigned int mIpsPerCell = 0;
    std::vector<T> mHistoryData;

public:
    //! @brief Returns the history data for a specific integration point
    //! @param cellNum Number of the cell containing the integration point
    //! @param ipNum Number of the integration point
    //! @return History data of the integration point
    const T& GetIpHistoryData(const unsigned int cellNum, const unsigned int ipNum) const
    {
        assert(cellNum < mHistoryData.size() / mIpsPerCell && "Have you initialized the history data?");
        assert(ipNum < mIpsPerCell);
        return mHistoryData[cellNum * mIpsPerCell + ipNum];
    }

    //! @brief Returns the history data for a specific integration point
    //! @param cellNum Number of the cell containing the integration point
    //! @param ipNum Number of the integration point
    //! @return History data of the integration point
    T& GetIpHistoryData(const unsigned int cellNum, const unsigned int ipNum)
    {
        assert(cellNum < mHistoryData.size() / mIpsPerCell && "Have you initialized the history data?");
        assert(ipNum < mIpsPerCell);
        return mHistoryData[cellNum * mIpsPerCell + ipNum];
    }

    //! @brief Initializes the history data
    //! @param numCells Number of cells that access the history data
    //! @param ipsPerCell Number of integration points per cell
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
    //! @brief Ctor
    //! @param E Youngs modulus
    //! @param E_KC Kelvin chain stiffness
    //! @param T_KC Kelvin chain retardation time
    CreepLaw(double E, VectorXd E_KC, VectorXd T_KC)
        : mE(E)
        , mE_KC(E_KC)
        , mT_KC(T_KC)
    {
        assert(mE_KC.rows() == mT_KC.rows());
    }

    //! @brief Calculates the stress at an integration point
    //! @param strain Strain at integration point
    //! @param delta_t Time increment
    //! @param ids Number of currently evaluated cell and Number of currently evaluated integration point
    //! @return Stress at integration point
    EngineeringStress<1> Stress(EngineeringStrain<1> strain, double delta_t, CellIds ids) const override
    {
        // Get history data
        const auto& hisData = GetIpHistoryData(ids.cellId, ids.ipId);

        // Calc strain increment
        EngineeringStrain<1> deltaStrain = strain - hisData.prevStrain;

        // Calc creep strain increment
        EngineeringStrain<1> deltaCreep{DeltaCreep(hisData, delta_t)};


        // Calc Stress
        EngineeringStress<1> deltaStress = Tangent(strain, delta_t, ids) * (deltaStrain - deltaCreep);
        return hisData.prevStress + deltaStress;
    }

    //! @brief Calculates the mechanical tangent(stiffness) at an integration point
    //! @param delta_t Time increment
    //! @return Mechanical tangent(stiffness) at an integration point
    EngineeringTangent<1> Tangent(EngineeringStrain<1>, double delta_t, CellIds) const override
    {
        // Calc Kelvin Chain compliance
        double chainCompliance = 1. / mE;
        for (unsigned int i = 0; i < mE_KC.rows(); ++i)
            chainCompliance += (1. - Lambda(delta_t, i)) / mE_KC[i];

        // Calc Kelvin Chain stiffness
        return EngineeringTangent<1>::Constant(1. / chainCompliance);
    }

    //! @brief Updates the history data
    //! @param cellIpData IP related data
    //! @param dofType Dof type (needed to calculate strains)
    //! @param delta_t Time increment
    void UpdateHistoryData(const NuTo::CellIpData& cellIpData, DofType dofType, double delta_t)
    {
        // Get history data
        auto& hisData = GetIpHistoryData(cellIpData.Ids().cellId, cellIpData.Ids().ipId);

        // Calculate necessary values for update
        Eigen::MatrixXd B = cellIpData.B(dofType, Nabla::Strain());
        EngineeringStrain<1> deltaCreep{DeltaCreep(hisData, delta_t)};
        EngineeringStrain<1> strain = cellIpData.Apply(dofType, Nabla::Strain());
        EngineeringStress<1> stress = Stress(strain, delta_t, cellIpData.Ids());
        EngineeringTangent<1> E = Tangent(strain, delta_t, cellIpData.Ids());
        EngineeringStrain<1> deltaStrain = strain - hisData.prevStrain;

        // The actual update
        hisData.prevStrain = strain;
        hisData.prevStress = stress;
        for (unsigned int i = 0; i < mE_KC.rows(); ++i)
            hisData.gamma[i] = Lambda(delta_t, i) * E / mE_KC[i] * (deltaStrain - deltaCreep) +
                               Beta(delta_t, i) * hisData.gamma[i];
    }

private:
    //! @brief Calculates the algorithm specific parameter beta
    //! @param delta_t Time increment
    //! @param index Index of the Kelvin Unit
    //! @return Algorithm specific parameter beta
    double Beta(double delta_t, unsigned int index) const
    {
        assert(index < mT_KC.rows());
        return std::exp(-delta_t / mT_KC[index]);
    }

    //! @brief Calculates the algorithm specific parameter lambda
    //! @param delta_t Time increment
    //! @param index Index of the Kelvin Unit
    //! @return Algorithm specific parameter lambda
    double Lambda(double delta_t, unsigned int index) const
    {
        assert(index < mT_KC.rows());
        return mT_KC[index] / delta_t * (1 - Beta(delta_t, index));
    }

    //! @brief Calculates the creep strain increment
    //! @param hisData History data object
    //! @param delta_t Time increment
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
    GeometryMeshFem geoMesh = UnitMeshFem::CreateLines(numElements);
    MeshFem mesh(geoMesh);

    DofType displ("displacements", 1);
    const auto& interpolation = mesh.CreateInterpolation(InterpolationTrussLinear());
    AddDofInterpolation(&mesh, displ, interpolation);


    // Create constraints %%%%%%%%%%%%%%%%%%%%%%%
    Constraint::Constraints constraints;
    auto& nodeLeft = mesh.NodeAtCoordinate(Eigen::VectorXd::Zero(1), displ);
    constraints.Add(displ, Constraint::Component(nodeLeft, {eDirection::X}));


    // DOF numbering %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DofInfo dofInfo = DofNumbering::Build(mesh.NodesTotal(displ), displ, constraints);
    const int numDofs = dofInfo.numIndependentDofs[displ] + dofInfo.numDependentDofs[displ];


    // Create law %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    constexpr double E = 40000;
    const Vector2d E_KC(60000., 120000.);
    const Vector2d D_KC(1., 2.);
    CreepLaw creepLaw(E, E_KC, D_KC);


    // Create integrand %%%%%%%%%%%%%%%%%%%%%%%%%
    Integrands::MomentumBalance<1> momentumBalance(displ, creepLaw);
    double delta_t = 0.1;
    auto MomentumGradientF =
            std::bind(&Integrands::MomentumBalance<1>::Gradient, momentumBalance, _1, std::ref(delta_t));
    auto MomentumHessian0F =
            std::bind(&Integrands::MomentumBalance<1>::Hessian0, momentumBalance, _1, std::ref(delta_t));
    auto MomentumUpdateHistoryDataF =
            std::bind(&CreepLaw::UpdateHistoryData, std::ref(creepLaw), _1, displ, std::ref(delta_t));

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
    SimpleAssembler assembler(dofInfo);
    DofVector<double> solution;
    solution[displ].setZero(dofInfo.numIndependentDofs[displ] + dofInfo.numDependentDofs[displ]);

    // Build external Force %%%%%%%%%%%%%%%%%%%%%
    constexpr double rhsForce = 2000.;

    DofVector<double> extF;
    extF[displ].setZero(dofInfo.numIndependentDofs[displ] + dofInfo.numDependentDofs[displ]);
    DofNode& nodeRight = mesh.NodeAtCoordinate(Eigen::VectorXd::Ones(1) * SpecimenLength, displ);
    extF[displ][nodeRight.GetDofNumber(0)] = rhsForce;

    // Post processing stuff %%%%%%%%%%%%%%%%%%%%
    std::vector<double> rhsDispNumerical;

    // Time loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    constexpr double timeFinal = 10.;
    constexpr unsigned int maxIter = 10;
    double time = 0.;
    delta_t = 0.01;

    auto cMatUnit(constraints.BuildUnitConstraintMatrix(
            {displ}, dofInfo.numIndependentDofs[displ] + dofInfo.numDependentDofs[displ]));

    while (time < timeFinal)
    {
        unsigned int numIter = 0;
        time += delta_t;

        // Calculate residual %%%%%%%%%%%%%%%%%%%
        DofVector<double> gradient = assembler.BuildVector(momentumBalanceCells, {displ}, MomentumGradientF);
        Eigen::VectorXd residualMod = cMatUnit.transpose() * ((gradient - extF)[displ]);

        // Iterate for equilibrium %%%%%%%%%%%%%%
        while (numIter < maxIter && residualMod.lpNorm<Infinity>() > 1e-9)
        {
            numIter++;
            // Build and solve system %%%%%%%%%%%
            DofMatrixSparse<double> hessian = assembler.BuildMatrix(momentumBalanceCells, {displ}, MomentumHessian0F);

            // the following line should all go to the solve routine, no need to deal with modified vectors and matrices
            auto hessianMod = cMatUnit.transpose() * hessian(displ, displ) * cMatUnit;
            Eigen::VectorXd deltaDisplacementsMod = EigenSparseSolve(hessianMod, residualMod, std::string("MumpsLDLT"));

            // compute full solution vector (dependent and independent dofs)
            Eigen::VectorXd deltaDisplacements =
                    cMatUnit * -deltaDisplacementsMod + constraints.GetSparseGlobalRhs(displ, numDofs, time);

            // for correct size, copy the inactive dofs
            solution[displ] += deltaDisplacements;

            // Merge dof values %%%%%%%%%%%%%%%%%
<<<<<<< HEAD
<<<<<<< HEAD
            for (NodeSimple& node : mesh.NodesTotal(displ))
=======
            int numUnconstrainedDofs = dofInfo.numIndependentDofs[displ];
=======
>>>>>>> Merge current PDE reviewed
                for (DofNode& node : mesh.NodesTotal(displ))
>>>>>>> Rename Nodes to DofNode and CoordinateNode
                {
                    int dofNumber = node.GetDofNumber(0);
                    node.SetValue(0, solution[displ][dofNumber]);
                }

            // Calculate new residual %%%%%%%%%%%
            gradient = assembler.BuildVector(momentumBalanceCells, {displ}, MomentumGradientF);
            residualMod = cMatUnit.transpose() * ((gradient - extF)[displ]);
        }
        if (numIter >= maxIter)
        {
            std::cout << residualMod.lpNorm<Infinity>() << std::endl;
            std::cout << time << std::endl;
            throw Exception(__PRETTY_FUNCTION__, "No convergence");
        }
        // Update history data %%%%%%%%%%%%%%%%%%
        for (auto& cell : momentumBalanceCells)
            cell.Apply(MomentumUpdateHistoryDataF);

        // Store rhs displacement %%%%%%%%%%%%%%%
        if (std::abs(time - std::round(time)) < delta_t / 2.) // <--- store only if time is an integer
            rhsDispNumerical.push_back(solution[displ][nodeRight.GetDofNumber(0)]);
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
