#include "nuto/mechanics/timeintegration/RK4.h"
#include "nuto/mechanics/cell/CellIpData.h"
#include "nuto/mechanics/cell/SimpleAssembler.h"
#include "nuto/mechanics/constraints/Constraints.h"
#include "nuto/mechanics/constraints/ConstraintCompanion.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeTensorProduct.h"
#include "nuto/mechanics/interpolation/InterpolationTrussLobatto.h"
#include "nuto/mechanics/mesh/MeshFemDofConvert.h"
#include "nuto/mechanics/mesh/UnitMeshFem.h"
#include "nuto/mechanics/tools/CellStorage.h"
#include "nuto/visualize/Visualizer.h"
#include "nuto/visualize/VoronoiGeometries.h"

#include "boost/filesystem.hpp"

#include <iostream>

using namespace NuTo;

//! @brief Represents a 1D wave equation example
class Wave1D
{
public:
    Wave1D(int numElements, int order)
        : mMesh(UnitMeshFem::CreateLines(numElements))
        , mOrder(order)
        , mDof("Displacement", 1)
        , mIntegrationType(mOrder + 1, eIntegrationMethod::LOBATTO)
    {

        AddDofInterpolation(&mMesh, mDof, mMesh.CreateInterpolation(InterpolationTrussLobatto(order)));

        mMesh.AllocateDofInstances(mDof, 2);
        mDomain = mMesh.ElementsTotal();
        mDofNodes = mMesh.NodesTotal(mDof);
        mCellGroup = mCells.AddCells(mDomain, mIntegrationType);

        std::string resultDirectory = "/Wave1DResults/";
        bool overwriteResultDirectory = true;

        boost::filesystem::path rootPath = boost::filesystem::initial_path();
        mResultDirectoryFull = rootPath.append(resultDirectory);
        if (boost::filesystem::exists(mResultDirectoryFull)) // does p actually exist?
        {
            if (boost::filesystem::is_directory(mResultDirectoryFull))
            {
                if (overwriteResultDirectory)
                {
                    boost::filesystem::remove_all(mResultDirectoryFull);
                    boost::filesystem::create_directory(mResultDirectoryFull);
                }
            }
        }
        else
        {
            boost::filesystem::create_directory(mResultDirectoryFull);
        }
    }

    void SetValues(std::function<double(double)> func, int instance = 0)
    {
        for (NuTo::ElementCollectionFem& elmColl : mDomain)
        {
            NuTo::ElementFem& elmCoord = elmColl.CoordinateElement();
            NuTo::ElementFem& elmDof = elmColl.DofElement(mDof);
            for (int i = 0; i < elmDof.Interpolation().GetNumNodes(); i++)
            {
                elmDof.GetNode(i).SetValue(0, func(Interpolate(elmCoord, elmDof.Interpolation().GetLocalCoords(i))[0]),
                                           instance);
            }
        }
    }

    void SetDirichletBoundaryLeft(std::function<double(double)> func)
    {
        NodeSimple& leftDofNode = mMesh.NodeAtCoordinate(Eigen::VectorXd::Constant(1, 0.), mDof);
        mConstraints.Add(mDof, Constraint::Value(leftDofNode, func));
    }

    void SetDirichletBoundaryRight(std::function<double(double)> func)
    {
        NodeSimple& rightDofNode = mMesh.NodeAtCoordinate(Eigen::VectorXd::Constant(1, 1.), mDof);
        mConstraints.Add(mDof, Constraint::Value(rightDofNode, func));
    }

    void Solve(int numSteps, double stepSize)
    {
        DofInfo dofInfo = DofNumbering::Build(mDofNodes, mDof, mConstraints);
        int numDofsJ = dofInfo.numIndependentDofs[mDof];
        int numDofsK = dofInfo.numDependentDofs[mDof];
        int numDofs = numDofsJ + numDofsK;
        auto cmat = mConstraints.BuildUnitConstraintMatrix(mDof, numDofs);
        Eigen::PermutationMatrix<Eigen::Dynamic> JKtoGlobal(mConstraints.GetJKNumbering(mDof, numDofs).mIndices);
        Eigen::PermutationMatrix<Eigen::Dynamic> GlobalToJK = JKtoGlobal.inverse();

        SimpleAssembler asmbl = SimpleAssembler(dofInfo);

        // Compute mass matrix
        Eigen::VectorXd lumpedMassMx =
                asmbl.BuildDiagonallyLumpedMatrix(mCellGroup, {mDof}, [&](const CellIpData& cipd) {
                    Eigen::MatrixXd N = cipd.N(mDof);
                    DofMatrix<double> massLocal;
                    massLocal(mDof, mDof) = N.transpose() * N;
                    return massLocal;
                })[mDof];

        // Compute inverse modified mass matrix
        // (this will not work for arbitrary cmat since then the modified mass is in
        // general not diagonal)
        Eigen::VectorXd modMass = (cmat.transpose() * lumpedMassMx.asDiagonal() * cmat).eval().diagonal();

        // Define gradient function
        auto rightHandSide = [&](const CellIpData& cipd) {
            Eigen::MatrixXd B = cipd.B(mDof, Nabla::Gradient());
            DofVector<double> result;
            result[mDof] = -B.transpose() * B * cipd.NodeValueVector(mDof);
            return result;
        };

        // Set up equation system

        auto ODESystem1stOrder = [&](const Eigen::VectorXd& w, Eigen::VectorXd& dwdt, double t) {
            // unpack independent velocities and values
            Eigen::VectorXd valsJ = w.head(numDofsJ);
            Eigen::VectorXd veloJ = w.tail(numDofsJ);
            // Update constrained dofs, get all velocities and values
            auto B = mConstraints.GetSparseGlobalRhs(mDof, numDofs, t);
            Eigen::VectorXd vals = cmat * valsJ + B;
            // Below is actually a dot(B) type thing missing
            Eigen::VectorXd velo = cmat * veloJ;
            // Node Merge
            for (NodeSimple& nd : mDofNodes)
            {
                int dofNr = nd.GetDofNumber(0);
                nd.SetValue(0, vals[dofNr], 0);
                nd.SetValue(0, velo[dofNr], 1);
            }
            Eigen::VectorXd rhsFull = cmat.transpose() * asmbl.BuildVector(mCellGroup, {mDof}, rightHandSide)[mDof];
            dwdt.head(numDofsJ) = veloJ;
            dwdt.tail(numDofsJ) = rhsFull.cwiseQuotient(modMass);
        };

        Eigen::VectorXd state(2 * numDofsJ);
        // Extract values
        for (NodeSimple& nd : mDofNodes)
        {
            int dofNr = nd.GetDofNumber(0);
            int dofNrJK = GlobalToJK.indices()[dofNr];
            if (dofNrJK < numDofsJ)
            {
                state[dofNrJK] = nd.GetValues(0)[0];
                state[dofNrJK + numDofsJ] = nd.GetValues(1)[0];
            }
            else
            {
                // ignore dependent dofs
            }
        }

        TimeIntegration::RK4<Eigen::VectorXd> ti;

        double t = 0.;
        int plotcounter = 1;
        for (int i = 0; i < numSteps; i++)
        {
            t = i * stepSize;
            ti.do_step(ODESystem1stOrder, state, t, stepSize);

            std::cout << i + 1 << std::endl;
            if ((i * 100) % numSteps == 0)
            {
                NuTo::Visualize::Visualizer visualize(
                        mCellGroup, NuTo::Visualize::VoronoiHandler(
                                            Visualize::VoronoiGeometryLine(mOrder + 1, Visualize::LOBATTO)));
                visualize.DofValues(mDof);
                visualize.WriteVtuFile(mResultDirectoryFull.string() + std::string("Wave1D_") +
                                       std::to_string(plotcounter) + ".vtu");
                plotcounter++;
            }
        }
    }

private:
    MeshFem mMesh;
    Group<ElementCollectionFem> mDomain;
    Group<NodeSimple> mDofNodes;

    int mOrder;

    DofType mDof;
    Constraint::Constraints mConstraints;

    IntegrationTypeTensorProduct<1> mIntegrationType;
    CellStorage mCells;
    Group<CellInterface> mCellGroup;

    boost::filesystem::path mResultDirectoryFull;
};

int main(int argc, char* argv[])
{

    // ***************************
    // Some nice looking functions
    // ***************************

    auto cosineBump = [](double x) {
        double a = 0.3;
        double b = 0.7;
        if ((a < x) && (x < b))
        {
            return 0.5 * (1. - cos(2 * M_PI * (x - a) / (b - a)));
        }
        return 0.;
    };

    auto cosineBumpDerivative = [](double x) {
        double a = 0.3;
        double b = 0.7;
        if ((a < x) && (x < b))
        {
            return M_PI / (b - a) * sin(2 * M_PI * (x - a) / (b - a));
        }
        return 0.;
    };

    auto smearedStepFunction = [](double t, double tau) {
        if ((0 < t) && (t < tau))
        {
            return 0.5 * (1. - cos(2 * M_PI * t / tau));
        }
        return 0.;
    };

    // ****************************
    // Solving the 1D wave equation
    // ****************************

    int numElms = 100;
    int interpolationOrder = 3;

    int numSteps = 1000;
    double stepSize = 0.001;
    double tau = 0.3;

    Wave1D example1(numElms, interpolationOrder);
    // example1.SetValues(cosineBump, 0); // Initial data values
    // example1.SetValues(cosineBumpDerivative, 1); // Initial data velocities
    example1.SetDirichletBoundaryLeft([](double /*t*/) { return 0.; });
    example1.SetDirichletBoundaryRight([&](double t) { return smearedStepFunction(t, tau); });
    example1.Solve(numSteps, stepSize);
}
