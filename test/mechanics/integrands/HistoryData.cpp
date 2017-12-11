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
#include "mechanics/mesh/UnitMeshFem.h"
#include "mechanics/nodes/NodeSimple.h"

#include <cassert>
#include <functional>
#include <map>
#include <vector>

using namespace Eigen;
using namespace NuTo;
using namespace NuTo::Groups;

constexpr int cellId = 354;


// Custom Integrand %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

class CustomIntegrand
{
    DofType mDofType;

public:
    CustomIntegrand(DofType dofType)
        : mDofType(dofType)
    {
    }

    NuTo::DofVector<double> Vector(const NuTo::CellData& cellData, const NuTo::CellIpData& cellIpData)
    {
        static std::set<int> ipIds;

        BOOST_CHECK(cellData.GetCellId() == cellId);

        // Each IP Id should only be passed once per cell and values have to be 0 and 1 (linear truss)
        BOOST_CHECK(cellIpData.GetIpId() == 0 || cellIpData.GetIpId() == 1);
        BOOST_CHECK(ipIds.find(cellIpData.GetIpId()) == ipIds.end());
        ipIds.insert(cellIpData.GetIpId());


        NuTo::DofVector<double> vector;
        vector[mDofType] = Eigen::VectorXd::Zero(2);
        return vector;
    }

    NuTo::DofMatrix<double> Matrix(const NuTo::CellData& cellData, const NuTo::CellIpData& cellIpData)
    {
        static std::set<int> ipIds;

        BOOST_CHECK(cellData.GetCellId() == cellId);

        // Each IP Id should only be passed once per cell and values have to be 0 and 1 (linear truss)
        BOOST_CHECK(cellIpData.GetIpId() == 0 || cellIpData.GetIpId() == 1);
        BOOST_CHECK(ipIds.find(cellIpData.GetIpId()) == ipIds.end());
        ipIds.insert(cellIpData.GetIpId());

        NuTo::DofMatrix<double> matrix;
        matrix(mDofType, mDofType) = Eigen::MatrixXd::Zero(2, 2);
        return matrix;
    }
};


// %%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

using namespace std::placeholders;

BOOST_AUTO_TEST_CASE(Pass_Data_To_Integrand)
{

    // Create mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MeshFem mesh = UnitMeshFem::CreateTrusses(1);
    DofType someDof("SomeDof", 1);
    const auto& interpolation = mesh.CreateInterpolation(InterpolationTrussLinear(1));
    AddDofInterpolation(&mesh, someDof, interpolation);


    // Create constraints %%%%%%%%%%%%%%%%%%%%%%%
    ConstraintPde::Constraints constraints;

    // DOF numbering %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DofNumbering::DofInfo dofInfo = DofNumbering::Build(mesh.NodesTotal(someDof), someDof, constraints);


    // Create integrand %%%%%%%%%%%%%%%%%%%%%%%%%
    CustomIntegrand integrand(someDof);


    auto vectorF = std::bind(&CustomIntegrand::Vector, integrand, _1, _2);
    auto matrixF = std::bind(&CustomIntegrand::Matrix, integrand, _1, _2);


    //    // Create integration type %%%%%%%%%%%%%%%%%%
    IntegrationTypeTensorProduct<1> integrationType(2, eIntegrationMethod::GAUSS);


    // Create cell %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    boost::ptr_vector<CellInterface> cellContainer;
    cellContainer.push_back(new Cell(mesh.Elements[0], integrationType, cellId));

    Group<CellInterface> cell;
    cell.Add(cellContainer.back());


    // Build Vector and matrix  %%%%%%%%%%%%%%%%%
    SimpleAssembler assembler(dofInfo.numIndependentDofs, dofInfo.numDependentDofs);
    assembler.BuildVector(cell, {someDof}, vectorF);
    assembler.BuildMatrix(cell, {someDof}, matrixF);
}
