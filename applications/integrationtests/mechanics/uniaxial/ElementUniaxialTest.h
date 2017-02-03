/*
 * ElementUniaxialTest.h
 *
 *  Created on: 13 May 2015
 *      Author: ttitsche
 */

#pragma once

#include <iomanip>

#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>

#include "mechanics/mesh/MeshGenerator.h"

#include "mechanics/MechanicsEnums.h"

#include "math/SparseMatrixCSRVector2.h"
#include "mechanics/dofSubMatrixStorage/BlockScalar.h"
#include "mechanics/elements/ElementBase.h"
#include "mechanics/interpolationtypes/InterpolationBase.h"
#include "mechanics/interpolationtypes/InterpolationType.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/structures/StructureOutputBlockMatrix.h"
#include "mechanics/structures/StructureOutputBlockVector.h"
#include "math/SparseDirectSolverMUMPS.h"


#ifdef ENABLE_VISUALIZE
#include "visualize/VisualizeEnum.h"
#endif

/*
 *    Interface for displacement controlled uni-axial tensile test
 */
namespace NuToTest
{

class ElementUniaxialTest
{
public:

    double lX = 365;
    double lY = 13;
    double lZ = 37;

    double E = 42;
    double nu = 0.26;
    double deltaL = 0.6174;
    double rho = M_PI;

    //! @brief provide a structure that contains the geometry
    //! of a uniaxial test, starting from (0,0,0) to (lx, ly, lz)
    void Run(NuTo::Structure& rS, const std::string& rVisualizationDirectory)
    {

        CheckVolume(rS);

        rS.SetVerboseLevel(10);
        rS.SetShowTime(true);


        SetConstitutiveLaw(rS);
        SetBoundaryConditions(rS);

        rS.NodeBuildGlobalDofs();
        rS.CalculateMaximumIndependentSets();

        rS.NodeInfo(10);

        CheckStiffness(rS);
        Solve(rS);
        CheckSolution(rS);
        CheckMass(rS);
        Visualize(rS, rVisualizationDirectory);
    }

private:

    void CheckVolume(NuTo::Structure& rS)
    {
        double volume = rS.ElementGroupGetVolume(rS.GroupGetElementsTotal());
        if (rS.GetDimension() == 2)
            volume *= lZ;
        if (rS.GetDimension() == 1)
            volume *= lY * lZ;

        BOOST_CHECK_CLOSE(volume, lX * lY * lZ, 1.e-6);
    }

    void SetConstitutiveLaw(NuTo::Structure& rS)
    {
        rS.ConstitutiveLawCreate(0, NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
        rS.ConstitutiveLawSetParameterDouble(0,NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, E);
        rS.ConstitutiveLawSetParameterDouble(0,NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO,nu);
        rS.ConstitutiveLawSetParameterDouble(0,NuTo::Constitutive::eConstitutiveParameter::DENSITY, rho);
        rS.ElementTotalSetConstitutiveLaw(0);
    }

    void SetBoundaryConditions(NuTo::Structure& rS)
    {

        int dimension = rS.GetDimension();

        Eigen::MatrixXd directions = Eigen::MatrixXd::Identity(dimension, dimension);

        double eps = 1.e-3;

        // fix origin
        Eigen::VectorXd origin(dimension);
        origin.setZero();
        int nodeGroupOrigin = rS.GroupCreate("Nodes");
        rS.GroupAddNodeRadiusRange(nodeGroupOrigin, origin, 0, eps);

        if (rS.GroupGetNumMembers(nodeGroupOrigin) != 1)
            throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Node at origin (0,0,0) does not exist.");

        int nodeOrigin = rS.GroupGetMemberIds(nodeGroupOrigin)[0];
        for (int iDim = 1; iDim < dimension; ++iDim)
            rS.ConstraintLinearSetDisplacementNode(nodeOrigin, directions.col(iDim), 0.0);

        if (dimension == 3)
        {
        // fix rotation
            Eigen::VectorXd originRot(3);
            originRot << 0, 0, lZ;
            int nodeGroupOriginRot = rS.GroupCreate(NuTo::eGroupId::Nodes);
            rS.GroupAddNodeRadiusRange(nodeGroupOriginRot, originRot, 0, eps);

            int nodeOriginRot = rS.GroupGetMemberIds(nodeGroupOriginRot)[0];
            rS.ConstraintLinearSetDisplacementNode(nodeOriginRot, directions.col(1), 0.0);
        }

        // fix x = 0 plane
        int nodesX0 = rS.GroupCreate(NuTo::eGroupId::Nodes);
        rS.GroupAddNodeCoordinateRange(nodesX0, 0, -eps, eps);
        rS.ConstraintLinearSetDisplacementNodeGroup(nodesX0, directions.col(0), 0.);

        // apply displacement on x = lX plane
        int nodesXlX = rS.GroupCreate(NuTo::eGroupId::Nodes);
        rS.GroupAddNodeCoordinateRange(nodesXlX, 0, lX-eps, lX+eps);
        rS.ConstraintLinearSetDisplacementNodeGroup(nodesXlX, directions.col(0), deltaL);

    }

    void CheckStiffness(NuTo::Structure& rS)
    {
        if (rS.GetNumTotalDofs() > 50)
            return;
        BOOST_CHECK(rS.ElementCheckHessian0(1.e-8, 1.e-6, true));
        BOOST_CHECK(rS.CheckHessian0(1.e-8, 1.e-6, true));
    }

    void Solve(NuTo::Structure& rS)
    {
        rS.SolveGlobalSystemStaticElastic();
        auto internalGradient = rS.BuildGlobalInternalGradient();
        BOOST_CHECK_SMALL(internalGradient.J.CalculateNormL2()[NuTo::Node::eDof::DISPLACEMENTS], 1e-8);
    }

    void CheckSolution(NuTo::Structure& rS)
    {
        double analyticStrainX = deltaL / lX;
        double analyticStressX = analyticStrainX * E;
        double analyticForce = analyticStressX * lY * lZ;

        std::cout << "\n  ### analytical values ###" << std::endl;
        std::cout << "sigma_xx: " << analyticStressX << std::endl;
        std::cout << "force_x : " << analyticForce << std::endl;

        std::cout << rS.ElementGetEngineeringStrain(0) << std::endl;

        // get element stresses
        int allElements = rS.GroupCreate(NuTo::eGroupId::Elements);
        int allNodes    = rS.GroupCreate(NuTo::eGroupId::Nodes);
        rS.GroupAddNodeCoordinateRange(allNodes, 0, -0.1, lX+0.1);
        rS.GroupAddElementsFromNodes(allElements, allNodes, true);
        for (int elementId : rS.GroupGetMemberIds(allElements))
        {
            auto stress = rS.ElementGetEngineeringStress(elementId);
            for (int iIP = 0; iIP < stress.cols(); ++iIP)
            {
                double numericStress = stress(0,iIP);
                std::cout << "numeric stress in element " << elementId << " at IP " << iIP << ": " << numericStress << std::endl;

                BOOST_CHECK_CLOSE_FRACTION(numericStress, analyticStressX, 1.e-6);
            }
        }


        // sum up reaction forces in x direction of all nodes on x = 0
        double numericForce = 0;
        int nodesX0 = rS.GroupCreate(NuTo::eGroupId::Nodes);
        rS.GroupAddNodeCoordinateRange(nodesX0, 0, lX-1.e-6, lX+1.e-6);
        for (int nodeId : rS.GroupGetMemberIds(nodesX0))
        {
            Eigen::VectorXd force;
            rS.NodeInternalForce(nodeId, force);
            numericForce += force(0);
        }

        BOOST_CHECK_CLOSE(numericForce, analyticForce, 1.e-6);
    }

    void CheckMass(NuTo::Structure& rS)
    {
        double analyticMass = lX*lY*lZ*rho;

        auto hessian2 = rS.BuildGlobalHessian2();

        double numericMass = 0;
        numericMass += hessian2.JJ(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS).Sum();
        numericMass += hessian2.JK(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS).Sum();
        numericMass += hessian2.KJ(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS).Sum();
        numericMass += hessian2.KK(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS).Sum();

        numericMass /= rS.GetDimension(); // since the mass is added to nodes in every direction

        BOOST_CHECK_CLOSE(numericMass, analyticMass, 1.e-6);

        hessian2 = rS.BuildGlobalHessian2Lumped();

        numericMass = 0;
        numericMass += hessian2.JJ(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS).Sum();
        numericMass += hessian2.JK(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS).Sum();
        numericMass += hessian2.KJ(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS).Sum();
        numericMass += hessian2.KK(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS).Sum();

        numericMass /= rS.GetDimension(); // since the mass is added to nodes in every direction

        BOOST_CHECK_CLOSE(numericMass, analyticMass, 1.e-6);
    }

    void Visualize(NuTo::Structure& rS, const std::string& rVisualizationDirectory)
    {
#ifdef ENABLE_VISUALIZE
        if (rVisualizationDirectory == "")
            throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Provide a valid visualization directory!");

        boost::filesystem::path directory(rVisualizationDirectory);
        boost::filesystem::create_directory(directory);

        // get a random element
        const NuTo::ElementBase* element = nullptr;
        int elementID = 0;
        while (element == nullptr)
        {
            try
            {
                element = rS.ElementGetElementPtr(elementID);
            }
            catch (...){}
            elementID++;
        }
        std::string fileName = NuTo::Interpolation::ShapeTypeToString(element->GetInterpolationType().GetShapeType());
        fileName += NuTo::Interpolation::TypeOrderToString(element->GetInterpolationType().Get(NuTo::Node::eDof::DISPLACEMENTS).GetTypeOrder());
        fileName += ".vtu";
        directory /= fileName;

        int visualizationGroup = rS.GroupCreate(NuTo::eGroupId::Elements);
        rS.GroupAddElementsTotal(visualizationGroup);

        rS.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::DISPLACEMENTS);
        rS.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRAIN);
        rS.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRESS);

        rS.ExportVtkDataFileElements(directory.string(),true);

#endif
    }


};

}//namespace NuToTest

