/*
 * ElementUniaxialTest.h
 *
 *  Created on: 13 May 2015
 *      Author: ttitsche
 */

#pragma once

#include "math/SparseMatrixCSRVector2.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/dofSubMatrixStorage/BlockScalar.h"
#include "mechanics/elements/ElementBase.h"
#include "mechanics/groups/GroupEnum.h"
#include "mechanics/sections/SectionEnum.h"
#include "mechanics/interpolationtypes/InterpolationBase.h"
#include "mechanics/interpolationtypes/InterpolationType.h"
#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/structures/StructureOutputBlockMatrix.h"
#include "mechanics/structures/StructureOutputBlockVector.h"
#include "math/SparseDirectSolverMUMPS.h"
#include <boost/filesystem.hpp>
#include "math/FullMatrix.h"

#ifdef ENABLE_VISUALIZE
#include "visualize/VisualizeEnum.h"
#endif

#define DEBUG_PRINT true

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

    std::string visualizationDirectory;

    //! @brief provide a structure that contains the geometry
    //! of a uniaxial test, starting from (0,0,0) to (lx, ly, lz)
    void Run(NuTo::Structure& rStructure)
    {
        rStructure.SetVerboseLevel(10);

        SetConstitutiveLaw(rStructure);
        SetBoundaryConditions(rStructure);

        rStructure.NodeBuildGlobalDofs();
        rStructure.CalculateMaximumIndependentSets();

        if (DEBUG_PRINT)
            rStructure.NodeInfo(10);


        CheckStiffness(rStructure);
        Solve(rStructure);
        CheckSolution(rStructure);
        CheckMass(rStructure);
        Visualize(rStructure);
    }

private:

    void SetConstitutiveLaw(NuTo::Structure& rStructure)
    {
        rStructure.ConstitutiveLawCreate(0, NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
        rStructure.ConstitutiveLawSetParameterDouble(0,NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, E);
        rStructure.ConstitutiveLawSetParameterDouble(0,NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO,nu);
        rStructure.ConstitutiveLawSetParameterDouble(0,NuTo::Constitutive::eConstitutiveParameter::DENSITY, rho);
        rStructure.ElementTotalSetConstitutiveLaw(0);
    }

    void SetBoundaryConditions(NuTo::Structure& rStructure)
    {

        int dimension = rStructure.GetDimension();

        NuTo::FullMatrix<double, Eigen::Dynamic> directions = NuTo::FullMatrix<double, Eigen::Dynamic>::Identity(dimension, dimension);

        // fix origin
        NuTo::FullVector<double, Eigen::Dynamic> origin(dimension);
        origin.setZero();
        int nodeGroupOrigin = rStructure.GroupCreate("Nodes");
        rStructure.GroupAddNodeRadiusRange(nodeGroupOrigin, origin, 0, 1.e-5);

        if (rStructure.GroupGetNumMembers(nodeGroupOrigin) != 1)
            throw NuTo::MechanicsException("[NuToTest::ElementUniaxialTest::SetBoundaryConditions] Node at origin (0,0,0) does not exist.");

        int nodeOrigin = rStructure.GroupGetMemberIds(nodeGroupOrigin)[0];
        for (int iDim = 1; iDim < dimension; ++iDim)
            rStructure.ConstraintLinearSetDisplacementNode(nodeOrigin, directions.GetColumn(iDim), 0.0);

        if (dimension == 3)
        {
        // fix rotation
            NuTo::FullVector<double, Eigen::Dynamic> originRot(3);
            originRot << 0, 0, lZ;
            int nodeGroupOriginRot = rStructure.GroupCreate(NuTo::eGroupId::Nodes);
            rStructure.GroupAddNodeRadiusRange(nodeGroupOriginRot, originRot, 0, 1.e-5);

            int nodeOriginRot = rStructure.GroupGetMemberIds(nodeGroupOriginRot)[0];
            rStructure.ConstraintLinearSetDisplacementNode(nodeOriginRot, directions.GetColumn(1), 0.0);
        }

        // fix x = 0 plane
        int nodesX0 = rStructure.GroupCreate(NuTo::eGroupId::Nodes);
        rStructure.GroupAddNodeCoordinateRange(nodesX0, 0, -1.e-6, 1.e-6);
        rStructure.ConstraintLinearSetDisplacementNodeGroup(nodesX0, directions.GetColumn(0), 0.);

        // apply displacement on x = lX plane
        int nodesXlX = rStructure.GroupCreate(NuTo::eGroupId::Nodes);
        rStructure.GroupAddNodeCoordinateRange(nodesXlX, 0, lX-1.e-6, lX+1.e-6);
        rStructure.ConstraintLinearSetDisplacementNodeGroup(nodesXlX, directions.GetColumn(0), deltaL);

        if (DEBUG_PRINT)
        {
            std::cout << "Number of nodes at origin :" << rStructure.GroupGetNumMembers(nodeGroupOrigin) << std::endl;
            std::cout << "Number of nodes at x = 0  :" << rStructure.GroupGetNumMembers(nodesX0)<< std::endl;
            std::cout << "Number of nodes at x = lX :" << rStructure.GroupGetNumMembers(nodesXlX)<< std::endl;
        }

    }
    void CheckStiffness(NuTo::Structure& rStructure)
    {
        if (rStructure.GetNumTotalDofs() > 50)
            return;
        rStructure.ElementCheckHessian0(1.e-8, 1.e-6, true);
        rStructure.CheckHessian0(1.e-8, 1.e-6, true);
    }

    void Solve(NuTo::Structure& rStructure)
    {
        rStructure.SolveGlobalSystemStaticElastic();

        auto internalGradient = rStructure.BuildGlobalInternalGradient();

        if (internalGradient.J.CalculateNormL2()[NuTo::Node::eDof::DISPLACEMENTS] > 1e-8)
        {
            std::cout << internalGradient.J[NuTo::Node::eDof::DISPLACEMENTS];
            throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "residual force vector is not zero.");
        }
    }

    void CheckSolution(NuTo::Structure& rStructure)
    {
        double analyticStrainX = deltaL / lX;
        double analyticStressX = analyticStrainX * E;
        double analyticForce = analyticStressX * lY * lZ;
        if (DEBUG_PRINT)
        {
            std::cout << "\n  ### analytical values ###" << std::endl;
            std::cout << "sigma_xx: " << analyticStressX << std::endl;
            std::cout << "force_x : " << analyticForce << std::endl;
        }

        std::cout << rStructure.ElementGetEngineeringStrain(0) << std::endl;

        // get element stresses
        int allElements = rStructure.GroupCreate(NuTo::eGroupId::Elements);
        int allNodes    = rStructure.GroupCreate(NuTo::eGroupId::Nodes);
        rStructure.GroupAddNodeCoordinateRange(allNodes, 0, -0.1, lX+0.1);
        rStructure.GroupAddElementsFromNodes(allElements, allNodes, true);
        NuTo::FullVector<int, Eigen::Dynamic> elementIds = rStructure.GroupGetMemberIds(allElements);
        for (int iElement = 0; iElement < elementIds.GetNumRows(); ++iElement)
        {
            int elementId = elementIds(iElement);
            auto stress = rStructure.ElementGetEngineeringStress(elementId);
            for (int iIP = 0; iIP < stress.cols(); ++iIP)
            {
                double numericStress = stress(0,iIP);
                if (DEBUG_PRINT)
                    std::cout << "numeric stress in element " << elementId << " at IP " << iIP << ": " << numericStress << std::endl;
                if(std::abs(numericStress - analyticStressX) > 1.e-6 )
                {
                    std::cout << "sigma_xx analytical : " << analyticStressX << std::endl;
                    std::cout << "sigma_xx numerical  : " << numericStress << std::endl;
                    throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "wrong stress calculation.");
                }
            }
        }


        // sum up reaction forces in x direction of all nodes on x = 0
        double numericForce = 0;
        int nodesX0 = rStructure.GroupCreate(NuTo::eGroupId::Nodes);
        rStructure.GroupAddNodeCoordinateRange(nodesX0, 0, lX-1.e-6, lX+1.e-6);
        NuTo::FullVector<int, Eigen::Dynamic> nodeX0Indices = rStructure.GroupGetMemberIds(nodesX0);
        for (int iNodeX0 = 0; iNodeX0 < nodeX0Indices.GetNumRows(); ++iNodeX0)
        {
            int nodeX0Id = nodeX0Indices(iNodeX0);
            NuTo::FullVector<double, Eigen::Dynamic> force;
            rStructure.NodeInternalForce(nodeX0Id, force);
            numericForce += force(0);
        }

        if(std::abs(numericForce - analyticForce) > 1.e-6 )
        {
            std::cout << "force_x analytical : " << analyticForce << std::endl;
            std::cout << "force_x numerical  : " << numericForce << std::endl;
            throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "wrong reaction force calculation.");
        }

    }

    void CheckMass(NuTo::Structure& rStructure)
    {
        double analyticMass = lX*lY*lZ*rho;

        auto hessian2 = rStructure.BuildGlobalHessian2();

        double numericMass = 0;
        numericMass += hessian2.JJ(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS).Sum();
        numericMass += hessian2.JK(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS).Sum();
        numericMass += hessian2.KJ(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS).Sum();
        numericMass += hessian2.KK(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS).Sum();

        numericMass /= rStructure.GetDimension(); // since the mass is added to nodes in every direction

        if(std::abs(numericMass - analyticMass)/numericMass > 1.e-6 )
        {
            std::cout << "mass analytical : " << std::setprecision(10) << analyticMass << std::endl;
            std::cout << "mass numerical  : " << std::setprecision(10) << numericMass << std::endl;
            throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "wrong mass calculation.");
        }


        hessian2 = rStructure.BuildGlobalHessian2Lumped();

        numericMass = 0;
        numericMass += hessian2.JJ(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS).Sum();
        numericMass += hessian2.JK(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS).Sum();
        numericMass += hessian2.KJ(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS).Sum();
        numericMass += hessian2.KK(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS).Sum();

        numericMass /= rStructure.GetDimension(); // since the mass is added to nodes in every direction

        if(std::abs(numericMass - analyticMass)/numericMass > 1.e-6 )
        {
            std::cout << "mass analytical : " << std::setprecision(10) << analyticMass << std::endl;
            std::cout << "mass numerical  : " << std::setprecision(10) << numericMass << std::endl;
            throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "wrong lumped mass calculation.");
        }

    }

    void Visualize(NuTo::Structure& rStructure)
    {
#ifdef ENABLE_VISUALIZE
        if (visualizationDirectory == "")
            throw NuTo::MechanicsException("[NuToTest::ElementUniaxialTest::Visualize] Set a visualization directory (this->visualizationDirectory) first!");

        boost::filesystem::path directory(visualizationDirectory);
        directory /= "ResultsElementUniaxialTest";
        boost::filesystem::create_directory(directory);

        // get a random element
        const NuTo::ElementBase* element = nullptr;
        int elementID = 0;
        while (element == nullptr)
        {
            try
            {
                element = rStructure.ElementGetElementPtr(elementID);
            }
            catch (...){}
            elementID++;
        }
        std::string fileName = NuTo::Interpolation::ShapeTypeToString(element->GetInterpolationType().GetShapeType());
        fileName += NuTo::Interpolation::TypeOrderToString(element->GetInterpolationType().Get(NuTo::Node::eDof::DISPLACEMENTS).GetTypeOrder());
        fileName += ".vtu";
        directory /= fileName;

        int visualizationGroup = rStructure.GroupCreate(NuTo::eGroupId::Elements);
        rStructure.GroupAddElementsTotal(visualizationGroup);

        rStructure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::DISPLACEMENTS);
        rStructure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRAIN);
        rStructure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRESS);

        rStructure.ExportVtkDataFileElements(directory.string(),true);

#endif
    }

    double E = 42;
    double nu = 0.26;
    double deltaL = 0.6174;
    double rho = M_PI;

};

}//namespace NuToTest

