/*
 * ElementUniaxialTest.h
 *
 *  Created on: 13 May 2015
 *      Author: ttitsche
 */

#pragma once
#include "BoostUnitTest.h"

#include <iomanip>
#include <boost/filesystem.hpp>

#include "mechanics/MechanicsEnums.h"
#include "visualize/VisualizeEnum.h"

#include "math/SparseMatrixCSRVector2.h"
#include "mechanics/dofSubMatrixStorage/BlockScalar.h"
#include "mechanics/constraints/ConstraintCompanion.h"
#include "mechanics/elements/ElementBase.h"
#include "mechanics/groups/Group.h"
#include "mechanics/interpolationtypes/InterpolationBase.h"
#include "mechanics/interpolationtypes/InterpolationType.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/structures/StructureOutputBlockMatrix.h"
#include "mechanics/structures/StructureOutputBlockVector.h"
#include "mechanics/mesh/MeshGenerator.h"

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
    void Run(NuTo::Structure& s, const std::string& rVisualizationDirectory)
    {

        CheckVolume(s);

        s.SetVerboseLevel(10);
        s.SetShowTime(true);


        SetConstitutiveLaw(s);
        SetBoundaryConditions(s);

        s.NodeBuildGlobalDofs();
        s.CalculateMaximumIndependentSets();

        s.NodeInfo(10);

        CheckStiffness(s);
        Solve(s);
        CheckSolution(s);
        CheckMass(s);
        Visualize(s, rVisualizationDirectory);
    }

private:

    void CheckVolume(NuTo::Structure& s)
    {
        double volume = s.ElementGroupGetVolume(s.GroupGetElementsTotal());
        if (s.GetDimension() == 2)
            volume *= lZ;
        if (s.GetDimension() == 1)
            volume *= lY * lZ;

        BOOST_CHECK_CLOSE(volume, lX * lY * lZ, 1.e-6);
    }

    void SetConstitutiveLaw(NuTo::Structure& s)
    {
        using namespace NuTo::Constitutive;
        s.ConstitutiveLawCreate(0, eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
        s.ConstitutiveLawSetParameterDouble(0, eConstitutiveParameter::YOUNGS_MODULUS, E);
        s.ConstitutiveLawSetParameterDouble(0, eConstitutiveParameter::POISSONS_RATIO,nu);
        s.ConstitutiveLawSetParameterDouble(0, eConstitutiveParameter::DENSITY, rho);
        s.ElementTotalSetConstitutiveLaw(0);
    }

    void SetBoundaryConditions(NuTo::Structure& s)
    {
        int dimension = s.GetDimension();

        const NuTo::Node::eDof eDofDispl = NuTo::Node::eDof::DISPLACEMENTS;

        // fix origin
        const auto& nodeOrigin = s.NodeGetAtCoordinate(Eigen::VectorXd::Zero(dimension));
        for (int iDim = 1; iDim < dimension; ++iDim)
            s.Constraints().Add(eDofDispl, NuTo::Constraint::Component(nodeOrigin, {static_cast<NuTo::eDirection>(iDim)}));

        if (dimension == 3)
        {
        // fix rotation
            const auto& nodeRotation = s.NodeGetAtCoordinate(Eigen::Vector3d(0,0,lZ));
            s.Constraints().Add(eDofDispl, NuTo::Constraint::Component(nodeRotation, {NuTo::eDirection::Y}));
        }

        // fix x = 0 plane
        const auto& groupX0 = s.GroupGetNodesAtCoordinate(NuTo::eDirection::X, 0);
        s.Constraints().Add(eDofDispl, NuTo::Constraint::Component(groupX0, {NuTo::eDirection::X}));

        // apply displacement on x = lX plane
        const auto& groupXlX = s.GroupGetNodesAtCoordinate(NuTo::eDirection::X, lX);
        s.Constraints().Add(eDofDispl, NuTo::Constraint::Component(groupXlX, {NuTo::eDirection::X}, deltaL));
    }

    void CheckStiffness(NuTo::Structure& s)
    {
        if (s.GetNumTotalDofs() > 50)
            return;
        BOOST_CHECK(s.ElementCheckHessian0(1.e-8, 1.e-6, true));
        BOOST_CHECK(s.CheckHessian0(1.e-8, 1.e-6, true));
    }

    void Solve(NuTo::Structure& s)
    {
        s.SolveGlobalSystemStaticElastic();
        auto internalGradient = s.BuildGlobalInternalGradient();
        BOOST_CHECK_SMALL(internalGradient.J.CalculateNormL2()[NuTo::Node::eDof::DISPLACEMENTS], 1e-8);
    }

    void CheckSolution(NuTo::Structure& s)
    {
        double analyticStrainX = deltaL / lX;
        double analyticStressX = analyticStrainX * E;
        double analyticForce = analyticStressX * lY * lZ;

        std::cout << "\n  ### analytical values ###" << std::endl;
        std::cout << "sigma_xx: " << analyticStressX << std::endl;
        std::cout << "force_x : " << analyticForce << std::endl;

        // get element stresses
        int allElements = s.GroupGetElementsTotal();
        for (int elementId : s.GroupGetMemberIds(allElements))
        {
            auto stress = s.ElementGetEngineeringStress(elementId);
            for (int iIP = 0; iIP < stress.cols(); ++iIP)
            {
                double numericStress = stress(0,iIP);
                std::cout << "numeric stress in element " << elementId << " at IP " << iIP << ": " << numericStress << std::endl;

                BOOST_CHECK_CLOSE_FRACTION(numericStress, analyticStressX, 1.e-6);
            }
        }


        // sum up reaction forces in x direction of all nodes on x = 0
        double numericForce = 0;
        const auto& nodesX0 = s.GroupGetNodesAtCoordinate(NuTo::eDirection::X, lX);
        for (int nodeId : nodesX0.GetMemberIds())
        {
            Eigen::VectorXd force;
            s.NodeInternalForce(nodeId, force);
            numericForce += force(0);
        }

        BOOST_CHECK_CLOSE(numericForce, analyticForce, 1.e-6);
    }

    void CheckMass(NuTo::Structure& s)
    {
        double analyticMass = lX*lY*lZ*rho;

        auto hessian2 = s.BuildGlobalHessian2();

        double numericMass = 0;
        numericMass += hessian2.JJ(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS).Sum();
        numericMass += hessian2.JK(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS).Sum();
        numericMass += hessian2.KJ(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS).Sum();
        numericMass += hessian2.KK(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS).Sum();

        numericMass /= s.GetDimension(); // since the mass is added to nodes in every direction

        BOOST_CHECK_CLOSE(numericMass, analyticMass, 1.e-6);

        hessian2 = s.BuildGlobalHessian2Lumped();

        numericMass = 0;
        numericMass += hessian2.JJ(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS).Sum();
        numericMass += hessian2.JK(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS).Sum();
        numericMass += hessian2.KJ(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS).Sum();
        numericMass += hessian2.KK(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS).Sum();

        numericMass /= s.GetDimension(); // since the mass is added to nodes in every direction

        BOOST_CHECK_CLOSE(numericMass, analyticMass, 1.e-6);
    }

    void Visualize(NuTo::Structure& s, const std::string& rVisualizationDirectory)
    {
#ifdef ENABLE_VISUALIZE
        if (rVisualizationDirectory == "")
            throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Provide a valid visualization directory!");

        boost::filesystem::path directory(rVisualizationDirectory);
        boost::filesystem::create_directory(directory);

        const auto& interpolationType = *s.InterpolationTypeGet(0);
        std::string fileName = NuTo::Interpolation::ShapeTypeToString(interpolationType.GetShapeType());
        fileName += NuTo::Interpolation::TypeOrderToString(interpolationType.Get(NuTo::Node::eDof::DISPLACEMENTS).GetTypeOrder());
        fileName += ".vtu";
        directory /= fileName;

        int visualizationGroup = s.GroupCreate(NuTo::eGroupId::Elements);
        s.GroupAddElementsTotal(visualizationGroup);

        s.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::DISPLACEMENTS);
        s.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRAIN);
        s.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRESS);

        s.ExportVtkDataFileElements(directory.string(),true);
#endif
    }
};
}//namespace NuToTest

