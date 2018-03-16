#include "BoostUnitTest.h"

#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/structures/StructureOutputBlockMatrix.h"
#include "nuto/mechanics/constitutive/laws/AdditiveInputExplicit.h"
#include "nuto/mechanics/constitutive/laws/AdditiveOutput.h"
#include "nuto/mechanics/mesh/MeshGenerator.h"
#include "nuto/mechanics/MechanicsEnums.h"
#include "nuto/mechanics/sections/SectionTruss.h"
#include "nuto/mechanics/constraints/ConstraintCompanion.h"

using namespace NuTo;

class TrussStructure
{
public:
    TrussStructure()
        : structure(1)
    {
        // create section
        auto section = SectionTruss::Create(1.0);

        auto additive_output = SetConstitutiveLaws();

        // generate mesh
        int group, interpolationType;
        std::tie(group, interpolationType) =
                MeshGenerator::Grid(structure, {length}, {1}, NuTo::Interpolation::eShapeType::TRUSS1D);

        structure.InterpolationTypeAdd(interpolationType, NuTo::Node::eDof::DISPLACEMENTS,
                                       NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
        structure.InterpolationTypeAdd(interpolationType, NuTo::Node::eDof::TEMPERATURE,
                                       NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);

        structure.ElementTotalSetSection(section);
        structure.ElementTotalSetConstitutiveLaw(additive_output);

        structure.InterpolationTypeSetIntegrationType(interpolationType,
                                                      NuTo::eIntegrationType::IntegrationType1D2NGauss1Ip);

        structure.ElementTotalConvertToInterpolationType();

        // set uniform temperature
        int gAllNodes = structure.GroupGetNodesTotal();
        for (int nodeId : structure.GroupGetMemberIds(gAllNodes))
            structure.NodeSetTemperature(nodeId, temperature);

        structure.Constraints().Add(NuTo::Node::eDof::TEMPERATURE,
                                    NuTo::Constraint::Value(*structure.NodeGetNodePtr(0), temperature));
    }


    int SetConstitutiveLaws()
    {
        using namespace NuTo::Constitutive;
        int additive_input_id = structure.ConstitutiveLawCreate(eConstitutiveType::ADDITIVE_INPUT_EXPLICIT);
        int additive_output_id = structure.ConstitutiveLawCreate(eConstitutiveType::ADDITIVE_OUTPUT);

        int lin_elastic_id = structure.ConstitutiveLawCreate(eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
        structure.ConstitutiveLawSetParameterDouble(lin_elastic_id, eConstitutiveParameter::YOUNGS_MODULUS,
                                                    youngsModulus);

        double capacity = 1.0;
        double conductivity = 1.0;
        double density = 1.0;
        int heat_conduction_id = structure.ConstitutiveLawCreate(eConstitutiveType::HEAT_CONDUCTION);
        structure.ConstitutiveLawSetParameterDouble(heat_conduction_id, eConstitutiveParameter::HEAT_CAPACITY,
                                                    capacity);
        structure.ConstitutiveLawSetParameterDouble(heat_conduction_id, eConstitutiveParameter::THERMAL_CONDUCTIVITY,
                                                    conductivity);
        structure.ConstitutiveLawSetParameterDouble(heat_conduction_id, eConstitutiveParameter::DENSITY, density);

        int thermal_strains_id = structure.ConstitutiveLawCreate(eConstitutiveType::THERMAL_STRAINS);
        structure.ConstitutiveLawSetParameterDouble(
                thermal_strains_id, eConstitutiveParameter::THERMAL_EXPANSION_COEFFICIENT, expansionCoefficient);

        auto additive_input = static_cast<NuTo::AdditiveInputExplicit*>(
                structure.ConstitutiveLawGetConstitutiveLawPtr(additive_input_id));
        auto additive_output =
                static_cast<NuTo::AdditiveOutput*>(structure.ConstitutiveLawGetConstitutiveLawPtr(additive_output_id));
        NuTo::ConstitutiveBase* lin_elastic = structure.ConstitutiveLawGetConstitutiveLawPtr(lin_elastic_id);
        NuTo::ConstitutiveBase* thermal_strains = structure.ConstitutiveLawGetConstitutiveLawPtr(thermal_strains_id);
        NuTo::ConstitutiveBase* heat_conduction = structure.ConstitutiveLawGetConstitutiveLawPtr(heat_conduction_id);

        additive_input->AddConstitutiveLaw(*lin_elastic);
        additive_input->AddConstitutiveLaw(*thermal_strains, eInput::ENGINEERING_STRAIN);

        additive_output->AddConstitutiveLaw(*additive_input);
        additive_output->AddConstitutiveLaw(*heat_conduction);

        return additive_output_id;
    }

    NuTo::Structure structure;
    double temperature = 50.0;
    double expansionCoefficient = 23.1e-6;
    double youngsModulus = 20000.;
    double length = 10.0;
};

BOOST_FIXTURE_TEST_CASE(restrained, TrussStructure)
{
    // set boundary conditions and loads
    auto& nodeLeft = structure.NodeGetAtCoordinate(0);
    auto& nodeRight = structure.NodeGetAtCoordinate(length);
    structure.Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Constraint::Value(nodeLeft));
    structure.Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Constraint::Value(nodeRight));

    structure.SolveGlobalSystemStaticElastic();

    BOOST_CHECK_SMALL(nodeRight.Get(NuTo::Node::eDof::DISPLACEMENTS)[0], 1e-9);

    double analytic_strain = -temperature * expansionCoefficient;
    auto strain = structure.ElementGetEngineeringStrain(0);
    BOOST_CHECK_CLOSE(strain(0, 0), analytic_strain, 1e-9);

    double analytic_stress = youngsModulus * analytic_strain;
    auto stress = structure.ElementGetEngineeringStress(0);
    BOOST_CHECK_CLOSE(stress(0, 0), analytic_stress, 1e-9);
}


/// @TODO this doesn't work; need to investigate

// BOOST_FIXTURE_TEST_CASE(unrestrained, TrussStructure)
//{
//    // set boundary conditions and loads
//    auto direction = Eigen::Matrix<double, 1, 1>::Ones();
//    structure.ConstraintLinearSetDisplacementNode(0, direction, 0.0);
//
//    // solve
//    structure.SolveGlobalSystemStaticElastic();
//    std::cout << structure.BuildGlobalHessian0() << '\n';
//    std::cout << structure.BuildGlobalExternalLoadVector() << '\n';
//    std::cout << structure.BuildGlobalInternalGradient() << '\n';
//
//    Eigen::VectorXd displacement;
//    for (int i = 0; i < structure.GetNumNodes(); ++i)
//    {
//        structure.NodeGetDisplacements(i, displacement);
//        std::cout << displacement << std::endl;
//    }
//    auto lastNode = structure.GetNumNodes() - 1;
//    structure.NodeGetDisplacements(lastNode, displacement);
//    double analytic_displacement = length * expansionCoefficient * temperature;
//    BOOST_CHECK_CLOSE(displacement[0], analytic_displacement, 1e-9);
//
//    auto strain = structure.ElementGetEngineeringStrain(0);
//    BOOST_CHECK_CLOSE(strain(0, 0), 0.0, 1e-9);
//
//    auto stress = structure.ElementGetEngineeringStress(0);
//    BOOST_CHECK_CLOSE(stress(0, 0), 0.0, 1e-9);
//}
