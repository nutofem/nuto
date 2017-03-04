#include "BoostUnitTest.h"

#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/structures/StructureOutputBlockMatrix.h"
#include "mechanics/constitutive/laws/AdditiveInputExplicit.h"
#include "mechanics/constitutive/laws/AdditiveOutput.h"
#include "mechanics/mesh/MeshGenerator.h"
#include "mechanics/MechanicsEnums.h"

using namespace NuTo;

class TrussStructure
{
public:
    TrussStructure()
        : structure(1)
    {
        // create section
        int section = structure.SectionCreate("Truss");
        structure.SectionSetArea(section, 1.0);

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
        for (int node = 0; node < structure.GetNumNodes(); ++node)
        {
            structure.NodeSetTemperature(node, temperature);
        }
        structure.ConstraintLinearSetTemperatureNode(0, temperature);
    }


    int SetConstitutiveLaws()
    {
        int additive_input_id =
                structure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::ADDITIVE_INPUT_EXPLICIT);
        int additive_output_id =
                structure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::ADDITIVE_OUTPUT);

        int lin_elastic_id = structure.ConstitutiveLawCreate(
                NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
        structure.ConstitutiveLawSetParameterDouble(
                lin_elastic_id, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, youngsModulus);

        double capacity     = 1.0;
        double conductivity = 1.0;
        double density      = 1.0;
        int heat_conduction_id =
                structure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::HEAT_CONDUCTION);
        structure.ConstitutiveLawSetParameterDouble(
                heat_conduction_id, NuTo::Constitutive::eConstitutiveParameter::HEAT_CAPACITY, capacity);
        structure.ConstitutiveLawSetParameterDouble(
                heat_conduction_id, NuTo::Constitutive::eConstitutiveParameter::THERMAL_CONDUCTIVITY, conductivity);
        structure.ConstitutiveLawSetParameterDouble(heat_conduction_id,
                                                    NuTo::Constitutive::eConstitutiveParameter::DENSITY, density);

        int thermal_strains_id =
                structure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::THERMAL_STRAINS);
        structure.ConstitutiveLawSetParameterDouble(
                thermal_strains_id, NuTo::Constitutive::eConstitutiveParameter::THERMAL_EXPANSION_COEFFICIENT,
                expansionCoefficient);

        auto additive_input = static_cast<NuTo::AdditiveInputExplicit*>(
                structure.ConstitutiveLawGetConstitutiveLawPtr(additive_input_id));
        auto additive_output =
                static_cast<NuTo::AdditiveOutput*>(structure.ConstitutiveLawGetConstitutiveLawPtr(additive_output_id));
        NuTo::ConstitutiveBase* lin_elastic     = structure.ConstitutiveLawGetConstitutiveLawPtr(lin_elastic_id);
        NuTo::ConstitutiveBase* thermal_strains = structure.ConstitutiveLawGetConstitutiveLawPtr(thermal_strains_id);
        NuTo::ConstitutiveBase* heat_conduction = structure.ConstitutiveLawGetConstitutiveLawPtr(heat_conduction_id);

        additive_input->AddConstitutiveLaw(*lin_elastic);
        additive_input->AddConstitutiveLaw(*thermal_strains, NuTo::Constitutive::eInput::ENGINEERING_STRAIN);

        additive_output->AddConstitutiveLaw(*additive_input);
        additive_output->AddConstitutiveLaw(*heat_conduction);

        return additive_output_id;
    }

    NuTo::Structure structure;
    double temperature          = 50.0;
    double expansionCoefficient = 23.1e-6;
    double youngsModulus        = 20000.;
    double length               = 10.0;
};

BOOST_FIXTURE_TEST_CASE(restrained, TrussStructure)
{
    // set boundary conditions and loads
    auto direction = Eigen::Matrix<double, 1, 1>::Ones();
    structure.ConstraintLinearSetDisplacementNode(0, direction, 0.0);
    structure.ConstraintLinearSetDisplacementNode(structure.GetNumNodes() - 1, direction, 0.0);

    // solve
    structure.SolveGlobalSystemStaticElastic();

    auto lastNode = structure.GetNumNodes() - 1;
    Eigen::VectorXd displacement;
    structure.NodeGetDisplacements(lastNode, displacement);
    BOOST_CHECK_SMALL(displacement[0], 1e-9);

    double analytic_strain = -temperature * expansionCoefficient;
    auto strain            = structure.ElementGetEngineeringStrain(0);
    BOOST_CHECK_CLOSE(strain(0, 0), analytic_strain, 1e-9);

    double analytic_stress = youngsModulus * analytic_strain;
    auto stress            = structure.ElementGetEngineeringStress(0);
    BOOST_CHECK_CLOSE(stress(0, 0), analytic_stress, 1e-9);
}


/// @TODO this doesn't work; need to investigate

//BOOST_FIXTURE_TEST_CASE(unrestrained, TrussStructure)
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
