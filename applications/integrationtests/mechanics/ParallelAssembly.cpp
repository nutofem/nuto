/**
 * Test file to
 *  1) check the parallel assembly
 *  2) ensure a significant (*) speedup
 *
 *  (*) ... not sure yet. First idea: speedup = 0.8 * numProc
 *
 */
#include "BoostUnitTest.h"

#include "base/Timer.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/MechanicsEnums.h"
#include "mechanics/structures/StructureOutputBlockMatrix.h"
#include "mechanics/structures/StructureOutputBlockVector.h"
#include "math/SparseMatrixCSRVector2.h"
#include "mechanics/mesh/MeshGenerator.h"

void SetupStructure(NuTo::Structure& rStructure, int rNumElementsPerDimension)
{
    rStructure.SetShowTime(false);
    rStructure.SetVerboseLevel(0);
    // create nodes

    double lx = 2, ly = 3, lz = 4.;

    int interpolationType = NuTo::MeshGenerator::Grid(
        rStructure, {lx, ly, lz},
        {rNumElementsPerDimension, rNumElementsPerDimension, rNumElementsPerDimension}).second;

    rStructure.InterpolationTypeAdd(interpolationType,
                                    NuTo::Node::eDof::DISPLACEMENTS,
                                    NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);

    rStructure.ElementTotalConvertToInterpolationType();
    rStructure.ElementTotalSetSection(rStructure.SectionCreate(NuTo::eSectionType::VOLUME));

    rStructure.ConstitutiveLawCreate(0, NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    rStructure.ConstitutiveLawSetParameterDouble(0, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, 10024);
    rStructure.ConstitutiveLawSetParameterDouble(0, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, .12);
    rStructure.ConstitutiveLawSetParameterDouble(0, NuTo::Constitutive::eConstitutiveParameter::DENSITY, 125);
    rStructure.ElementTotalSetConstitutiveLaw(0);
}

void CompareHessiansAndInternalGradients(NuTo::Structure& rStructure1, NuTo::Structure& rStructure2)
{
    auto hessian1 = rStructure1.BuildGlobalHessian0();
    auto hessian2 = rStructure2.BuildGlobalHessian0();

    auto intGrad1 = rStructure1.BuildGlobalInternalGradient();
    auto intGrad2 = rStructure2.BuildGlobalInternalGradient();

    hessian1.AddScal(hessian2, -1);
    intGrad1 -= intGrad2;

    BOOST_CHECK_SMALL(hessian1.JJ(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS).AbsMax(), 1.e-10);
    BOOST_CHECK_SMALL(intGrad1.J[NuTo::Node::eDof::DISPLACEMENTS].norm(), 1.e-10);
}

double MeasureHessianAndInternalGradient(NuTo::Structure& rStructure)
{
    NuTo::Timer timer("", false);
    auto hessian = rStructure.BuildGlobalHessian0();
    auto intGrad = rStructure.BuildGlobalInternalGradient();
    return timer.GetTimeDifference();
}

constexpr int numProc = 4;

BOOST_AUTO_TEST_CASE(ParallelAssemblyCorrectnes)
{
    NuTo::Structure smallStructureSerial(3);
    NuTo::Structure smallStructureParallel(3);

    smallStructureSerial.SetNumProcessors(1);
    smallStructureParallel.SetNumProcessors(numProc);

    SetupStructure(smallStructureSerial, 2);
    SetupStructure(smallStructureParallel, 2);

    smallStructureSerial.CalculateMaximumIndependentSets();
    smallStructureParallel.CalculateMaximumIndependentSets();

    CompareHessiansAndInternalGradients(smallStructureSerial, smallStructureParallel);
}

BOOST_AUTO_TEST_CASE(ParallelAssemblyPerformance)
{
    NuTo::Structure bigStructureSerial(3);
    NuTo::Structure bigStructureParallel(3);

    bigStructureSerial.SetNumProcessors(1);
    bigStructureParallel.SetNumProcessors(numProc);

    SetupStructure(bigStructureSerial, 15);
    SetupStructure(bigStructureParallel, 15);

    bigStructureSerial.CalculateMaximumIndependentSets();
    bigStructureParallel.CalculateMaximumIndependentSets();

    double timeSerial = MeasureHessianAndInternalGradient(bigStructureSerial);
    std::cout << "Serial:   " << timeSerial << std::endl;

    double timeParallel = MeasureHessianAndInternalGradient(bigStructureParallel);
    std::cout << "Parallel: " << timeParallel << std::endl;

    double speedup = timeSerial / timeParallel;
    std::cout << "Speedup:  " << speedup << std::endl;

    constexpr double requiredSpeedup = 0.8 * numProc;
    std::cout << "required: " << requiredSpeedup << std::endl;

    BOOST_CHECK_GT(speedup, requiredSpeedup);
}
