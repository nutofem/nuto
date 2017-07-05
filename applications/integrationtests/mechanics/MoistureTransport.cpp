#include "base/Timer.h"
#include "mechanics/MechanicsException.h"
#include "mechanics/constitutive/staticData/DataMoistureTransport.h"
#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "mechanics/timeIntegration/NewmarkDirect.h"
#include "mechanics/mesh/MeshGenerator.h"

#include <array>


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  TESTS
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#include "MoistureTransport_ConstitutiveOutputTest.h"
#include "MoistureTransport_SimulationTest.h"

void RunSimulationTest_AllDiemnsions(std::map<NuTo::Node::eDof, NuTo::Interpolation::eTypeOrder> rDofIPTMap,
                                     int rNumElementsXYZ)
{
    SimulationTest<1>({rNumElementsXYZ}, {0.01}, {true}, rDofIPTMap, true);
    SimulationTest<2>({rNumElementsXYZ, rNumElementsXYZ}, {0.01, 0.01}, {true, true}, rDofIPTMap, true);
    SimulationTest<3>({rNumElementsXYZ, rNumElementsXYZ, rNumElementsXYZ}, {0.01, 0.01, 0.01}, {true, true, true},
                      rDofIPTMap, true);
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  Main
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int main()
{
    std::map<NuTo::Node::eDof, NuTo::Interpolation::eTypeOrder> dofIPTMap;
    dofIPTMap[NuTo::Node::eDof::RELATIVEHUMIDITY] = NuTo::Interpolation::eTypeOrder::EQUIDISTANT1;
    dofIPTMap[NuTo::Node::eDof::WATERVOLUMEFRACTION] = NuTo::Interpolation::eTypeOrder::EQUIDISTANT1;

    // Constitutive output tests

    ConstitutiveOutputTests<1>(dofIPTMap);
    ConstitutiveOutputTests<2>(dofIPTMap);
    ConstitutiveOutputTests<3>(dofIPTMap);

    // Simulation Tests --- Compare results to paper

    SimulationTest<1>({16}, {0.16}, {true}, dofIPTMap);

    SimulationTest<2>({16, 1}, {0.16, 0.1}, {true, false}, dofIPTMap);

    SimulationTest<2>({1, 16}, {0.1, 0.16}, {false, true}, dofIPTMap);


    SimulationTest<3>({16, 1, 1}, {0.16, 0.1, 0.1}, {true, false, false}, dofIPTMap);

    SimulationTest<3>({1, 16, 1}, {0.1, 0.16, 0.1}, {false, true, false}, dofIPTMap);

    SimulationTest<3>({1, 1, 16}, {0.1, 0.1, 0.16}, {false, false, true}, dofIPTMap);

    // Simulation Tests --- Interpolation type


    dofIPTMap[NuTo::Node::eDof::RELATIVEHUMIDITY] = NuTo::Interpolation::eTypeOrder::EQUIDISTANT2;
    dofIPTMap[NuTo::Node::eDof::WATERVOLUMEFRACTION] = NuTo::Interpolation::eTypeOrder::EQUIDISTANT1;
    RunSimulationTest_AllDiemnsions(dofIPTMap, 2);


    dofIPTMap[NuTo::Node::eDof::RELATIVEHUMIDITY] = NuTo::Interpolation::eTypeOrder::EQUIDISTANT1;
    dofIPTMap[NuTo::Node::eDof::WATERVOLUMEFRACTION] = NuTo::Interpolation::eTypeOrder::EQUIDISTANT2;
    RunSimulationTest_AllDiemnsions(dofIPTMap, 2);


    dofIPTMap[NuTo::Node::eDof::RELATIVEHUMIDITY] = NuTo::Interpolation::eTypeOrder::EQUIDISTANT2;
    dofIPTMap[NuTo::Node::eDof::WATERVOLUMEFRACTION] = NuTo::Interpolation::eTypeOrder::EQUIDISTANT2;
    RunSimulationTest_AllDiemnsions(dofIPTMap, 2);
    return 0;
}
