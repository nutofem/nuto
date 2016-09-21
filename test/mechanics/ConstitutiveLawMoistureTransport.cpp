#include "nuto/base/Timer.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/constitutive/staticData/DataMoistureTransport.h"
#include "nuto/mechanics/timeIntegration/NewmarkDirect.h"
#include "nuto/mechanics/tools/MeshGenerator.h"

#include <iostream>
#include <array>



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  TESTS
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#include "ConstitutiveLawMoistureTransport_ConstitutiveOutputTest.h"
#include "ConstitutiveLawMoistureTransport_SimulationTest.h"

void RunSimulationTest_AllDiemnsions(std::map<NuTo::Node::eDof,NuTo::Interpolation::eTypeOrder> rDofIPTMap,
                                     int rNumElementsXYZ)
{
    SimulationTest<1>({rNumElementsXYZ},
                      {0.01},
                      {true},
                      rDofIPTMap,
                      true);
    SimulationTest<2>({rNumElementsXYZ,rNumElementsXYZ},
                      {0.01,0.01},
                      {true,true},
                      rDofIPTMap,
                      true);
    SimulationTest<3>({rNumElementsXYZ,rNumElementsXYZ,rNumElementsXYZ},
                      {0.01,0.01,0.01},
                      {true,true,true},
                      rDofIPTMap,
                      true);
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  Main
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int main()
{
    try
    {
        std::map<NuTo::Node::eDof,NuTo::Interpolation::eTypeOrder> dofIPTMap;
        dofIPTMap[NuTo::Node::eDof::COORDINATES]            = NuTo::Interpolation::EQUIDISTANT1;
        dofIPTMap[NuTo::Node::eDof::RELATIVEHUMIDITY]       = NuTo::Interpolation::EQUIDISTANT1;
        dofIPTMap[NuTo::Node::eDof::WATERVOLUMEFRACTION]    = NuTo::Interpolation::EQUIDISTANT1;

        // Constitutive output tests

        ConstitutiveOutputTests<1>(dofIPTMap);
        ConstitutiveOutputTests<2>(dofIPTMap);
        ConstitutiveOutputTests<3>(dofIPTMap);

        // Simulation Tests --- Compare results to paper

        SimulationTest<1>({16},
                          {0.16},
                          {true},
                          dofIPTMap);

        SimulationTest<2>({16,1},
                          {0.16,0.1},
                          {true,false},
                          dofIPTMap);

        SimulationTest<2>({1,16},
                          {0.1,0.16},
                          {false,true},
                          dofIPTMap);


        SimulationTest<3>({16,1,1},
                          {0.16,0.1,0.1},
                          {true,false,false},
                          dofIPTMap);

        SimulationTest<3>({1,16,1},
                          {0.1,0.16,0.1},
                          {false,true,false},
                          dofIPTMap);

        SimulationTest<3>({1,1,16},
                          {0.1,0.1,0.16},
                          {false,false,true},
                          dofIPTMap);

        // Simulation Tests --- Interpolation type


        dofIPTMap[NuTo::Node::eDof::RELATIVEHUMIDITY]       = NuTo::Interpolation::EQUIDISTANT2;
        dofIPTMap[NuTo::Node::eDof::WATERVOLUMEFRACTION]    = NuTo::Interpolation::EQUIDISTANT1;
        RunSimulationTest_AllDiemnsions(dofIPTMap,2);


        dofIPTMap[NuTo::Node::eDof::RELATIVEHUMIDITY]       = NuTo::Interpolation::EQUIDISTANT1;
        dofIPTMap[NuTo::Node::eDof::WATERVOLUMEFRACTION]    = NuTo::Interpolation::EQUIDISTANT2;
        RunSimulationTest_AllDiemnsions(dofIPTMap,2);


        dofIPTMap[NuTo::Node::eDof::RELATIVEHUMIDITY]       = NuTo::Interpolation::EQUIDISTANT2;
        dofIPTMap[NuTo::Node::eDof::WATERVOLUMEFRACTION]    = NuTo::Interpolation::EQUIDISTANT2;
        RunSimulationTest_AllDiemnsions(dofIPTMap,2);

    }
    catch (NuTo::Exception& e)
    {
        std::cout << e.ErrorMessage() << std::endl;
        std::cout << "\n\n\n Errors occured! \n\n\n" << std::endl;
        return 1;
    }

    std::cout << "Test finished" << std::endl;
    return 0;
}
