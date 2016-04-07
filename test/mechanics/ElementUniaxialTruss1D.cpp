/*
 * ElementUniaxialTruss1D.cpp
 *
 *  Created on: 13 May 2015
 *      Author: ttitsche
 */

#include "../test/mechanics/ElementUniaxialTest.h"

std::string directory = "";

void Run(NuTo::Interpolation::eTypeOrder rTypeOrder)
{
    NuToTest::ElementUniaxialTest test;

#ifdef ENABLE_VISUALIZE
    test.visualizationDirectory = directory;
#endif

    NuTo::Structure myStructure(1);
    myStructure.SetShowTime(false);

    int numElementsX = 3;

    //create nodes
    int numNodesX = numElementsX+1;
    double deltaX = test.lX/(numElementsX);

    int nodeNum = 0;
    for (int countX=0; countX<numNodesX; countX++)
    {
        NuTo::FullVector<double,Eigen::Dynamic> coordinates(1);
        coordinates(0) = countX*deltaX;
        myStructure.NodeCreate(nodeNum,coordinates);
        nodeNum++;
    }

    int myInterpolationType = myStructure.InterpolationTypeCreate("Truss1D");
    myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::DISPLACEMENTS, rTypeOrder);

    //create elements
    for (int countX=0; countX<numElementsX; countX++)
    {
        NuTo::FullVector<int,Eigen::Dynamic> nodes(2);
        nodes(0) = countX;
        nodes(1) = countX+1;
        myStructure.ElementCreate(myInterpolationType, nodes);
    }

    myStructure.SetVerboseLevel(10);
    myStructure.ElementTotalConvertToInterpolationType();

    int mySection = myStructure.SectionCreate("Truss");
    myStructure.SectionSetArea(mySection, test.lZ*test.lY);
    myStructure.ElementTotalSetSection(mySection);
    myStructure.SetNumProcessors(2);

    test.Run(myStructure);
}

int main(int argc, char* argv[])
{

    boost::filesystem::path path = boost::filesystem::system_complete(boost::filesystem::path( argv[0] ));
    directory = path.parent_path().string();

    try
    {
        Run(NuTo::Interpolation::EQUIDISTANT1);
//        Run(NuTo::Interpolation::EQUIDISTANT2);
//        Run(NuTo::Interpolation::EQUIDISTANT3);
//        Run(NuTo::Interpolation::EQUIDISTANT4);
//        Run(NuTo::Interpolation::LOBATTO2);
//        Run(NuTo::Interpolation::LOBATTO3);
//        Run(NuTo::Interpolation::LOBATTO4);
    }
    catch (NuTo::Exception& e)
    {
        std::cout << "## Test failed ##" << std::endl;
        std::cout << e.ErrorMessage();
        return -1;
    }

    return 0;
}

