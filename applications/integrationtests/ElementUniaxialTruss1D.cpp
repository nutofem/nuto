/*
 * ElementUniaxialTruss1D.cpp
 *
 *  Created on: 13 May 2015
 *      Author: ttitsche
 */

#include "ElementUniaxialTest.h"

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
        Eigen::VectorXd coordinates(1);
        coordinates(0) = countX*deltaX;
        myStructure.NodeCreate(nodeNum,coordinates);
        nodeNum++;
    }

    int myInterpolationType = myStructure.InterpolationTypeCreate("Truss1D");
    myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::eDof::DISPLACEMENTS, rTypeOrder);

    //create elements
    std::vector<int> nodes(2);
    for (int countX=0; countX<numElementsX; countX++)
    {
        nodes[0] = countX;
        nodes[1] = countX+1;
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
        Run(NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
        Run(NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
        Run(NuTo::Interpolation::eTypeOrder::EQUIDISTANT3);
        Run(NuTo::Interpolation::eTypeOrder::EQUIDISTANT4);
        Run(NuTo::Interpolation::eTypeOrder::LOBATTO2);
        Run(NuTo::Interpolation::eTypeOrder::LOBATTO3);
        Run(NuTo::Interpolation::eTypeOrder::LOBATTO4);
    }
    catch (NuTo::Exception& e)
    {
        std::cout << "## Test failed ##" << std::endl;
        std::cout << e.ErrorMessage();
        return -1;
    }

    return 0;
}

