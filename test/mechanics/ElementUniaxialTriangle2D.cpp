/*
 * ElementUniaxialTriangle2D.cpp
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

    NuTo::Structure myStructure(2);
    myStructure.SetShowTime(false);

    int numElementsX = 3;
    int numElementsY = 5;

    //create nodes
    int numNodesX = numElementsX+1;
    int numNodesY = numElementsY+1;
    double deltaX = test.lX/(numElementsX);
    double deltaY = test.lY/(numElementsY);

    int nodeNum = 0;
    for (int countY=0; countY<numNodesY; countY++)
    {
        for (int countX=0; countX<numNodesX; countX++)
        {
            NuTo::FullVector<double,Eigen::Dynamic> coordinates(2);
            coordinates(0) = countX*deltaX;
            coordinates(1) = countY*deltaY;
            myStructure.NodeCreate(nodeNum, coordinates);
            nodeNum++;
        }
    }

    int myInterpolationType = myStructure.InterpolationTypeCreate("Triangle2D");
    myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::DISPLACEMENTS, rTypeOrder);

    //create elements
    for (int countY=0; countY<numElementsY; countY++)
    {
        for (int countX=0; countX<numElementsX; countX++)
        {
            NuTo::FullVector<int,Eigen::Dynamic> nodes(3);
            nodes(0) = countX  +  countY   *numNodesX;
            nodes(1) = countX+1+  countY   *numNodesX;
            nodes(2) = countX+1+ (countY+1)*numNodesX;
            myStructure.ElementCreate(myInterpolationType, nodes);

            nodes(0) = countX  +  countY   *numNodesX;
            nodes(2) = countX+1+ (countY+1)*numNodesX;
            nodes(1) = countX  + (countY+1)*numNodesX;
            myStructure.ElementCreate(myInterpolationType, nodes);

        }
    }

    myStructure.SetVerboseLevel(10);
    myStructure.ElementTotalConvertToInterpolationType();

    int mySection = myStructure.SectionCreate("Plane_Stress");
    myStructure.SectionSetThickness(mySection, test.lZ);
    myStructure.ElementTotalSetSection(mySection);


    test.Run(myStructure);
}

int main(int argc, char* argv[])
{

    boost::filesystem::path path = boost::filesystem::system_complete(boost::filesystem::path( argv[0] ));
    directory = path.parent_path().string();

    try
    {
        Run(NuTo::Interpolation::EQUIDISTANT1);
        Run(NuTo::Interpolation::EQUIDISTANT2);
        Run(NuTo::Interpolation::EQUIDISTANT3);
        Run(NuTo::Interpolation::EQUIDISTANT4);
    }
    catch (NuTo::Exception& e)
    {
        std::cout << "## Test failed ##" << std::endl;
        std::cout << e.ErrorMessage();
        return -1;
    }

    return 0;
}

