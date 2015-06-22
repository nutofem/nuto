/*
 * ElementUniaxialBrick3D.cpp
 *
 *  Created on: 20 May 2015
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

    NuTo::Structure myStructure(3);
    myStructure.SetShowTime(false);


    int numElementsX = 10;
    int numElementsY = 2;
    int numElementsZ = 2;


    //create nodes
    int numNodesX = numElementsX+1;
    int numNodesY = numElementsY+1;
    int numNodesZ = numElementsZ+1;

    double deltaX = test.lX/(numElementsX);
    double deltaY = test.lY/(numElementsY);
    double deltaZ = test.lZ/(numElementsZ);

    int nodeNum = 0;
    for (int iZ=0; iZ<numNodesZ; iZ++)
        for (int iY=0; iY<numNodesY; iY++)
            for (int iX=0; iX<numNodesX; iX++)
            {
                NuTo::FullVector<double,Eigen::Dynamic> coordinates(3);
                coordinates(0) = iX*deltaX;
                coordinates(1) = iY*deltaY;
                coordinates(2) = iZ*deltaZ;
                myStructure.NodeCreate(nodeNum, coordinates);
                nodeNum++;
            }

    int myInterpolationType = myStructure.InterpolationTypeCreate("Brick3D");
    myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::DISPLACEMENTS, rTypeOrder);

    for (int iZ=0; iZ<numElementsZ; iZ++)
        for (int iY=0; iY<numElementsY; iY++)
            for (int iX=0; iX<numElementsX; iX++)
            {
                NuTo::FullVector<int,Eigen::Dynamic> nodes(8);
                nodes(0) = iX   +  iY    * numNodesX +  iZ    * numNodesX * numNodesY;
                nodes(1) = iX+1 +  iY    * numNodesX +  iZ    * numNodesX * numNodesY;
                nodes(2) = iX+1 + (iY+1) * numNodesX +  iZ    * numNodesX * numNodesY;
                nodes(3) = iX   + (iY+1) * numNodesX +  iZ    * numNodesX * numNodesY;
                nodes(4) = iX   +  iY    * numNodesX + (iZ+1) * numNodesX * numNodesY;
                nodes(5) = iX+1 +  iY    * numNodesX + (iZ+1) * numNodesX * numNodesY;
                nodes(6) = iX+1 + (iY+1) * numNodesX + (iZ+1) * numNodesX * numNodesY;
                nodes(7) = iX   + (iY+1) * numNodesX + (iZ+1) * numNodesX * numNodesY;

                myStructure.ElementCreate(myInterpolationType, nodes);

            }

    int allElements = myStructure.GroupCreate("Elements");
    myStructure.GroupAddElementFromType(allElements, myInterpolationType);
    double volume = myStructure.ElementGroupGetVolume(allElements);
    std::cout << "######### VOLUME: " << volume << std::endl;

    myStructure.SetVerboseLevel(10);
    myStructure.ElementTotalConvertToInterpolationType(1.e-6, 3);

    int mySection = myStructure.SectionCreate("VOLUME");
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
    }
    catch (NuTo::Exception& e)
    {
        std::cout << "## Test failed ##" << std::endl;
        std::cout << e.ErrorMessage();
        return -1;
    }

    return 0;
}

