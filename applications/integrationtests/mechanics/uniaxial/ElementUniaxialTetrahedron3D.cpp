/*
 * ElementUniaxialTetrahedron3D.cpp
 *
 *  Created on: 19 May 2015
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
                myStructure.NodeCreate(nodeNum, Eigen::Vector3d({iX*deltaX, iY*deltaY, iZ*deltaZ}));
                nodeNum++;
            }

    int myInterpolationType = myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TETRAHEDRON3D);
    myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::eDof::DISPLACEMENTS, rTypeOrder);

    for (int iZ=0; iZ<numElementsZ; iZ++)
        for (int iY=0; iY<numElementsY; iY++)
            for (int iX=0; iX<numElementsX; iX++)
            {
                std::vector<int> nodes(8);
                nodes[0] = iX   +  iY    * numNodesX +  iZ    * numNodesX * numNodesY;
                nodes[1] = iX+1 +  iY    * numNodesX +  iZ    * numNodesX * numNodesY;
                nodes[2] = iX+1 + (iY+1) * numNodesX +  iZ    * numNodesX * numNodesY;
                nodes[3] = iX   + (iY+1) * numNodesX +  iZ    * numNodesX * numNodesY;
                nodes[4] = iX   +  iY    * numNodesX + (iZ+1) * numNodesX * numNodesY;
                nodes[5] = iX+1 +  iY    * numNodesX + (iZ+1) * numNodesX * numNodesY;
                nodes[6] = iX+1 + (iY+1) * numNodesX + (iZ+1) * numNodesX * numNodesY;
                nodes[7] = iX   + (iY+1) * numNodesX + (iZ+1) * numNodesX * numNodesY;

                std::vector<int> nodesTet0({nodes[0], nodes[1], nodes[3], nodes[7]});
                std::vector<int> nodesTet1({nodes[0], nodes[1], nodes[7], nodes[4]});
                std::vector<int> nodesTet2({nodes[5], nodes[4], nodes[7], nodes[1]});
                std::vector<int> nodesTet3({nodes[6], nodes[5], nodes[7], nodes[1]});
                std::vector<int> nodesTet4({nodes[2], nodes[7], nodes[1], nodes[6]});
                std::vector<int> nodesTet5({nodes[2], nodes[3], nodes[1], nodes[7]});

                myStructure.ElementCreate(myInterpolationType, nodesTet0);
                myStructure.ElementCreate(myInterpolationType, nodesTet1);
                myStructure.ElementCreate(myInterpolationType, nodesTet2);
                myStructure.ElementCreate(myInterpolationType, nodesTet3);
                myStructure.ElementCreate(myInterpolationType, nodesTet4);
                myStructure.ElementCreate(myInterpolationType, nodesTet5);

            }

    int allElements = myStructure.GroupCreate(NuTo::eGroupId::Elements);
    myStructure.GroupAddElementFromType(allElements, myInterpolationType);
    double volume = myStructure.ElementGroupGetVolume(allElements);
    std::cout << "######### VOLUME: " << volume << std::endl;

    myStructure.SetVerboseLevel(10);
    myStructure.ElementTotalConvertToInterpolationType(1.e-6, 3);

    int mySection = myStructure.SectionCreate(NuTo::eSectionType::VOLUME);
    myStructure.ElementTotalSetSection(mySection);


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
    }
    catch (NuTo::Exception& e)
    {
        std::cout << "## Test failed ##" << std::endl;
        std::cout << e.ErrorMessage();
        return -1;
    }

    return 0;
}

