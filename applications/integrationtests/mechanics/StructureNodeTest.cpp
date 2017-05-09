#include <stdlib.h>
#include <iostream>
#include "BoostUnitTest.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/MechanicsEnums.h"

void CheckElementsBelongingToNode(int node, std::vector<int> expected)
{

    NuTo::Structure s(1);
    s.SetVerboseLevel(0);

    int numElements = 3;
    int numNodes = numElements+1;

    for (int i=0; i<numNodes ; i++) {
        Eigen::VectorXd coordinates(1);
        coordinates(0) =  i;
        s.NodeCreate(i,coordinates);
    }

    int myInterpolationType = s.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TRUSS1D);
    s.InterpolationTypeAdd(myInterpolationType, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);

    for (int i=0; i< numElements; i++) {
        s.ElementCreate(myInterpolationType, {i,i+1});
    }

    std::vector<int> elementsBelongingToNode;
    s.NodeGetElements(node,elementsBelongingToNode);

    BOOST_CHECK_EQUAL_COLLECTIONS(elementsBelongingToNode.begin(), elementsBelongingToNode.end(), expected.begin(), expected.end());
}

BOOST_AUTO_TEST_CASE(CheckNodeGetElement)
{
    CheckElementsBelongingToNode(0, {0});
    CheckElementsBelongingToNode(1, {0, 1});
    CheckElementsBelongingToNode(2, {1 ,2});
}
