#include <iostream>

#include "TestClasses/Numbering.h"

#include "TestFiles/MyTruss1D2N.h"
#include "TestFiles/MyHeatedPlate.h"

#include "mechanics/nodes/NodeEnum.h"



void testNumbering()
{
    Numbering numb;

    numb.addGlobalDof(NuTo::Node::eDof::COORDINATES, 1, 1);
    numb.addGlobalDof(NuTo::Node::eDof::COORDINATES, 1, 2);
    numb.addGlobalDof(NuTo::Node::eDof::TEMPERATURE, 0, 1);
    numb.setGlobalDofValue(NuTo::Node::eDof::COORDINATES, 1, 1, 1);

    numb.printAllDofTypes();
}


void testMyTruss1D2N()
{
    MyTruss1D2N truss;

    truss.Run();
}

void testMyHeatedPlate()
{
    MyHeatedPlate plate;
    plate.Run();
}


int main()
{
    testNumbering();

    //testMyTruss1D2N();

    //testMyHeatedPlate();

    return 0;
}
