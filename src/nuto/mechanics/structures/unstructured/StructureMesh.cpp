// $Id: StructureElement.cpp 539 2011-05-16 14:00:36Z unger3 $

#include <assert.h>
#include "nuto/mechanics/structures/unstructured/Structure.h"


//! @brief creates a lattice mesh from the positions of the circles and the bounding box
//! @parameters rBoundingBox (min and max for x and y)
//! @parameters rCircles (coordinates x,y and radius)
void NuTo::Structure::MeshCreateLattice2D(FullMatrix<double>& rBoundingBox, NuTo::FullMatrix<double>& rCircles, NuTo::FullMatrix<double>& rTriangles)
{
    if (mDimension!=2)
    	throw MechanicsException("[NuTo::Structure::MeshCreateLattice2D] structure is not 2D.");

    //create nodes
    std::string dofs("COORDINATES DISPLACEMENTS ROTATIONS RADIUS");
	for (int count=0; count<rCircles.GetNumRows(); count++)
	{
		NuTo::FullMatrix<double> coordinates(2,1);
		coordinates(0,0) = rCircles(count,0);
		coordinates(1,0) = rCircles(count,1);
		int theNode = NodeCreate(dofs, coordinates);
		NodeBase* nodePtr = NodeGetNodePtr(theNode);
		nodePtr->SetRadius(&(rCircles(count,2)));
	}

	//create elements
	std::vector<NodeBase*> nodeVector(3);
	for (int count=0; count<rTriangles.GetNumRows(); count++)
	{
		nodeVector[0] = NodeGetNodePtr(rTriangles(count,0));
		nodeVector[1] = NodeGetNodePtr(rTriangles(count,1));
		nodeVector[2] = NodeGetNodePtr(rTriangles(count,2));
		ElementCreate(NuTo::Element::LATTICE2D, nodeVector, NuTo::ElementData::CONSTITUTIVELAWIP, NuTo::IpData::STATICDATAWEIGHTCOORDINATES2D);
	}
}

//! @brief creates a lattice mesh from the positions of the spheres and the bounding box
//! @parameters rBoundingBox (min and max for x and y)
//! @parameters rSpheres (coordinates x,y,z and radius)
void NuTo::Structure::MeshCreateLattice3D(FullMatrix<double>& rBoundingBox, NuTo::FullMatrix<double>& rSpheres, NuTo::FullMatrix<double>& rTetraeders)
{
    if (mDimension!=3)
    	throw MechanicsException("[NuTo::Structure::MeshCreateLattice3D] structure is not 3D.");
	std::string dofs("COORDINATES DISPLACEMENTS RADIUS");
	for (int count=0; count<rSpheres.GetNumRows(); count++)
	{
		NuTo::FullMatrix<double> coordinates(3,1);
		coordinates(0,0) = rSpheres(count,0);
		coordinates(1,0) = rSpheres(count,1);
		coordinates(2,0) = rSpheres(count,2);
		int theNode = NodeCreate(dofs, coordinates);
		NodeBase* nodePtr = NodeGetNodePtr(theNode);
		nodePtr->SetRadius(&(rSpheres(count,3)));
	}
}
