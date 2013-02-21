// $Id: StructureElement.cpp 539 2011-05-16 14:00:36Z unger3 $

#include <assert.h>
#include "nuto/mechanics/structures/unstructured/Structure.h"


//! @brief creates a lattice mesh from the positions of the circles
//! @parameters rTypeOfSpecimen 0 box, 1 dogbone
//! @parameters rBoundingBox box for the spheres (3*2 matrix)
//! @parameters rCircles (coordinates x,y and radius)
//! @parameters rTriangles (triangles connecting the circle centers)
void NuTo::Structure::MeshCreateLattice2D(int rTypeOfSpecimen, FullMatrix<double>& rBoundingBox, NuTo::FullMatrix<double>& rCircles, NuTo::FullMatrix<double>& rTriangles)
{
/*    if (mDimension!=2)
    	throw MechanicsException("[NuTo::Structure::MeshCreateLattice2D] structure is not 2D.");

    //create nodes
    std::string dofs("COORDINATES DISPLACEMENTS VELOCITIES ACCELERATIONS ROTATIONS ANGULAR_VELOCITIES ANGULAR_ACCELERATIONS RADIUS");
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
	//only relevant for dog bone specimens
	double D(rBoundingBox(0,1)-rBoundingBox(0,0));
	double X1(rBoundingBox(0,0)-0.525*D);
	double X2(rBoundingBox(0,1)+0.525*D);
	double Y(0.5*(rBoundingBox(1,0) + rBoundingBox(1,1)));
	double R2(0.725*D*0.725*D);
	std::vector<NodeBase*> nodeVector(3);
	for (int count=0; count<rTriangles.GetNumRows(); count++)
	{
		nodeVector[0] = NodeGetNodePtr(rTriangles(count,0));
		nodeVector[1] = NodeGetNodePtr(rTriangles(count,1));
		nodeVector[2] = NodeGetNodePtr(rTriangles(count,2));
		if (rTypeOfSpecimen==1)
		{
			//for the dogbone specimen, check if the center of gravity is inside the circles that are cut out,
			//because qhull is meshing this as well
			double coordX = (nodeVector[0]->GetCoordinate(0)+nodeVector[1]->GetCoordinate(0)+nodeVector[2]->GetCoordinate(0))/3.;
			double coordY = (nodeVector[0]->GetCoordinate(1)+nodeVector[1]->GetCoordinate(1)+nodeVector[2]->GetCoordinate(1))/3.;
			if ((coordX-X1)*(coordX-X1)+(coordY-Y)*(coordY-Y)<R2)
				continue;
			if ((coordX-X2)*(coordX-X2)+(coordY-Y)*(coordY-Y)<R2)
				continue;
		}
//		ElementCreate(NuTo::Element::LATTICE2D, nodeVector, NuTo::ElementData::CONSTITUTIVELAWIP, NuTo::IpData::STATICDATAWEIGHTCOORDINATES2D);
	}
*/
}

//! @brief creates a lattice mesh from the positions of the spheres and the bounding box
//! @parameters rBoundingBox (min and max for x and y)
//! @parameters rSpheres (coordinates x,y,z and radius)
void NuTo::Structure::MeshCreateLattice3D(int rTypeOfSpecimen, FullMatrix<double>& rBoundingBox, NuTo::FullMatrix<double>& rSpheres, NuTo::FullMatrix<double>& rTetraeders)
{
/*
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
	*/
}
