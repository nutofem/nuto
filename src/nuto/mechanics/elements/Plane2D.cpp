// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include <algorithm>

#include "nuto/mechanics/constitutive/mechanics/DeformationGradient1D.h"
#include "nuto/mechanics/elements/Plane2D.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include <assert.h>


#ifdef HAVE_CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2_algorithms.h>
#endif //HAVE_CGAL


NuTo::Plane2D::Plane2D(const NuTo::StructureBase* rStructure, ElementData::eElementDataType rElementDataType,
		IntegrationType::eIntegrationType rIntegrationType, IpData::eIpDataType rIpDataType) :
        Plane(rStructure, rElementDataType, rIntegrationType, rIpDataType)
{}

//! @brief calculates the local coordinates of the nodes
//! @param localCoordinates vector with already correct size allocated
//! this can be checked with an assertation
void NuTo::Plane2D::CalculateLocalCoordinates(std::vector<double>& rLocalCoordinates)const
{
	assert((int)rLocalCoordinates.size()==2*GetNumNodes());
    for (int theNode=0; theNode<GetNumNodes(); theNode++)
    {
        GetNode(theNode)->GetCoordinates2D(&(rLocalCoordinates[2*theNode]));
    }
}

//! @brief calculates the local displacements of the nodes
//! @param localDisplacements vector with already correct size allocated
//! this can be checked with an assertation
void NuTo::Plane2D::CalculateLocalDisplacements(std::vector<double>& rLocalDisplacements)const
{
	assert((int)rLocalDisplacements.size()==2*GetNumNodes());
    for (int theNode=0; theNode<GetNumNodes(); theNode++)
    {
        GetNode(theNode)->GetDisplacements2D(&(rLocalDisplacements[2*theNode]));
    }
}

// build global row dofs
void NuTo::Plane2D::CalculateGlobalRowDofs(std::vector<int>& rGlobalRowDofs, int numDispDofs, int numTempDofs) const
{
    rGlobalRowDofs.resize(numDispDofs+numTempDofs);
    for (int nodeCount = 0; nodeCount < this->GetNumNodes(); nodeCount++)
    {
        const NodeBase * nodePtr(GetNode(nodeCount));
        if (nodePtr->GetNumDisplacements()>0 && numDispDofs>0)
        {
            rGlobalRowDofs[2 * nodeCount    ] = nodePtr->GetDofDisplacement(0);
            rGlobalRowDofs[2 * nodeCount + 1] = nodePtr->GetDofDisplacement(1);
        }
        if (nodePtr->GetNumTemperatures()>0 && numTempDofs>0)
        {
            rGlobalRowDofs[numDispDofs + nodeCount ] = nodePtr->GetDofTemperature();
        }
    }
}

// build global column dof
void NuTo::Plane2D::CalculateGlobalColumnDofs(std::vector<int>& rGlobalColumnDofs, int numDispDofs, int numTempDofs) const
{
    int NumNonlocalElements(GetNumNonlocalElements());
	if (GetNumNonlocalElements()==0)
	    this->CalculateGlobalRowDofs(rGlobalColumnDofs,numDispDofs,numTempDofs);
    else
    {
        rGlobalColumnDofs.resize(numDispDofs+numTempDofs);
	    const std::vector<const ElementBase*>& nonlocalElements(GetNonlocalElements());
        int shift(0);
        for (int theNonlocalElement=0; theNonlocalElement<NumNonlocalElements; theNonlocalElement++)
        {
            const ElementBase* nonlocalElement(nonlocalElements[theNonlocalElement]);
        	for (int nodeCount = 0; nodeCount < nonlocalElement->GetNumNodes(); nodeCount++)
            {
                const NodeBase *nodePtr = nonlocalElement->GetNode(nodeCount);
                if (nodePtr->GetNumDisplacements()>0 && numDispDofs>0)
                {
					rGlobalColumnDofs[shift + 2 * nodeCount    ] = nodePtr->GetDofDisplacement(0);
					rGlobalColumnDofs[shift + 2 * nodeCount + 1] = nodePtr->GetDofDisplacement(1);
                }
            }
            shift+=2*nonlocalElement->GetNumNodes();
        }
        assert(shift==numDispDofs);
        for (int nodeCount = 0; nodeCount < this->GetNumNodes(); nodeCount++)
        {
            const NodeBase * nodePtr(GetNode(nodeCount));
            if (nodePtr->GetNumTemperatures()>0 && numTempDofs>0)
            {
            	rGlobalColumnDofs[numDispDofs + nodeCount ] = nodePtr->GetDofTemperature();
            }
        }
    }
}

//! @brief calculates the area of a plane element via the nodes (probably faster than sum over integration points)
//! @return Area
double NuTo::Plane2D::CalculateArea()const
{
    double coordinates1[2];
    double coordinates2[2];
    GetNode(GetNumNodes()-1)->GetCoordinates2D(coordinates2);
	double area(0);
	for (int theNode = 0; theNode<GetNumNodes(); theNode++)
	{
		coordinates1[0] = coordinates2[0];
		coordinates1[1] = coordinates2[1];
		GetNode(theNode)->GetCoordinates2D(coordinates2);
		area += coordinates1[0]*coordinates2[1] - coordinates1[1]*coordinates2[0];
	}
    return 0.5*area;
}


//! @brief checks if a node is inside a polygon
//! @param rPoint (input) ... a pointer to a 2D tuple containing the coordinates
//! @param rPoint (input) ... a pointer to a vector of 2D tuples representing the polyline
//! @return True if coordinates are within the element, False otherwise
//! @todo move to geometry class
bool NuTo::Plane2D::CheckPointInsidePolygon( const std::tuple<double,double> *rPoint, const std::vector<std::tuple<double,double> > * rPolygon)const
{
#ifdef HAVE_CGAL
	/// this is from CGAL examples
	
	typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
	typedef K::Point_2 Point;
	
	Point pt(std::get<0>(*rPoint),std::get<1>(*rPoint));
	
	std::vector<Point> polygon;
	
	for (size_t i=0; i<rPolygon->size();++i){
		//std::cout << "The Polygon point " << std::get<0>((*rPolygon)[i]) << " " << std::get<1>((*rPolygon)[i]) << std::endl;
		polygon.push_back(Point( std::get<0>((*rPolygon)[i]), std::get<1>((*rPolygon)[i])));
	}
	polygon.push_back(Point( std::get<0>((*rPolygon)[0]), std::get<1>((*rPolygon)[0])));
	//~ Point points[] = { Point(0,0), Point(5.1,0), Point(1,1), Point(0.5,6)};

	bool inside=false;
	//std::cout << "The point " << pt;
	switch(CGAL::bounded_side_2(polygon.begin(), polygon.end(),pt, K())) {
		case CGAL::ON_BOUNDED_SIDE :
			//std::cout << " is inside the polygon.\n";
			inside=true;
			break;
		case CGAL::ON_BOUNDARY:
			//std::cout << " is on the polygon boundary.\n";
			inside=true;
			break;
		case CGAL::ON_UNBOUNDED_SIDE:
			//std::cout << " is outside the polygon.\n";
			inside=false;
			break;
	}
	return inside;
#else //HAVE_CGAL
/*! This function is a modified function of
 *
 *  public domain function by Darel Rex Finley
 *  http://alienryderflex.com/polygon/
 *
 *  double fuzzy  (internal constant)  =  fuzzy control
 *  bool   inside (internal variable)  =  check result
 *  double x, y   (internal constants) =  point to be tested
 *
 *  The function will return TRUE if the point rPoint is inside the polygon, or
 *  NO if it is not.  If the point is exactly on the edge of the polygon,
 *  then the function will return TRUE.
 *
 *  Note that division by zero is avoided because the division is protected
 *  by the "if" clause which surrounds it.
 */
	double fuzzy=1e-14;
	bool  inside=false ;

	const double x=std::get<0>(*rPoint);
	const double y=std::get<1>(*rPoint);

	unsigned int j=rPolygon->size()-1 ;
	for (unsigned int i=0; i<rPolygon->size(); i++)
	{
		if ( (std::get<1>((*rPolygon)[i])-y<fuzzy && std::get<1>((*rPolygon)[j])-y>fuzzy)
		||   (std::get<1>((*rPolygon)[j])-y<fuzzy && std::get<1>((*rPolygon)[i])-y>fuzzy) )
		{
			if (std::get<0>((*rPolygon)[i])
				+(y-std::get<1>((*rPolygon)[i]))
					/(std::get<1>((*rPolygon)[j])-std::get<1>((*rPolygon)[i]))
					*(std::get<0>((*rPolygon)[j])
				-std::get<0>((*rPolygon)[i]))-x<fuzzy)
			{
				inside=!inside;
			}
		}
		j=i;
	}

	return inside;
#endif //HAVE_CGAL
}

//! @brief help function used in the Cracklength calculation to sort
//! @parameter array1 (coordinates of point1 (x,y) and then the angle wrt mean
bool CompareFunctionCrackAngle (const boost::array<double,3>&  array1, const boost::array<double,3>&  array2)
{
	return (array1[2]<array2[2]);
}

//! @brief calculate the length of a crack passing through the center of gravity and intersecting two edges
//! @parameter alpha... angle of the crack
double NuTo::Plane2D::CalculateCrackLength2D(double rAlpha)const
{
	//get all nodes and calculate mean
	std::vector<boost::array<double,3> > coordinates;
	coordinates.resize(GetNumNodes());
	boost::array<double,2> mean = {{0.,0.}};
	for (int theNode = 0; theNode<GetNumNodes(); theNode++)
	{
		GetNode(theNode)->GetCoordinates2D(&(coordinates[theNode][0]));
		mean[0]+=coordinates[theNode][0];
		mean[1]+=coordinates[theNode][1];
	}
	mean[0]/=(double)GetNumNodes();
	mean[1]/=(double)GetNumNodes();

	//calculate angle wrt arithmetic center
	for (int theNode = 0; theNode<GetNumNodes(); theNode++)
	{
		coordinates[theNode][2]=atan2(coordinates[theNode][1]-mean[1],coordinates[theNode][0]-mean[0]);
	}

	//print unsorted
//	std::cout << "before sort" << "\n";
//	for (int theNode = 0; theNode<GetNumNodes(); theNode++)
//	{
//		std::cout << " coordinates " << coordinates[theNode][0] << " " << coordinates[theNode][1] << ", angle " << coordinates[theNode][2] << "\n";
//	}

	//sort wrt to angle
    std::sort (coordinates.begin(), coordinates.end(), CompareFunctionCrackAngle);

//	std::cout << "after sort" << "\n";
//	for (int theNode = 0; theNode<GetNumNodes(); theNode++)
//	{
//		std::cout << " coordinates " << coordinates[theNode][0] << " " << coordinates[theNode][1] << ", angle " << coordinates[theNode][2] << "\n";
//	}

	//calculate center of gravity
	coordinates.resize(GetNumNodes()+1);
	coordinates[GetNumNodes()] = coordinates[0];
	boost::array<double,2> centerOfGravity = {{0.,0.}};
	double area(0);
	for (int theNode = 0; theNode<GetNumNodes(); theNode++)
	{
		double tmp(coordinates[theNode][0]*coordinates[theNode+1][1]-coordinates[theNode+1][0]*coordinates[theNode][1]);
		area += tmp;
		centerOfGravity[0] += (coordinates[theNode][0]+coordinates[theNode+1][0])*tmp;
		centerOfGravity[1] += (coordinates[theNode][1]+coordinates[theNode+1][1])*tmp;
	}
    area*=0.5;
    centerOfGravity[0]/=6.*area;
    centerOfGravity[1]/=6.*area;

    //std::cout << "area "<< area << "\n";
    //std::cout << "center of gravity "<< centerOfGravity[0] << " " << centerOfGravity[1] << "\n";
    //std::cout << "global angle "<< rAlpha*180./M_PI << "\n";

    //calculate intersection with all edges
    double u[2];
    double v[2];
    double w[2];
    v[0] = cos(rAlpha);
    v[1] = sin(rAlpha);
    double s;
    double denom;
    double intersectionPoints[2][2];
    int numIntersectionPoints(0);
	for (int theNode = 0; theNode<GetNumNodes(); theNode++)
	{
        u[0] = coordinates[theNode+1][0]-coordinates[theNode][0];
        u[1] = coordinates[theNode+1][1]-coordinates[theNode][1];
        w[0] = coordinates[theNode][0] - mean[0];
        w[1] = coordinates[theNode][1] - mean[1];
        denom = v[0]*u[1]-v[1]*u[0];
        if (fabs(denom)>1e-10)
        {
            s = (v[1]*w[0]-v[0]*w[1])/denom;
            if (s>=0 && s<1)
            {
            	if (numIntersectionPoints>1)
            	{
            		throw MechanicsException("[NuTo::Plane2D::CalculateCrackLength2D] Found more than two intersection points, there is something wrong.");
            	}
      //      	std::cout << "s " << s << " theNode " << theNode << "\n";
      //      	std::cout << "coordinates 1 " << coordinates[theNode][0] << " " << coordinates[theNode][1] << "\n";
      //      	std::cout << "coordinates 2 " << coordinates[theNode+1][0] << " " << coordinates[theNode+1][1] << "\n";
      //      	std::cout << "u " << u[0] << " " << u[1] << "\n";
            	intersectionPoints[numIntersectionPoints][0] = coordinates[theNode][0] + s*u[0];
            	intersectionPoints[numIntersectionPoints][1] = coordinates[theNode][1] + s*u[1];
      //      	std::cout << "intersection " << intersectionPoints[numIntersectionPoints][0] << " " << intersectionPoints[numIntersectionPoints][1] << "\n";

            	numIntersectionPoints++;
            }
        }
	}
	if (numIntersectionPoints!=2)
	{
		throw MechanicsException("[NuTo::Plane2D::CalculateCrackLength2D] Did not found two intersection points, there is something wrong.");
	}
    //std::cout << "intersection point 1 "<< intersectionPoints[0][0] << " " << intersectionPoints[0][1] << "\n";
    //std::cout << "intersection point 2 "<< intersectionPoints[1][0] << " " << intersectionPoints[1][1] << "\n";

	//std::cout << "crack length " << sqrt((intersectionPoints[0][0]-intersectionPoints[1][0])*(intersectionPoints[0][0]-intersectionPoints[1][0])+
	//	    (intersectionPoints[0][1]-intersectionPoints[1][1])*(intersectionPoints[0][1]-intersectionPoints[1][1]))<< "\n";
	//calculate length of crack and return
	return sqrt((intersectionPoints[0][0]-intersectionPoints[1][0])*(intersectionPoints[0][0]-intersectionPoints[1][0])+
			    (intersectionPoints[0][1]-intersectionPoints[1][1])*(intersectionPoints[0][1]-intersectionPoints[1][1]));
}

//! @brief cast the base pointer to an ElementPlane, otherwise throws an exception
const NuTo::Plane2D* NuTo::Plane2D::AsPlane2D()const
{
    return this;
}

//! @brief cast the base pointer to an Plane, otherwise throws an exception
NuTo::Plane2D* NuTo::Plane2D::AsPlane2D()
{
    return this;
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::Plane2D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::Plane2D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::Plane2D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::Plane2D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::Plane2D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::Plane2D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::Plane2D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize Plane2D" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Plane);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize Plane2D" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::Plane2D)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::Plane2D)
#endif // ENABLE_SERIALIZATION

