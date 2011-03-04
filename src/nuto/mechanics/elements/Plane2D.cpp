// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/mechanics/DeformationGradient1D.h"
#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal1x1.h"
#include "nuto/mechanics/elements/Plane2D.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include <assert.h>


NuTo::Plane2D::Plane2D(NuTo::StructureBase* rStructure, ElementData::eElementDataType rElementDataType,
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
void NuTo::Plane2D::CalculateGlobalRowDofs(std::vector<int>& rGlobalRowDofs) const
{
	rGlobalRowDofs.resize(2 * this->GetNumNodes());
    for (int nodeCount = 0; nodeCount < this->GetNumNodes(); nodeCount++)
    {
        const NodeBase *nodePtr = this->GetNode(nodeCount);
        if (nodePtr->GetNodeType()==Node::NodeCoordinatesDisplacementsMultiscale2D)
        {
            rGlobalRowDofs[2 * nodeCount    ] = nodePtr->GetDofFineScaleDisplacement(0);
            rGlobalRowDofs[2 * nodeCount + 1] = nodePtr->GetDofFineScaleDisplacement(1);
        }
        else
        {
            rGlobalRowDofs[2 * nodeCount    ] = nodePtr->GetDofDisplacement(0);
            rGlobalRowDofs[2 * nodeCount + 1] = nodePtr->GetDofDisplacement(1);
        }
    }
}

// build global column dof
void NuTo::Plane2D::CalculateGlobalColumnDofs(std::vector<int>& rGlobalColumnDofs) const
{
    int NumNonlocalElements(GetNumNonlocalElements());
	if (GetNumNonlocalElements()==0)
	    this->CalculateGlobalRowDofs(rGlobalColumnDofs);
    else
    {
	    const std::vector<const ElementBase*>& nonlocalElements(GetNonlocalElements());
        int NumCols(0);
        for (int theNonlocalElement=0; theNonlocalElement<NumNonlocalElements; theNonlocalElement++)
        {
        	NumCols += nonlocalElements[theNonlocalElement]->AsPlane()->GetNumLocalDofs();
        }
        rGlobalColumnDofs.resize(NumCols);
        int shift(0);
        for (int theNonlocalElement=0; theNonlocalElement<NumNonlocalElements; theNonlocalElement++)
        {
            const ElementBase* nonlocalElement(nonlocalElements[theNonlocalElement]);
        	for (int nodeCount = 0; nodeCount < nonlocalElement->GetNumNodes(); nodeCount++)
            {
                const NodeBase *nodePtr = nonlocalElement->GetNode(nodeCount);
                if (nodePtr->GetNodeType()==Node::NodeCoordinatesDisplacementsMultiscale2D)
                {
                    rGlobalColumnDofs[shift + 2 * nodeCount    ] = nodePtr->GetDofFineScaleDisplacement(0);
                    rGlobalColumnDofs[shift + 2 * nodeCount + 1] = nodePtr->GetDofFineScaleDisplacement(1);
                }
                else
                {
                    rGlobalColumnDofs[shift + 2 * nodeCount    ] = nodePtr->GetDofDisplacement(0);
                    rGlobalColumnDofs[shift + 2 * nodeCount + 1] = nodePtr->GetDofDisplacement(1);
                }
            }
            shift+=2*nonlocalElement->GetNumNodes();
        }
        assert(shift==NumCols);
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
//! @todo move to geometry class
bool NuTo::Plane2D::CheckPointInsidePolygon( const std::tuple<double,double> *rPoint, const std::vector<std::tuple<double,double> > * rPolygon)const
{
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

