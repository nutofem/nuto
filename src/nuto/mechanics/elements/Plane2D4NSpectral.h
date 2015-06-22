// $Id: Plane2D4NSpectralOrder2.h 647 2013-10-31 13:23:00Z unger3 $
#ifndef Plane2D4NSpectral_H
#define Plane2D4NSpectral_H

#ifdef ENABLE_SERIALIZATION
// serialize
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif  // ENABLE_SERIALIZATION

#include <vector>
#include "nuto/mechanics/elements/Plane2D.h"
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/elements/Truss1D2NSpectral.h"

namespace NuTo
{

template <int TOrder>
class Plane2D4NSpectral : public Plane2D
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION

public:
    Plane2D4NSpectral(NuTo::StructureBase* rStructure, const std::vector<NuTo::NodeBase* >& rNodes,
    		ElementData::eElementDataType rElementDataType, IpData::eIpDataType rIpDataType)
    : NuTo::Plane2D::Plane2D(rStructure, rElementDataType, GetStandardIntegrationType(),rIpDataType)
    {
		if ((int)rNodes.size()!=GetNumNodes())
		{
			std::cout << "rNodes.size() = " << rNodes.size() << std::endl;
			std::cout << "GetNumNodes() = " << GetNumNodes() << std::endl;
			throw MechanicsException("[NuTo::Plane2D4NSpectral::Plane2D4NSpectral] Number of nodes is (Order+1)^2.");
		}
		for (int count=0; count<GetNumNodes(); count++)
			mNodes[count] = rNodes[count];

		this->CheckElement(); //Geometry

		//check the position of the nodes on the edges
		if (GetNumIntegrationPoints()!=GetNumNodes())
		{
			throw MechanicsException("[NuTo::Plane2D4NSpectral::Plane2D4NSpectral] Number of nodes and integration points must be identical for spectral elements.");
		}

		//calculate coordinates
		int numLocalCoordinates(2*GetNumNodesGeometry());
		std::vector<double> localNodeCoord(numLocalCoordinates);
		CalculateLocalCoordinates(localNodeCoord);

		double naturalIPCoord[2];
		std::vector<double> shapeFunctionsGeometry(GetNumNodesGeometry());    //allocate space for shape functions
		for (int countNode=0; countNode<GetNumNodes(); countNode++)
		{
			GetLocalIntegrationPointCoordinates(countNode, naturalIPCoord);
			CalculateShapeFunctionsGeometry(naturalIPCoord, shapeFunctionsGeometry);

			double globalCalc[2];
			globalCalc[0] = 0.;
			globalCalc[1] = 0.;
			for (int count=0; count<GetNumNodesGeometry(); count++)
			{
				globalCalc[0] += shapeFunctionsGeometry[count]*localNodeCoord[2*count];
				globalCalc[1] += shapeFunctionsGeometry[count]*localNodeCoord[2*count+1];
			}

			//check the coordinates
			double globalReal[2];
			GetNode(countNode)->GetCoordinates2D(globalReal);
			double norm2 = (globalCalc[0]-globalReal[0])*(globalCalc[0]-globalReal[0])+
					(globalCalc[1]-globalReal[1])*(globalCalc[1]-globalReal[1]);
			if (norm2>1e-10)
			{
				std::cout <<"natural coordinates " << naturalIPCoord[0] << " " <<  naturalIPCoord[1] << std::endl;
				std::cout <<"local node coordinates 0 " << localNodeCoord[0] << " " <<  localNodeCoord[1] << std::endl;
				std::cout <<"local node coordinates 1 " << localNodeCoord[2] << " " <<  localNodeCoord[3] << std::endl;
				std::cout <<"local node coordinates 2 " << localNodeCoord[4] << " " <<  localNodeCoord[5] << std::endl;
				std::cout <<"local node coordinates 3 " << localNodeCoord[6] << " " <<  localNodeCoord[7] << std::endl;

				std::cout <<"shape functions " << shapeFunctionsGeometry[0] << " " <<  shapeFunctionsGeometry[1] << " " <<  shapeFunctionsGeometry[2]<< " " <<  shapeFunctionsGeometry[3]<< std::endl;
				std::cout <<"real global coordinates " << globalReal[0] << " " <<  globalReal[1] << std::endl;
				std::cout <<"calc global coordinates " << globalCalc[0] << " " <<  globalCalc[1] << std::endl;
				throw MechanicsException("[NuTo::Plane2D4NSpectral::Plane2D4NSpectral] The coordinates of local node " + std::to_string(countNode)
				+ " for the geometry interpolation are not correct. (error " + std::to_string(std::sqrt(norm2)) + ").");
			}
		}
/*		//check shape functions at the nodes
		std::vector<double> shapeFunctions(GetNumNodesField1D());    //allocate space for shape functions
		for (int count=0; count<GetNumNodesField1D(); count++)
		{
			double naturalIPCoord[2];
			GetLocalIntegrationPointCoordinates(count, naturalIPCoord);
			CalculateShapeFunctionsField1D(naturalIPCoord[0], shapeFunctions);
			std::cout << "shape function " << count << " with ";
			for (int count2=0; count2<GetNumNodesField1D(); count2++)
			{
				std::cout << shapeFunctions[count2] << " ";
			}
			std::cout << std::endl;
		}
*/
    }

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
	#ifdef DEBUG_SERIALIZATION
		std::cout << "start serialize Plane2D4NSpectral" << std::endl;
	#endif
		ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Plane2D)
		   & BOOST_SERIALIZATION_NVP(mNodes);
	#ifdef DEBUG_SERIALIZATION
		std::cout << "finish serialize Plane2D4NSpectral" << std::endl;
	#endif
    }
#endif // ENABLE_SERIALIZATION

    //! @brief returns the number of nodes in this element
    //! @return number of nodes
    inline int GetNumNodes()const
    {
        return (TOrder+1)*(TOrder+1);
    }

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @return pointer to the node
    inline NodeBase* GetNode(int rLocalNodeNumber)
    {
        assert(rLocalNodeNumber>=0 && rLocalNodeNumber<GetNumNodes());
        return mNodes[rLocalNodeNumber];
    }

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @return pointer to the node
    inline const NodeBase* GetNode(int rLocalNodeNumber)const
    {
    	assert(rLocalNodeNumber>=0 && rLocalNodeNumber<GetNumNodes());
        return mNodes[rLocalNodeNumber];
    }


    //! @brief sets the rLocalNodeNumber-th node of the element
    //! @param local node number
    //! @param pointer to the node
    inline void SetNode(int rLocalNodeNumber, NodeBase* rNode)
    {
    	assert(rLocalNodeNumber>=0 && rLocalNodeNumber<GetNumNodes());
        mNodes[rLocalNodeNumber] = rNode;
    }

    //! @brief returns the number of nodes in this element(geometry interpolation)
    //! @return number of nodes
    inline int GetNumNodesGeometry()const
    {
    	return 4;
    }

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @return pointer to the node
    inline const NodeBase* GetNodeGeometry(int rLocalNodeNumber)const
    {
        assert(rLocalNodeNumber>=0 && rLocalNodeNumber<4);
        switch(rLocalNodeNumber)
        {
        case 0:
        	return mNodes[0];
        	break;
        case 1:
        	return mNodes[TOrder];
        	break;
        case 2:
            return mNodes[GetNumNodes()-1];
            break;
        case 3:
        	return mNodes[GetNumNodes()-1-TOrder]; //numnodes-1-(TOrder)
        	break;
        default:
        	throw MechanicsException("[Plane2D4NSpectralOrder2::GetNodeGeometry] element has only 4 nodes for geometry interpolation.");
        }
    }

    //! @brief returns a pointer to the i-th node of the element (geometry interpolation)
    //! @param local node number
    //! @return pointer to the node
    inline NodeBase* GetNodeGeometry(int rLocalNodeNumber)
    {
        assert(rLocalNodeNumber>=0 && rLocalNodeNumber<4);
        switch(rLocalNodeNumber)
        {
        case 0:
        	return mNodes[0];
        	break;
        case 1:
        	return mNodes[2];
        	break;
        case 2:
            return mNodes[8];
            break;
        case 3:
        	return mNodes[6];
        	break;
        default:
        	throw MechanicsException("[Plane2D4NSpectralOrder2::GetNodeGeometry] element has only 4 nodes for geometry interpolation.");
        }
    }

    //number of nodes in one dimension (for 2D, it's a tensor product)
    inline int GetNumNodesField1D()const
    {
    	return TOrder+1;
    }

    //! @brief returns the number of nodes in this element(geometry interpolation)
    //! @return number of nodes
    inline int GetNumNodesField()const
    {
    	return GetNumNodes();
    }

    //! @brief returns a pointer to the i-th node of the element (geometry interpolation)
    //! @param local node number
    //! @return pointer to the node
    const NodeBase* GetNodeField(int rLocalNodeNumber)const
    {
        assert(rLocalNodeNumber>=0 && rLocalNodeNumber<GetNumNodes());
        return mNodes[rLocalNodeNumber];
    }

    //! @brief returns a pointer to the i-th node of the element (field interpolation)
    //! @param local node number
    //! @return pointer to the node
    NodeBase* GetNodeField(int rLocalNodeNumber)
    {
        assert(rLocalNodeNumber>=0 && rLocalNodeNumber<GetNumNodes());
        return mNodes[rLocalNodeNumber];
    }

    //! @brief returns the enum (type of the element)
    //! @return enum
    NuTo::Element::eElementType GetEnumType()const
    {
        switch (TOrder) {
        case 2:
            return NuTo::Element::PLANE2D4NSPECTRALORDER2;
        case 3:
            return NuTo::Element::PLANE2D4NSPECTRALORDER3;
        case 4:
            return NuTo::Element::PLANE2D4NSPECTRALORDER4;
        default:
             throw MechanicsException("[NuTo::Truss2D4NSpectral] The specified order is not implemented yet.");
            break;
        }
    }

    //! brief exchanges the node ptr in the full data set (elements, groups, loads, constraints etc.)
    //! this routine is used, if e.g. the data type of a node has changed, but the restraints, elements etc. are still identical
    //! @param rOldPtr old node ptr
    //! @param rNewPtr new node ptr
    void ExchangeNodePtr(NodeBase* rOldPtr, NodeBase* rNewPtr)
    {
        for (int count=0; count<GetNumNodes(); count++)
        {
            if (this->mNodes[count]==rOldPtr)
            {
                this->mNodes[count]=rNewPtr;
                break;
            }
        }
    }

    //! @brief calculates the shape functions
    //! @param rNaturalCoordinates natural coordinates of the integration point
    //! @param shape functions for all the nodes
    void CalculateShapeFunctionsGeometry(const double rNaturalCoordinates[2], std::vector<double>& rShapeFunctions)const
    {
        NuTo::ShapeFunctions2D::ShapeFunctionsQuadOrder1(rNaturalCoordinates, rShapeFunctions);
    }

    //! @brief calculates the derivative of the shape functions
    //! @param rNaturalCoordinates natural coordinates (-1,1) of the integration point
    //! @param derivative of the shape functions for all the nodes,
    //! first all the directions for a single node, and then for the next node
    void CalculateDerivativeShapeFunctionsGeometryNatural(const double rNaturalCoordinates[2], std::vector<double>& rDerivativeShapeFunctions)const
    {
        NuTo::ShapeFunctions2D::DerivativeShapeFunctionsQuadOrder1(rNaturalCoordinates, rDerivativeShapeFunctions);
    }

    //! @brief calculates the shape functions
    //! @param rNaturalCoordinates natural coordinates of the integration point
    //! @param shape functions for all the nodes
    void CalculateShapeFunctionsField(const double rNaturalCoordinates[2], std::vector<double>& rShapeFunctions)const
    {
		assert(((int)rShapeFunctions.size())==GetNumNodes());

        std::vector<double> shapeFunctions1Dx(GetNumNodesField1D());
		std::vector<double> shapeFunctions1Dy(GetNumNodesField1D());

        switch (TOrder) {
        case 2:
            NuTo::ShapeFunctions1D::ShapeFunctionsTrussOrder2(rNaturalCoordinates[0],shapeFunctions1Dx);
            NuTo::ShapeFunctions1D::ShapeFunctionsTrussOrder2(rNaturalCoordinates[1],shapeFunctions1Dy);
            break;
        case 3:
            NuTo::ShapeFunctions1D::ShapeFunctionsTrussSpectralOrder3(rNaturalCoordinates[0],shapeFunctions1Dx);
            NuTo::ShapeFunctions1D::ShapeFunctionsTrussSpectralOrder3(rNaturalCoordinates[1],shapeFunctions1Dy);
            break;
        case 4:
            NuTo::ShapeFunctions1D::ShapeFunctionsTrussSpectralOrder4(rNaturalCoordinates[0],shapeFunctions1Dx);
            NuTo::ShapeFunctions1D::ShapeFunctionsTrussSpectralOrder4(rNaturalCoordinates[1],shapeFunctions1Dy);
            break;
        default:
            break;
        }

		int theShapeFunction(0);
		for (int county=0; county<GetNumNodesField1D(); county++)
		{
			for (int countx=0; countx<GetNumNodesField1D(); countx++)
			{
				rShapeFunctions[theShapeFunction] = shapeFunctions1Dx[countx]*shapeFunctions1Dy[county];
				theShapeFunction++;
			}
		}
    }

    //! @brief calculates the 1D shape functions (2D is just the tensor product)
    //! @param rNaturalCoordinate natural coordinates of the integration point
    //! @param shape functions for all the nodes
    void CalculateShapeFunctionsField1D(double rNaturalCoordinate, std::vector<double>& rShapeFunctions)const;

    //! @brief calculates the derivative of the shape functions
    //! @param rNaturalCoordinates natural coordinates (-1,1) of the integration point
    //! @param derivative of the shape functions for all the nodes,
    //! first all the directions for a single node, and then for the next node
    void CalculateDerivativeShapeFunctionsFieldNatural(const double rNaturalCoordinates[2], std::vector<double>& rDerivativeShapeFunctions)const
    {
		assert(((int)rDerivativeShapeFunctions.size())==GetNumNodes()*2);

		std::vector<double> shapeFunctions1Dx(GetNumNodesField1D());
        std::vector<double> shapeFunctions1Dy(GetNumNodesField1D());
        std::vector<double> derShapeFunctions1Dx(GetNumNodesField1D());
		std::vector<double> derShapeFunctions1Dy(GetNumNodesField1D());

        switch (TOrder) {
        case 2:
            NuTo::ShapeFunctions1D::ShapeFunctionsTrussOrder2(rNaturalCoordinates[0],shapeFunctions1Dx);
            NuTo::ShapeFunctions1D::ShapeFunctionsTrussOrder2(rNaturalCoordinates[1],shapeFunctions1Dy);
            NuTo::ShapeFunctions1D::DerivativeShapeFunctionsTrussOrder2(rNaturalCoordinates[0],derShapeFunctions1Dx);
            NuTo::ShapeFunctions1D::DerivativeShapeFunctionsTrussOrder2(rNaturalCoordinates[1],derShapeFunctions1Dy);
            break;
        case 3:
            NuTo::ShapeFunctions1D::ShapeFunctionsTrussSpectralOrder3(rNaturalCoordinates[0],shapeFunctions1Dx);
            NuTo::ShapeFunctions1D::ShapeFunctionsTrussSpectralOrder3(rNaturalCoordinates[1],shapeFunctions1Dy);
            NuTo::ShapeFunctions1D::DerivativeShapeFunctionsTrussSpectralOrder3(rNaturalCoordinates[0],derShapeFunctions1Dx);
            NuTo::ShapeFunctions1D::DerivativeShapeFunctionsTrussSpectralOrder3(rNaturalCoordinates[1],derShapeFunctions1Dy);
            break;
        case 4:
            NuTo::ShapeFunctions1D::ShapeFunctionsTrussSpectralOrder4(rNaturalCoordinates[0],shapeFunctions1Dx);
            NuTo::ShapeFunctions1D::ShapeFunctionsTrussSpectralOrder4(rNaturalCoordinates[1],shapeFunctions1Dy);
            NuTo::ShapeFunctions1D::DerivativeShapeFunctionsTrussSpectralOrder4(rNaturalCoordinates[0],derShapeFunctions1Dx);
            NuTo::ShapeFunctions1D::DerivativeShapeFunctionsTrussSpectralOrder4(rNaturalCoordinates[1],derShapeFunctions1Dy);
            break;
        default:
            break;
        }

		int theDerShapeFunction(0);
		for (int county=0; county<GetNumNodesField1D(); county++)
		{
			for (int countx=0; countx<GetNumNodesField1D(); countx++)
			{
				rDerivativeShapeFunctions[theDerShapeFunction] = derShapeFunctions1Dx[countx]*shapeFunctions1Dy[county];
				theDerShapeFunction++;
				rDerivativeShapeFunctions[theDerShapeFunction] = shapeFunctions1Dx[countx]*derShapeFunctions1Dy[county];
				theDerShapeFunction++;
			}
		}
    }

    //! @brief calculates the derivative of the shape functions in 1D
    //! @param rNaturalCoordinate natural coordinates (-1,1) of the integration point
    //! @param derivative of the shape functions for all the nodes,
    //! first all the directions for a single node, and then for the next node
    void CalculateDerivativeShapeFunctionsFieldNatural1D(double rNaturalCoordinate, std::vector<double>& rDerivativeShapeFunctions)const;

    //! @brief calculates the shape functions for the surfaces (required for surface loads)
    //! @param rLocalCoordinates local coordinates of the integration point (in the local surface coordinate system)
    //! @param shape functions for all the nodes, size should already be correct, but can be checked with an assert
    void CalculateShapeFunctionsSurface(double rLocalCoordinates, std::vector<double>& rShapeFunctions)const override
    {
        assert(((int)rShapeFunctions.size())==GetNumNodesField1D());

        switch (TOrder) {
        case 2:
            NuTo::ShapeFunctions1D::ShapeFunctionsTrussOrder2(rLocalCoordinates,rShapeFunctions);
            break;
        case 3:
            NuTo::ShapeFunctions1D::ShapeFunctionsTrussSpectralOrder3(rLocalCoordinates,rShapeFunctions);
            break;
        case 4:
            NuTo::ShapeFunctions1D::ShapeFunctionsTrussSpectralOrder4(rLocalCoordinates,rShapeFunctions);
            break;
        default:
            break;
        }
    }

    //! @brief calculates the derivative of the shape functions with respect to local coordinatesfor the surfaces (required for surface loads)
    //! @param rLocalCoordinates local coordinates of the integration point
    //! @param derivative of the shape functions for all the nodes, size should already be correct, but can be checked with an assert
    //! first all the directions for a single node, and then for the next node
    void CalculateDerivativeShapeFunctionsLocalSurface(double rLocalCoordinates, std::vector<double>& rDerivativeShapeFunctions)const override
    {
        assert(((int)rDerivativeShapeFunctions.size())==GetNumNodesField1D());

        switch (TOrder) {
        case 2:
            NuTo::ShapeFunctions1D::DerivativeShapeFunctionsTrussOrder2(rLocalCoordinates,rDerivativeShapeFunctions);
            break;
        case 3:
            NuTo::ShapeFunctions1D::DerivativeShapeFunctionsTrussSpectralOrder3(rLocalCoordinates,rDerivativeShapeFunctions);
            break;
        case 4:
            NuTo::ShapeFunctions1D::DerivativeShapeFunctionsTrussSpectralOrder4(rLocalCoordinates,rDerivativeShapeFunctions);
            break;
        default:
            break;
        }
    }

    //! @brief returns the surface nodes
    //! @param surface (numbering so that the normal (right hand /thumb rule) is pointing outwards)
    //! @param surface nodes
    void GetSurfaceNodes(int rSurface, std::vector<const NodeBase*>& rSurfaceNodes)const override
    {
    	rSurfaceNodes.resize(GetNumNodesField1D());
    	switch (rSurface)
    	{
    	case 0:
			{
				for (int count=0; count<GetNumNodesField1D(); count++)
				{
					rSurfaceNodes[count] = mNodes[count];
				}
			}
		break;
    	case 1:
			{
				int theNode=GetNumNodesField1D()-1;
				for (int count=0; count<GetNumNodesField1D(); count++)
				{
					rSurfaceNodes[count] = mNodes[theNode];
					theNode+=GetNumNodesField1D();
				}
			}
			break;
    	case 2:
			{
				int theNode=GetNumNodesField()-1;
				for (int count=0; count<GetNumNodesField1D(); count++)
				{
					rSurfaceNodes[count] = mNodes[theNode];
					theNode--;
				}
			}
			break;
    	case 3:
			{
				int theNode=GetNumNodesField()-GetNumNodesField1D();
				for (int count=0; count<GetNumNodesField1D(); count++)
				{
					rSurfaceNodes[count] = mNodes[theNode];
                    theNode-=GetNumNodesField1D();
				}
			}
			break;
    	default:
    		throw MechanicsException("[Plane2D4NSpectral::GetSurfaceNodes] error.");
    	}
    }

    //! @brief returns the number of external surfaces
    //! @param surface (numbering so that the normal (right hand /thumb rule) is pointing outwards)
    //! @param surface nodes
    int GetNumSurfaces()const override
    {
        return 4;
    }

    int GetOrder()
    {
    	return TOrder;
    }

    //! @brief returns the enum of the standard integration type for this element
    NuTo::IntegrationType::eIntegrationType GetStandardIntegrationType()
    {
        switch (TOrder) {
        case 2:
            return NuTo::IntegrationType::IntegrationType2D4NLobatto9Ip;
        case 3:
            return NuTo::IntegrationType::IntegrationType2D4NLobatto16Ip;
        case 4:
            return NuTo::IntegrationType::IntegrationType2D4NLobatto25Ip;
        default:
            throw MechanicsException("[NuTo::Truss2D4NSpectral] The specified order is not implemented yet.");
            break;
        }
    }

protected:
    //! @brief ... just for serialization
    Plane2D4NSpectral(){}
    //! @brief ... reorder nodes such that the sign of the length of the element changes
    void ReorderNodes()
    {
    	throw MechanicsException("[NuTo::Plane2D4NSpectral::ReorderNodes] not implemented.");
    }

    //! @brief element nodes
    NodeBase* mNodes[(TOrder+1)*(TOrder+1)];
}; //class definition
} //namespace nuto

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::Plane2D4NSpectral<2>)
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::Plane2D4NSpectral<3>)
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::Plane2D4NSpectral<4>)
#ifndef SWIG
BOOST_CLASS_EXPORT_KEY(BOOST_IDENTITY_TYPE((NuTo::Plane2D4NSpectral<2>)))
BOOST_CLASS_EXPORT_KEY(BOOST_IDENTITY_TYPE((NuTo::Plane2D4NSpectral<3>)))
BOOST_CLASS_EXPORT_KEY(BOOST_IDENTITY_TYPE((NuTo::Plane2D4NSpectral<4>)))
#endif // SWIG
#endif // ENABLE_SERIALIZATION

#endif //PLANE2D4NSPECTRAL_H
