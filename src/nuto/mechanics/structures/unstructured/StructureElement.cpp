// $Id$

#include <assert.h>
#include <boost/tokenizer.hpp>
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/elements/IpDataEnum.h"
#include "nuto/mechanics/elements/Truss1D2N.h"
#include "nuto/mechanics/elements/BoundaryGradientDamage1D.h"
#include "nuto/mechanics/elements/Brick8N.h"
#include "nuto/mechanics/elements/Plane2D10N.h"
#include "nuto/mechanics/elements/Plane2D3N.h"
#include "nuto/mechanics/elements/Plane2D6N.h"
#include "nuto/mechanics/elements/Plane2D4N.h"
#include "nuto/mechanics/elements/Plane2D4NSpectral.h"
#include "nuto/mechanics/elements/Truss1D3N.h"
#include "nuto/mechanics/elements/Tetrahedron4N.h"
#include "nuto/mechanics/elements/Tetrahedron10N.h"
#include "nuto/mechanics/groups/Group.h"
#include "nuto/mechanics/nodes/NodeDof.h"
#include "nuto/mechanics/nodes/NodeCoordinates.h"
#include "nuto/mechanics/nodes/NodeCoordinatesDof.h"
//! @brief returns the number of nodes
//! @return number of nodes
int NuTo::Structure::GetNumElements() const
{
    return mElementMap.size();
}

//! @brief a reference to an element
//! @param identifier
//! @return reference to an element
NuTo::ElementBase* NuTo::Structure::ElementGetElementPtr(int rIdent)
{
    boost::ptr_map<int,ElementBase>::iterator it = mElementMap.find(rIdent);
    if (it!=mElementMap.end())
        return it->second;
    else
    {
    	std::stringstream message;
    	message << "[NuTo::Structure::ElementGetElementPtr] Element with identifier " << rIdent << " does not exist." << std::endl;
    	throw MechanicsException(message.str());
    }
}

//! @brief a reference to an element
//! @param identifier
//! @return reference to an element
const NuTo::ElementBase* NuTo::Structure::ElementGetElementPtr(int rIdent)const
{
    boost::ptr_map<int,ElementBase>::const_iterator it = mElementMap.find(rIdent);
    if (it!=mElementMap.end())
        return it->second;
    else
    {
    	std::stringstream message;
    	message << "[NuTo::Structure::ElementGetElementPtr] Element with identifier " << rIdent << " does not exist." << std::endl;
    	throw MechanicsException(message.str());
    }
}

//! @brief gives the identifier of an element
//! @param reference to an element
//! @return element number
int NuTo::Structure::ElementGetId(const ElementBase* rElement)const
{
    for (boost::ptr_map<int,ElementBase>::const_iterator
            it = mElementMap.begin(); it!= mElementMap.end(); it++)
    {
        if (it->second==rElement)
            return it->first;
    }
    throw MechanicsException("[NuTo::Structure::GetElementId] Element does not exist.");
}

//! @brief info about one single element
//! @param rElement (Input) ... pointer to the element
//! @param rVerboseLevel (Input) ... level of verbosity
void NuTo::Structure::ElementInfo(const ElementBase* rElement, int rVerboseLevel)const
{
	std::cout << "element : " << rElement->ElementGetId() << std::endl;
	if (rVerboseLevel>2)
	{
		std::cout << "\tenum::type=" << rElement->GetEnumType() << std::endl;
		if (rVerboseLevel>3)
		{
			std::cout << "\tNodes:";
			for(unsigned short iNode=0; iNode<rElement->GetNumNodes(); ++iNode)
				std::cout << "\t" << this->NodeGetId(rElement->GetNode(iNode));
			std::cout << std::endl;
			if (rVerboseLevel>4)
			{
				std::cout << "\tintegration points :" << std::endl;
				for(int iIp=0; iIp<rElement->GetNumIntegrationPoints(); ++iIp)
				{
					double coor[3];
					rElement->GetGlobalIntegrationPointCoordinates(iIp,coor);
					std::cout << "\t\t" << iIp << ": [" << coor[0] << ";" << coor[1] << ";" << coor[2] << "]" << std::endl;
				}
			}
		}
	}
}

//! @brief info about the elements in the Structure
void NuTo::Structure::ElementInfo(int rVerboseLevel)const
{
    std::cout<<"number of elements: " << mElementMap.size() <<std::endl;
    if (rVerboseLevel>3)
    {
    	std::cout << "\t\telements :" <<std::endl;
    	for (boost::ptr_map<int,ElementBase>::const_iterator it = mElementMap.begin(); it!= mElementMap.end(); it++)
    	{
    		std::cout << "\t\t" << it->first;
    	    if (rVerboseLevel>4)
    	    {
				std::cout << "\t:";
				for(unsigned short iNode=0; iNode<it->second->GetNumNodes(); ++iNode)
					std::cout << "\t" << this->NodeGetId(it->second->GetNode(iNode));
    	    }
        	std::cout << std::endl;
    	}
    }
}



//! @brief Creates an element
//! @param rElementIdent identifier for the element
//! @param rElementType element type
//! @param rNodeIdents Identifier for the corresponding nodes
int NuTo::Structure::ElementCreate (const std::string& rElementType,
        const NuTo::FullVector<int,Eigen::Dynamic>& rNodeNumbers)
{
	return ElementCreate(rElementType,rNodeNumbers,std::string("CONSTITUTIVELAWIP"),std::string("NOIPDATA") );
}


//! @param rGroupNumberElements group for elements on the real boundary
//! @param rGroupNumberBoundaryNodes nodes on the boundary
//! @param rElementType element type
//! @param rNodeIdents Identifier for the corresponding nodes
void NuTo::Structure::BoundaryElementsCreate (const std::string& rElementType,
		int rGroupNumberElements, int rGroupNumberBoundaryNodes,
		int rOrder, double rVirtualBoundary,
		const std::string& rElementDataType, const std::string& rIpDataType)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    // get element type
    std::string upperCaseElementType;
    std::transform(rElementType.begin(), rElementType.end(), std::back_inserter(upperCaseElementType), (int(*)(int)) toupper);

    Element::eElementType elementType;
    if (upperCaseElementType=="BOUNDARYGRADIENTDAMAGE1D")
    {
    	elementType = NuTo::Element::BOUNDARYGRADIENTDAMAGE1D;
    }
    else
    {
    	throw MechanicsException("[NuTo::Structure::BoundaryElementsCreate] Element type "+upperCaseElementType +" does not exist.");
    }

    //find groups
	boost::ptr_map<int,GroupBase>::iterator itGroupElements = mGroupMap.find(rGroupNumberElements);
    if (itGroupElements==mGroupMap.end())
        throw MechanicsException("[NuTo::Structure::BoundaryElementsCreate] Group with the given identifier does not exist.");
    if (itGroupElements->second->GetType()!=NuTo::Groups::Elements)
    	throw MechanicsException("[NuTo::Structure::BoundaryElementsCreate] Group is not an element group.");

	boost::ptr_map<int,GroupBase>::iterator itGroupBoundaryNodes = mGroupMap.find(rGroupNumberBoundaryNodes);
    if (itGroupBoundaryNodes==mGroupMap.end())
        throw MechanicsException("[NuTo::Structure::BoundaryElementsCreate] Group with the given identifier does not exist.");
    if (itGroupBoundaryNodes->second->GetType()!=NuTo::Groups::Nodes)
    	throw MechanicsException("[NuTo::Structure::BoundaryElementsCreate] Group is not a node group.");


    // check element data
    std::string upperCaseElementDataType;
    std::transform(rElementDataType.begin(), rElementDataType.end(), std::back_inserter(upperCaseElementDataType), (int(*)(int)) toupper);

    ElementData::eElementDataType elementDataType;
    if (upperCaseElementDataType=="CONSTITUTIVELAWIP")
    {
    	elementDataType = NuTo::ElementData::CONSTITUTIVELAWIP;
    }
    else
	{
		throw MechanicsException("[NuTo::Structure::BoundaryElementsCreate] Element data type "+upperCaseElementDataType +" does not exist.");
	}

    // check ip data
    std::string upperCaseIpDataType;
    std::transform(rIpDataType.begin(), rIpDataType.end(), std::back_inserter(upperCaseIpDataType), (int(*)(int)) toupper);

    IpData::eIpDataType ipDataType;
    if (upperCaseIpDataType=="STATICDATA")
   	{
    	ipDataType = NuTo::IpData::STATICDATA;
   	}
	else
	{
		throw MechanicsException("[NuTo::Structure::BoundaryElementsCreate] Ip data type "+upperCaseIpDataType +" does not exist.");
	}

    BoundaryElementsCreate(elementType,
    		itGroupElements->second->AsGroupElement(), itGroupBoundaryNodes->second->AsGroupNode(),
    		rOrder, rVirtualBoundary,
    		elementDataType, ipDataType);

#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::Structure::BoundaryElementsCreate] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}

//! @brief Create boundary elements defined by all boundary elements and the nodes characterizing the edges
//! @param rGroupNumberElements group for elements on the real boundary
//! @param rGroupNumberBoundaryNodes nodes on the boundary
//! @param rElementType element type
//! @param rNodeIdents Identifier for the corresponding nodes
void NuTo::Structure::BoundaryElementsCreate (Element::eElementType rType,
		const Group<ElementBase>* rGroupElements, const Group<NodeBase>* rGroupBoundaryNodes,
		int rOrder, double rVirtualBoundary,
		ElementData::eElementDataType rElementDataType, IpData::eIpDataType rIpDataType)
{
	double deltaL(rVirtualBoundary/rOrder);

    for (Group<ElementBase>::const_iterator itElement=rGroupElements->begin(); itElement!=rGroupElements->end();itElement++)
    {
        try
        {
        	std::map<NodeBase*,std::vector<NodeBase*> > copiedNodes;

        	//better check local dimension, but this is in general identical
    		switch(itElement->second->GetGlobalDimension())
    		{
    		case 1:
    		{
    			//in 1D it is either the first or the last node of the element
      			NodeBase* nodePtr [2];
      			nodePtr[0] = itElement->second->GetNode(0);
      			nodePtr[1] = itElement->second->GetNode(itElement->second->GetNumNodes()-1);

       			int nodeNumber[2];
       			double coordinate[2];

       			for (int count=0; count<2; count++)
       			{
       				nodeNumber[count] = NodeGetId(nodePtr[count]);
       				nodePtr[count]->GetCoordinates1D(&(coordinate[count]));
       			}

       			//loop over first/last node of element
       			for (int countNode=0; countNode<2; countNode++)
       			{
					if (rGroupBoundaryNodes->Contain(nodeNumber[countNode]))
					{
						int numTimeDerivatives = nodePtr[countNode]->GetNumTimeDerivatives();
						int numNonlocalNonlocalTotalStrain = nodePtr[countNode]->GetNumNonlocalTotalStrain();
						//check, if the virtual boundary extends to the left or to the right
						std::vector<NodeBase*> newNodes(rOrder+1);
						newNodes[0] = nodePtr[countNode];
						for (int countOrder=0; countOrder<rOrder; countOrder++)
						{
							NodeBase* newNodePtr(0);
							double NewCoordinate;
							if (countNode==0)
							{
								if (coordinate[0]<coordinate[1])
									//boundary extends to the left
									NewCoordinate = coordinate[0]-(countOrder+1)*deltaL;
								else
									//boundary extends to the right
									NewCoordinate = coordinate[0]+(countOrder+1)*deltaL;
							}
							else
							{
								if (coordinate[0]<coordinate[1])
									//boundary extends to the right
									NewCoordinate = coordinate[1]+(countOrder+1)*deltaL;
								else
									//boundary extends to the right
									NewCoordinate = coordinate[1]-(countOrder+1)*deltaL;
							}

							switch(numTimeDerivatives)
							{
							case 0:
								if (numNonlocalNonlocalTotalStrain==1)
									newNodePtr = new NuTo::NodeCoordinatesDof<1,0,0,0,0,0,1>();
								break;
							case 2:
								if (numNonlocalNonlocalTotalStrain==1)
									newNodePtr = new NuTo::NodeCoordinatesDof<1,1,0,0,0,0,1>();
								break;
							default:
								throw MechanicsException("[BoundaryElementsCreate] time derivative should be either 0 or 2.");
							}

							newNodePtr->SetCoordinates1D(&NewCoordinate);

							//find unused integer id
							int newNodeId(mNodeMap.size());
							boost::ptr_map<int,NodeBase>::iterator it = mNodeMap.find(newNodeId);
							while (it!=mNodeMap.end())
							{
								newNodeId++;
								it = mNodeMap.find(newNodeId);
							}
							// add new node to map
							this->mNodeMap.insert(newNodeId, newNodePtr);

							newNodes[countOrder+1] = newNodePtr;
						}
						//add to map (copied nodes with the original boundary node
						//this has to be done for 2D and 3D, but for 1D, boundary nodes belong only to a single element
						//copiedNodes[nodePtr1] = newNodes;

						//create the boundary element
						ElementBase* newElementPtr = new NuTo::BoundaryGradientDamage1D(this,newNodes,
								itElement->second->AsTruss(),countNode,
								NuTo::ElementData::CONSTITUTIVELAWIP,
								IntegrationType::IntegrationType1D2NBoundaryGauss3Ip,
								IpData::STATICDATA);

						//find unused integer id
						int elementNumber(mElementMap.size());
						boost::ptr_map<int,ElementBase>::iterator it = mElementMap.find(elementNumber);
						while (it!=mElementMap.end())
						{
							elementNumber++;
							it = mElementMap.find(elementNumber);
						}

						mElementMap.insert(elementNumber, newElementPtr);

						//modify the material of the element, use the constitutive law of the first ip in the neigboring element
						newElementPtr->SetConstitutiveLaw(itElement->second->GetConstitutiveLaw(0));


					}//if (rGroupBoundaryNodes->Contain(nodeNumber1))
				}//for (int count=0; count<2; count++)
    		}
    		break;
    		case 2:
    			// determine a list with the extended nodes that determine the virtual boundary elements
    			throw MechanicsException("[NuTo::Structure::BoundaryElementsCreate] not yet implemented for 2D.");
   			break;
    		case 3:
    			// determine a list with the extended nodes that determine the virtual boundary elements
    			throw MechanicsException("[NuTo::Structure::BoundaryElementsCreate] not yet implemented for 3D.");
    		break;
    		default:
    			throw MechanicsException("[NuTo::Structure::BoundaryElementsCreate] Dimension of the real boundary elements is neither 1,2 or 3.");
    		}

        	this->mNodeNumberingRequired  = true;
        }
        catch(NuTo::MechanicsException &e)
        {
            std::stringstream ss;
            assert(ElementGetId(itElement->second)==itElement->first);
            ss << itElement->first;
            e.AddMessage("[NuTo::StructureBase::BoundaryElementsCreate] Error creating boundary element for element "
            	+ ss.str() + ".");
            throw e;
        }
        catch(...)
        {
            std::stringstream ss;
            assert(ElementGetId(itElement->second)==itElement->first);
            ss << itElement->first;
        	throw NuTo::MechanicsException
        	   ("[NuTo::StructureBase::BoundaryElementsCreate] Error creating boundary element for element " + ss.str() + ".");
        }
    }




}




//! @brief Creates an element
//! @param rElementIdent identifier for the element
//! @param rElementType element type
//! @param rNodeIdents Identifier for the corresponding nodes
int NuTo::Structure::ElementCreate (const std::string& rElementType,
        const NuTo::FullVector<int,Eigen::Dynamic>& rNodeNumbers, const std::string& rElementDataType, const std::string& rIpDataType)
{
	//find unused integer id
	int elementNumber(mElementMap.size());
	boost::ptr_map<int,ElementBase>::iterator it = mElementMap.find(elementNumber);
	while (it!=mElementMap.end())
	{
		elementNumber++;
		it = mElementMap.find(elementNumber);
	}

	// create element
	this->ElementCreate(elementNumber, rElementType, rNodeNumbers,rElementDataType,rIpDataType);

	// return element number
	return elementNumber;
}
void NuTo::Structure::ElementCreate (int rElementNumber, const std::string& rElementType,
        const NuTo::FullVector<int,Eigen::Dynamic>& rNodeNumbers)
{
	ElementCreate(rElementNumber,rElementType,rNodeNumbers,std::string("CONSTITUTIVELAWIP"),std::string("NOIPDATA"));
}

void NuTo::Structure::ElementCreate (int rElementNumber, const std::string& rElementType,
        const NuTo::FullVector<int,Eigen::Dynamic>& rNodeNumbers, const std::string& rElementDataType, const std::string& rIpDataType)
{
	// check node number
	boost::ptr_map<int,ElementBase>::iterator it = mElementMap.find(rElementNumber);
	if(it != this->mElementMap.end())
	{
    	throw MechanicsException("[NuTo::Structure::ElementCreate] Element already exists.");
	}

	// convert node numbers to pointer
	if (rNodeNumbers.GetNumColumns()!=1)
    	throw MechanicsException("[NuTo::Structure::ElementCreate] Matrix with node numbers should have a single column.");
	std::vector<NodeBase*> nodeVector;
    for (int count=0; count<rNodeNumbers.GetNumRows(); count++)
    {
        nodeVector.push_back(NodeGetNodePtr(rNodeNumbers(count,0)));
    }

    //! @brief Returns the internal enum (number) for the element types
    //! @param Element name in Nuto
    //! @return enum
    Element::eElementType elementType = this->ElementTypeGetEnum(rElementType);

    // check element data
    std::string upperCaseElementDataType;
    std::transform(rElementDataType.begin(), rElementDataType.end(), std::back_inserter(upperCaseElementDataType), (int(*)(int)) toupper);

    ElementData::eElementDataType elementDataType;
    if (upperCaseElementDataType=="CONSTITUTIVELAWIP")
    {
    	elementDataType = NuTo::ElementData::CONSTITUTIVELAWIP;
    }
    else if (upperCaseElementDataType=="CONSTITUTIVELAWIPCRACK")
   	{
		elementDataType = NuTo::ElementData::CONSTITUTIVELAWIPCRACK;
   	}
    else if (upperCaseElementDataType=="CONSTITUTIVELAWIPNONLOCAL")
   	{
		elementDataType = NuTo::ElementData::CONSTITUTIVELAWIPNONLOCAL;
   	}
    else
	{
		throw MechanicsException("[NuTo::Structure::ElementCreate] Element data type "+upperCaseElementDataType +" does not exist.");
	}

    // check ip data
    std::string upperCaseIpDataType;
    std::transform(rIpDataType.begin(), rIpDataType.end(), std::back_inserter(upperCaseIpDataType), (int(*)(int)) toupper);

    IpData::eIpDataType ipDataType;
    if (upperCaseIpDataType=="NOIPDATA")
    {
    	ipDataType = NuTo::IpData::NOIPDATA;
    }
    else if (upperCaseIpDataType=="STATICDATA")
   	{
    	ipDataType = NuTo::IpData::STATICDATA;
   	}
    else if (upperCaseIpDataType=="STATICDATANONLOCAL")
   	{
    	ipDataType = NuTo::IpData::STATICDATANONLOCAL;
   	}
	else
	{
		throw MechanicsException("[NuTo::Structure::ElementCreate] Ip data type "+upperCaseIpDataType +" does not exist.");
	}

    // create element
    this->ElementCreate(rElementNumber, elementType, nodeVector, elementDataType, ipDataType);

}

//! @brief Returns the internal enum (number) for the element types
//! @param Element name in Nuto
//! @return enum
NuTo::Element::eElementType NuTo::Structure::ElementTypeGetEnum(const std::string& rElementType)
{
	// get element type
	std::string upperCaseElementType;
	std::transform(rElementType.begin(), rElementType.end(), std::back_inserter(upperCaseElementType), (int(*)(int)) toupper);

	Element::eElementType elementType;
	if (upperCaseElementType=="TRUSS1D2N")
	{
		elementType = NuTo::Element::TRUSS1D2N;
	}
	else if (upperCaseElementType=="TRUSS1D3N")
	{
		elementType = NuTo::Element::TRUSS1D3N;
	}
	else if (upperCaseElementType=="BRICK8N")
	{
		elementType = NuTo::Element::BRICK8N;
	}
	else if (upperCaseElementType=="PLANE2D10N")
	{
		elementType = NuTo::Element::PLANE2D10N;
	}
	else if (upperCaseElementType=="PLANE2D3N")
	{
		elementType = NuTo::Element::PLANE2D3N;
	}
	else if (upperCaseElementType=="PLANE2D4N")
	{
		elementType = NuTo::Element::PLANE2D4N;
	}
	else if (upperCaseElementType=="PLANE2D4NSPECTRALORDER2")
	{
		elementType = NuTo::Element::PLANE2D4NSPECTRALORDER2;
	}
	else if (upperCaseElementType=="PLANE2D4NSPECTRALORDER3")
	{
		elementType = NuTo::Element::PLANE2D4NSPECTRALORDER3;
	}
	else if (upperCaseElementType=="PLANE2D4NSPECTRALORDER4")
	{
		elementType = NuTo::Element::PLANE2D4NSPECTRALORDER4;
	}
	else if (upperCaseElementType=="PLANE2D6N")
	{
		elementType = NuTo::Element::PLANE2D6N;
	}
	else if (upperCaseElementType=="TETRAHEDRON4N")
	{
		elementType = NuTo::Element::TETRAHEDRON4N;
	}
	else if (upperCaseElementType=="TETRAHEDRON10N")
	{
		elementType = NuTo::Element::TETRAHEDRON10N;
	}
	else
	{
		throw MechanicsException("[NuTo::Structure::ElementCreate] Element type "+upperCaseElementType +" does not exist.");
	}
	return elementType;
}

//! @brief Creates an element
//! @param rElementNumber element number
//! @param rElementType element type
//! @param rNodeIdents pointers to the corresponding nodes
//! @return int rElementNumber
int NuTo::Structure::ElementCreate(Element::eElementType rType,
        const std::vector<NodeBase*>& rNodeVector, ElementData::eElementDataType rElementDataType, IpData::eIpDataType rIpDataType)
{
	//find unused integer id
	int elementNumber(mElementMap.size());
	boost::ptr_map<int,ElementBase>::iterator it = mElementMap.find(elementNumber);
	while (it!=mElementMap.end())
	{
		elementNumber++;
		it = mElementMap.find(elementNumber);
	}

	// create element
	this->ElementCreate(elementNumber,rType,rNodeVector,rElementDataType,rIpDataType);

	// return element number
	return elementNumber;
}
//! @brief Creates an element
//! @param rElementIdent identifier for the element
//! @param rElementType element type
//! @param rNodeIdents pointers to the corresponding nodes
//! @return element number
void NuTo::Structure::ElementCreate(int rElementNumber, Element::eElementType rType,
        const std::vector<NodeBase*>& rNodeVector, ElementData::eElementDataType rElementDataType, IpData::eIpDataType rIpDataType)
{
    try
    {
		ElementBase* ptrElement(0);
		switch (rType)
		{
		case NuTo::Element::TRUSS1D2N:
			if (1!=mDimension)
				throw MechanicsException("[NuTo::Structure::ElementCreate] TRUSS1D2N is only a 1D element, either change the dimension of the structure to one or use TRUSS3D2N.");
			ptrElement = new NuTo::Truss1D2N(this, rNodeVector, rElementDataType, rIpDataType);
			break;
		case NuTo::Element::TRUSS1D3N:
			if (1!=mDimension)
				throw MechanicsException("[NuTo::Structure::ElementCreate] TRUSS1D3N is only a 1D element, either change the dimension of the structure to one or use TRUSS3D3N.");
			ptrElement = new NuTo::Truss1D3N(this, rNodeVector, rElementDataType, rIpDataType);
			break;
		case NuTo::Element::BRICK8N:
			if (this->mDimension != 3)
			{
				throw MechanicsException("[NuTo::Structure::ElementCreate] Brick8N is a 3D element.");
			}
			ptrElement = new NuTo::Brick8N(this, rNodeVector, rElementDataType, rIpDataType);
			break;
		case NuTo::Element::PLANE2D10N:
			if (this->mDimension != 2)
			{
				throw MechanicsException("[NuTo::Structure::ElementCreate] PLANE2D10N is a 2D element.");
			}
			ptrElement = new NuTo::Plane2D10N(this, rNodeVector, rElementDataType, rIpDataType);
			break;
		case NuTo::Element::PLANE2D3N:
			if (this->mDimension != 2)
			{
				throw MechanicsException("[NuTo::Structure::ElementCreate] PLANE2D3N is a 2D element.");
			}
			ptrElement = new NuTo::Plane2D3N(this, rNodeVector, rElementDataType, rIpDataType);
			break;
		case NuTo::Element::PLANE2D4N:
			if (this->mDimension != 2)
			{
				throw MechanicsException("[NuTo::Structure::ElementCreate] PLANE2D4N is a 2D element.");
			}
			ptrElement = new NuTo::Plane2D4N(this, rNodeVector, rElementDataType, rIpDataType);
			break;
		case NuTo::Element::PLANE2D4NSPECTRALORDER2:
			if (this->mDimension != 2)
			{
				throw MechanicsException("[NuTo::Structure::ElementCreate] PLANE2D4NSPECTRALORDER2 is a 2D element.");
			}
			ptrElement = new NuTo::Plane2D4NSpectral<2>(this, rNodeVector, rElementDataType, rIpDataType);
			break;
		case NuTo::Element::PLANE2D4NSPECTRALORDER3:
			if (this->mDimension != 2)
			{
				throw MechanicsException("[NuTo::Structure::ElementCreate] PLANE2D4NSPECTRALORDER2 is a 2D element.");
			}
			ptrElement = new NuTo::Plane2D4NSpectral<3>(this, rNodeVector, rElementDataType, rIpDataType);
			break;
		case NuTo::Element::PLANE2D4NSPECTRALORDER4:
			if (this->mDimension != 2)
			{
				throw MechanicsException("[NuTo::Structure::ElementCreate] PLANE2D4NSPECTRALORDER2 is a 2D element.");
			}
			ptrElement = new NuTo::Plane2D4NSpectral<4>(this, rNodeVector, rElementDataType, rIpDataType);
			break;
		case NuTo::Element::PLANE2D6N:
			if (this->mDimension != 2)
			{
				throw MechanicsException("[NuTo::Structure::ElementCreate] PLANE2D6N is a 2D element.");
			}
			ptrElement = new NuTo::Plane2D6N(this, rNodeVector, rElementDataType, rIpDataType);
			break;
		case NuTo::Element::TETRAHEDRON4N:
			if (this->mDimension != 3)
			{
				throw MechanicsException("[NuTo::Structure::ElementCreate] Tetrahedron4N is a 3D element.");
			}
			ptrElement = new NuTo::Tetrahedron4N(this, rNodeVector, rElementDataType, rIpDataType);
			break;
		case NuTo::Element::TETRAHEDRON10N:
			if (this->mDimension != 3)
			{
				throw MechanicsException("[NuTo::Structure::ElementCreate] Tetrahedron10N is a 3D element.");
			}
			ptrElement = new NuTo::Tetrahedron10N(this, rNodeVector, rElementDataType, rIpDataType);
			break;
		default:
			throw NuTo::MechanicsException("[NuTo::Structure::ElementCreate] Invalid element type.");
		}
		mElementMap.insert(rElementNumber, ptrElement);
    }
    catch(NuTo::MechanicsException &e)
    {
        std::stringstream ss;
        ss << rElementNumber;
        e.AddMessage("[NuTo::Structure::ElementCreate] Error creating element "
        	+ ss.str() + ".");
        throw e;
    }
    catch(...)
    {
        std::stringstream ss;
        ss << rElementNumber;
    	throw NuTo::MechanicsException
    	   ("[NuTo::Structure::ElementCreate] Error creating element " + ss.str() + ".");
    }
}

//! @brief creates multiple elements
//! @param rElementType element type
//! @param rNodeIdents Identifier for the corresponding nodes (Incidences have to be stored column-wise)
//! @return a NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic> containing the element numbers
NuTo::FullVector<int,Eigen::Dynamic> NuTo::Structure::ElementsCreate (const std::string& rElementType, NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic> & rNodeNumbers)
{
	std::vector<int> idVec;
	/// go through the elements
	for(size_t i=0 ; i<(size_t)rNodeNumbers.GetNumColumns(); ++i)
	{
		const NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic> incidence(rNodeNumbers.GetColumn(i));
		idVec.push_back(this-> ElementCreate(rElementType,incidence ));
	}

    //return int identifiers of the new elements as FullMatrix
	NuTo::FullVector<int,Eigen::Dynamic> ids(idVec);
    return ids;
}

//! @brief creates multiple elements
//! @param rElementType element type
//! @param rNodeIdents Identifier for the corresponding nodes (Incidences have to be stored column-wise)
//! @param rElementDataType Element data for the elements
//! @param rIpDataType Integration point data for the elements
//! @return a NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic> containing the element numbers
NuTo::FullVector<int,Eigen::Dynamic> NuTo::Structure::ElementsCreate (const std::string& rElementType, NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic> & rNodeNumbers, const std::string& rElementDataType, const std::string& rIpDataType)
{
	std::vector<int> idVec;
	/// go through the elements
	for(size_t i=0 ; i<(size_t)rNodeNumbers.GetNumColumns(); ++i)
	{
		const NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic> incidence(rNodeNumbers.GetColumn(i));
		idVec.push_back(this-> ElementCreate(rElementType, incidence, rElementDataType, rIpDataType));
	}

    //return int identifiers of the new elements as FullMatrix
	NuTo::FullVector<int,Eigen::Dynamic> ids(idVec);
    return ids;
}

//! @brief Deletes an element
//! @param rElementIdent identifier for the element
void NuTo::Structure::ElementDelete(int rElementNumber)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
   	ElementDeleteInternal(rElementNumber);
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::Structure::ElementDelete] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}

//! @brief Deletes a group of elements element
//! @param rGroupNumber group number
void NuTo::Structure::ElementGroupDelete (int rGroupNumber, bool deleteNodes)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
	boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(rGroupNumber);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::ElementGroupDelete] Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::Groups::Elements)
    	throw MechanicsException("[NuTo::StructureBase::ElementGroupDelete] Group is not an element group.");

    //the group has to be copied, since the elements are removed from this group, which invalidates the iterators
    Group<ElementBase> copyOfElementGroup = *(dynamic_cast<Group<ElementBase>*>(itGroup->second));

    std::set<NodeBase*> potentialNodesToBeRemoved;
    for (Group<ElementBase>::iterator itElement=copyOfElementGroup.begin(); itElement!=copyOfElementGroup.end();itElement++)
    {
        try
        {
        	//save the nodes, which are eventually to be removed
        	if (deleteNodes)
        	{
        		for (int countNode=0; countNode<itElement->second->GetNumNodes(); countNode++)
        		{
        			NodeBase* nodePtr = itElement->second->GetNode(countNode);
        			potentialNodesToBeRemoved.insert(nodePtr);

        		}
        	}
        	ElementDeleteInternal(itElement->second->ElementGetId());
        }
        catch(NuTo::MechanicsException &e)
        {
            std::stringstream ss;
            assert(ElementGetId(itElement->second)==itElement->first);
            ss << itElement->first;
            e.AddMessage("[NuTo::StructureBase::ElementGroupDelete] Error deleting element "
            	+ ss.str() + ".");
            throw e;
        }
        catch(...)
        {
            std::stringstream ss;
            assert(ElementGetId(itElement->second)==itElement->first);
            ss << itElement->first;
        	throw NuTo::MechanicsException
        	   ("[NuTo::StructureBase::ElementGroupDelete] Error deleting element " + ss.str() + ".");
        }
    }

    //check all the other elements and see, if they have one of the potential Nodes To Be Removed as valid node
    for (boost::ptr_map<int,ElementBase>::iterator itElement = mElementMap.begin(); itElement != mElementMap.end(); itElement++)
    {
    	for (int countNode=0; countNode<(itElement->second)->GetNumNodes(); countNode++)
		{
			NodeBase* nodePtr = (itElement->second)->GetNode(countNode);
			//int numRemoved = potentialNodesToBeRemoved.erase(nodePtr);
			potentialNodesToBeRemoved.erase(nodePtr);
		}
    }

    for (std::set<NodeBase*>::iterator itNode = potentialNodesToBeRemoved.begin(); itNode != potentialNodesToBeRemoved.end(); itNode++)
    {
    	int nodeId(NodeGetId(*itNode));
    	NodeDelete(nodeId,false);
    }
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::ElementGroupDelete] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}

//! @brief Deletes an element
//! @param rItElement iterator of the map
void NuTo::Structure::ElementDeleteInternal(int rElementId)
{
	// find element
	boost::ptr_map<int,ElementBase>::iterator itElement = mElementMap.find(rElementId);
    if (itElement == this->mElementMap.end())
    {
        throw MechanicsException("[NuTo::Structure::ElementDeleteInternal] Element does not exist.");
    }
    else
    {
		// Search for elements in groups: using a loop over all groups
		for(boost::ptr_map<int,GroupBase>::iterator groupIt=mGroupMap.begin();groupIt!=mGroupMap.end(); ++groupIt)
		{
			if(groupIt->second->GetType()==NuTo::Groups::Elements)
			{
				if(groupIt->second->Contain(rElementId))
				{
					groupIt->second->RemoveMember(rElementId);
				}
			}
		}

		// delete element from map
		this->mElementMap.erase(itElement);
    }
}

// store all elements of a structure in a vector
void NuTo::Structure::GetElementsTotal(std::vector<const ElementBase*>& rElements) const
{
    rElements.reserve(mElementMap.size());
    rElements.resize(0);
	boost::ptr_map<int,ElementBase>::const_iterator ElementIter = this->mElementMap.begin();
    while (ElementIter != this->mElementMap.end())
    {
        rElements.push_back(ElementIter->second);
        ElementIter++;
    }
}

// store all elements of a structure in a vector
void NuTo::Structure::GetElementsTotal(std::vector<std::pair<int, const ElementBase*> >& rElements) const
{
	rElements.reserve(mElementMap.size());
	rElements.resize(0);
	boost::ptr_map<int,ElementBase>::const_iterator ElementIter = this->mElementMap.begin();
    while (ElementIter != this->mElementMap.end())
    {
    	rElements.push_back(std::pair<int, const ElementBase*>(ElementIter->first,ElementIter->second));
        ElementIter++;
    }
}

// store all elements of a structure in a vector
void NuTo::Structure::GetElementsTotal(std::vector<ElementBase*>& rElements)
{
    rElements.reserve(mElementMap.size());
    rElements.resize(0);
    boost::ptr_map<int,ElementBase>::iterator ElementIter = this->mElementMap.begin();
    while (ElementIter != this->mElementMap.end())
    {
        rElements.push_back(ElementIter->second);
        ElementIter++;
    }
}

// store all elements of a structure in a vector
void NuTo::Structure::GetElementsTotal(std::vector<std::pair<int, ElementBase*> >& rElements)
{
	rElements.reserve(mElementMap.size());
	rElements.resize(0);
	boost::ptr_map<int,ElementBase>::iterator ElementIter = this->mElementMap.begin();
    while (ElementIter != this->mElementMap.end())
    {
    	rElements.push_back(std::pair<int, ElementBase*>(ElementIter->first,ElementIter->second));
        ElementIter++;
    }
}

//! @param rGroupNumberElements group for elements (Plane2D4N) to be converted
//! @param rNodeDistanceMerge Distance of nodes to be joined (should be significantly smaller than the node distance in the mesh)
//! @param number of boxes to increase the search speed for the neighbor search
void NuTo::Structure::ElementConvertPlane2D4NToPlane2D4NSpectral (int rGroupNumberElements,
		int rOrder, double rNodeDistanceMerge, double rMeshSize)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
	double nodeDistanceMerge2 = rNodeDistanceMerge*rNodeDistanceMerge;
    boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(rGroupNumberElements);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::Structure::ElementConvertPlane2D4NToPlane2D4NSpectral] Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::Groups::Elements)
    	throw MechanicsException("[NuTo::Structure::ElementConvertPlane2D4NToPlane2D4NSpectral] Group is not an element group.");
    Group<ElementBase> *elementGroup = itGroup->second->AsGroupElement();
    assert(elementGroup!=0);

	//calculate bounding box for all existing nodes of the elements
	double bb[2][2];
	double coordinates[2];
	if (elementGroup->GetNumMembers()==0)
	{
		throw MechanicsException("[NuTo::Structure::ElementConvertPlane2D4NToPlane2D4NSpectral] number of elements in group to convert is zero");
	}
	if (elementGroup->begin()->second->GetNumNodesGeometry()==0)
	{
		throw MechanicsException("[NuTo::Structure::ElementConvertPlane2D4NToPlane2D4NSpectral] number of nodes for geometry interpolation in first element to convert is zero");
	}

	NodeBase* theNode = elementGroup->begin()->second->GetNodeGeometry(0);
	if (theNode->GetNumCoordinates()!=2)
		throw MechanicsException("[NuTo::Structure::ElementConvertPlane2D4NToPlane2D4NSpectral] coordinate dimension of node is not 2.)");
	theNode->GetCoordinates2D(coordinates);
	bb[0][0] = coordinates[0];
	bb[0][1] = coordinates[1];
	bb[1][0] = coordinates[0];
	bb[1][1] = coordinates[1];
	for (Group<ElementBase>::iterator itElement=elementGroup->begin(); itElement!=elementGroup->end();itElement++)
	{
		for (int count=0; count<itElement->second->GetNumNodes(); count++)
		{
			theNode = itElement->second->GetNode(count);
			if (theNode->GetNumCoordinates()!=2)
				throw MechanicsException("[NuTo::Structure::ElementConvertPlane2D4NToPlane2D4NSpectral] coordinate dimension of node is not 2.)");
			theNode->GetCoordinates2D(coordinates);
			if (coordinates[0]<bb[0][0])
				bb[0][0] = coordinates[0];
			if (coordinates[0]>bb[0][1])
				bb[0][1] = coordinates[0];
			if (coordinates[1]<bb[1][0])
				bb[1][0] = coordinates[1];
			if (coordinates[1]>bb[1][1])
				bb[1][1] = coordinates[1];
		}
	}

	//create boxes for the additional nodes
	double deltaBox[2];
	int numBoxes1D[2];
	numBoxes1D[0] = (bb[0][1]-bb[0][0])/rMeshSize+1;
	numBoxes1D[1] = (bb[1][1]-bb[1][0])/rMeshSize+1;
	deltaBox[0] = (bb[0][1]-bb[0][0])/numBoxes1D[0];
	deltaBox[1] = (bb[1][1]-bb[1][0])/numBoxes1D[1];

	//std::cout << "bounding box x " << bb[0][0] << " " << bb[1][0] << std::endl;
	//std::cout << "bounding box y " << bb[0][1] << " " << bb[1][1] << std::endl;
	//std::cout << "number of boxes " << numBoxes1D[0] << " " << numBoxes1D[1] << std::endl;
	std::vector<std::vector<NodeBase*> >nodeBoxes(numBoxes1D[0]*numBoxes1D[1]);

	if (rNodeDistanceMerge>deltaBox[0] || rNodeDistanceMerge>deltaBox[1])
	{
		throw MechanicsException("[NuTo::Structure::ElementConvertPlane2D4NToPlane2D4NSpectral] The merge distance is larger than the mesh size, that should not happen.");
	}

	Element::eElementType elementType;
	NuTo::IntegrationTypeBase* integrationType1D;
	switch (rOrder)
	{
	case 2:
		integrationType1D = this->GetPtrIntegrationType(NuTo::IntegrationType::IntegrationType1D2NLobatto3Ip);
		elementType = NuTo::Element::PLANE2D4NSPECTRALORDER2;
		break;
	case 3:
		integrationType1D = this->GetPtrIntegrationType(NuTo::IntegrationType::IntegrationType1D2NLobatto4Ip);
		elementType = NuTo::Element::PLANE2D4NSPECTRALORDER3;
		break;
	case 4:
		integrationType1D = this->GetPtrIntegrationType(NuTo::IntegrationType::IntegrationType1D2NLobatto5Ip);
		elementType = NuTo::Element::PLANE2D4NSPECTRALORDER4;
		break;
	default:
		throw MechanicsException("[NuTo::Structure::ElementConvertPlane2D4NToPlane2D4NSpectral] only order 2,3 and 4 implemented.");
	}

	//Get natural coordinates (1D)
	std::vector<double> newNodeCoordinatesNatural1D(integrationType1D->GetNumIntegrationPoints());

	//Get local coordinates (1D, 2D is just the tensor product) of the integration points/new node coordinates
	for (int theIP=0; theIP<integrationType1D->GetNumIntegrationPoints(); theIP++)
	{
		integrationType1D->GetLocalIntegrationPointCoordinates1D(theIP, newNodeCoordinatesNatural1D[theIP]);
	}

	//copy the "old" element pointer, since they are removed in the following loop
	// without the copy, the iterator over the group is no longer valid
	std::vector<std::pair<int,ElementBase*> > oldElements(elementGroup->GetNumMembers());
	unsigned int countElement(0);
	for (Group<ElementBase>::iterator itElement=elementGroup->begin(); itElement!=elementGroup->end();itElement++)
	{
		oldElements[countElement].first = itElement->first;
		oldElements[countElement].second = itElement->second;
		countElement++;
	}

	//add the additional nodes and create the elements
	for (countElement=0; countElement<oldElements.size(); countElement++)
	{
		try
		{
	        ElementBase* theElement = oldElements[countElement].second;
	        //std::cout << "start converting element " << oldElements[countElement].first << std::endl;

			std::vector<NuTo::NodeBase* > newNodeVector(integrationType1D->GetNumIntegrationPoints()*integrationType1D->GetNumIntegrationPoints());
            int curNode(0);
			for (int countY=0; countY<integrationType1D->GetNumIntegrationPoints(); countY++)
			{
				for (int countX=0; countX<integrationType1D->GetNumIntegrationPoints(); countX++)
				{
					if ((countX%(integrationType1D->GetNumIntegrationPoints()-1)==0) && (countY%(integrationType1D->GetNumIntegrationPoints()-1)==0))
					{
						//corner nodes already exist
						if (countX==0)
						{
							if (countY==0)
							{
								newNodeVector[curNode] = theElement->GetNode(0);
							}
							else
							{
								newNodeVector[curNode] = theElement->GetNode(3);

							}
						}
						else
						{
							if (countY==0)
							{
								newNodeVector[curNode] = theElement->GetNode(1);
							}
							else
							{
								newNodeVector[curNode] = theElement->GetNode(2);

							}
						}
						//std::cout << "use new node (local) " << curNode << " with global id " <<  NodeGetId(newNodeVector[curNode]) << std::endl;
					}
					else
					{
						double newNodeCoordinatesNatural[2];
						newNodeCoordinatesNatural[0] = newNodeCoordinatesNatural1D[countX];
						newNodeCoordinatesNatural[1] = newNodeCoordinatesNatural1D[countY];

						//calculate global coordinates
						double newNodeCoordinatesGlobal[3]; //this is a global routine for 3D, though z should always be zero
						theElement->InterpolateCoordinatesFrom2D(newNodeCoordinatesNatural,newNodeCoordinatesGlobal);

						//check if node already exists
						//calculate corresponding box
						int theBoxX = ((int)((newNodeCoordinatesGlobal[0]-bb[0][0])/deltaBox[0]))%(numBoxes1D[0]-1);
						int theBoxY = ((int)((newNodeCoordinatesGlobal[1]-bb[1][0])/deltaBox[1]))%(numBoxes1D[1]-1);
						int theBox = theBoxX + theBoxY*numBoxes1D[0];
						//std::cout << "theBox "<< theBox << " x "<< theBoxX << " y " << theBoxY << std::endl;

						//check distance to all nodes in this box and the neighboring boxes
						//(assuming that the merge distance is less than the box dimensions
						newNodeVector[curNode] = 0;
						for (int countX=-1; countX<2 && newNodeVector[curNode] == 0; countX++)
						{
							int curBoxX = theBoxX+countX;
							if (curBoxX<0 || curBoxX==numBoxes1D[0])
								continue;

							for (int countY=-1; countY<2 && newNodeVector[curNode] == 0; countY++)
							{
								int curBoxY = theBoxY+countY;
								if (curBoxY<0 || curBoxY==numBoxes1D[1])
									continue;

								int curBox = curBoxX+ curBoxY*numBoxes1D[0];
								//std::cout << "curBox "<< curBox << " x "<< curBoxX << " y " << curBoxY << std::endl;
								//loop over all nodes in the box and get the distance of those nodes to the new node position
								for (unsigned int countNode=0; countNode<nodeBoxes[curBox].size() && newNodeVector[curNode] == 0;countNode++)
								{
									double coordinates[2];
									nodeBoxes[curBox][countNode]->GetCoordinates2D(coordinates);
									double deltaX = coordinates[0]-newNodeCoordinatesGlobal[0];
									double deltaY = coordinates[1]-newNodeCoordinatesGlobal[1];
									double dist2 = deltaX*deltaX+deltaY*deltaY;
									if (dist2<nodeDistanceMerge2)
									{
										//node already exists, use the existing node
										assert(newNodeVector[curNode]==0);
										//std::cout << "merge new node (local) " << curNode << " with global node " <<  NodeGetId(nodeBoxes[curBox][countNode]) << std::endl;
										newNodeVector[curNode] = nodeBoxes[curBox][countNode];
									}
								}
							}
						}

						//if node does not exist, create new node by copying the first node of the element
						// using the same dofs etc.
						if (newNodeVector[curNode] == 0)
						{
							newNodeVector[curNode] = theElement->GetNode(0)->Clone();

							//set coordinates
							newNodeVector[curNode]->SetCoordinates2D(newNodeCoordinatesGlobal);

							//find unused integer id
						    int id(mNodeMap.size());
						    boost::ptr_map<int,NodeBase>::iterator it = mNodeMap.find(id);
						    while (it!=mNodeMap.end())
						    {
						        id++;
						        it = mNodeMap.find(id);
						    }

						    //insert into map
						    this->mNodeMap.insert(id, newNodeVector[curNode]);

						    //std::cout << "create new node (local) " << curNode << " with global id " <<  NodeGetId(newNodeVector[curNode]);
						    //std::cout << " with coordinates " << newNodeCoordinatesGlobal[0] << " " << newNodeCoordinatesGlobal[1] << std::endl;


							nodeBoxes[theBox].push_back(newNodeVector[curNode]);
						}
					}
					curNode++;
				}
			}
			//delete old element, but store the information related to id, element data and ipdata first
			int id = oldElements[countElement].first;
			ElementData::eElementDataType elementDataType(theElement->GetElementDataType());
			IpData::eIpDataType ipDataType(theElement->GetIpDataType(0));

			this->ElementDeleteInternal(id);

			//create new element with same id
			//std::cout << "create element " << id << std::endl;
			this->ElementCreate (id, elementType,newNodeVector , elementDataType, ipDataType);
		}
		catch(NuTo::MechanicsException &e)
		{
			e.AddMessage("[NuTo::Structure::ElementConvertPlane2D4NToPlane2D4NSpectral] Error converting element.");
			throw e;
		}
		catch(...)
		{
			throw NuTo::MechanicsException
			   ("[NuTo::Structure::ElementConvertPlane2D4NToPlane2D4NSpectral] error converting element.");
		}
	}

#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::ElementGroupDelete] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}


//! @param rGroupNumberElements group for elements (Plane2D4N) to be converted
//! @param rElementTypeStr new element type (must be a triangular type)
//! @param rNodeDistanceMerge distance between nodes for merging
//! @param rMeshSize for creation of a more efficient neighbor search
void NuTo::Structure::ElementConvertPlane2D3N (int rGroupNumberElements,
		std::string rElementTypeStr, double rNodeDistanceMerge, double rMeshSize)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
	double nodeDistanceMerge2 = rNodeDistanceMerge*rNodeDistanceMerge;
    boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(rGroupNumberElements);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::Structure::ElementConvertPlane2D3N] Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::Groups::Elements)
    	throw MechanicsException("[NuTo::Structure::ElementConvertPlane2D3N] Group is not an element group.");
    Group<ElementBase> *elementGroup = itGroup->second->AsGroupElement();
    assert(elementGroup!=0);

	// get element type
	std::string upperCaseElementType;
	std::transform(rElementTypeStr.begin(), rElementTypeStr.end(), std::back_inserter(upperCaseElementType), (int(*)(int)) toupper);

	Element::eElementType elementType;
	std::vector< std::array<double,2> > naturalNodeCoordinates;
	if (upperCaseElementType=="PLANE2D10N")
	{
		elementType = NuTo::Element::PLANE2D10N;
		NuTo::Plane2D10N::CalculateNaturalNodeCoordinates(naturalNodeCoordinates);
	}
	else if (upperCaseElementType=="PLANE2D6N")
	{
		elementType = NuTo::Element::PLANE2D6N;
		NuTo::Plane2D6N::CalculateNaturalNodeCoordinates(naturalNodeCoordinates);
	}
	else
	{
		throw MechanicsException("[NuTo::Structure::ElementConvertPlane2D3N] only conversion from plane3N to 6N and 10N implemented.");
	}

	//calculate bounding box for all existing nodes of the elements
	double bb[2][2];
	double coordinates[2];
	if (elementGroup->GetNumMembers()==0)
	{
		throw MechanicsException("[NuTo::Structure::ElementConvertPlane2D3N] number of elements in group to convert is zero");
	}
	if (elementGroup->begin()->second->GetNumNodesGeometry()==0)
	{
		throw MechanicsException("[NuTo::Structure::ElementConvertPlane2D3N] number of nodes for geometry interpolation in first element to convert is zero");
	}

	NodeBase* theNode = elementGroup->begin()->second->GetNodeGeometry(0);
	if (theNode->GetNumCoordinates()!=2)
		throw MechanicsException("[NuTo::Structure::ElementConvertPlane2D3N] coordinate dimension of node is not 2.)");
	theNode->GetCoordinates2D(coordinates);
	bb[0][0] = coordinates[0];
	bb[0][1] = coordinates[1];
	bb[1][0] = coordinates[0];
	bb[1][1] = coordinates[1];
	for (Group<ElementBase>::iterator itElement=elementGroup->begin(); itElement!=elementGroup->end();itElement++)
	{
		for (int count=0; count<itElement->second->GetNumNodes(); count++)
		{
			theNode = itElement->second->GetNode(count);
			if (theNode->GetNumCoordinates()!=2)
				throw MechanicsException("[NuTo::Structure::ElementConvertPlane2D3N] coordinate dimension of node is not 2.)");
			theNode->GetCoordinates2D(coordinates);
			if (coordinates[0]<bb[0][0])
				bb[0][0] = coordinates[0];
			if (coordinates[0]>bb[0][1])
				bb[0][1] = coordinates[0];
			if (coordinates[1]<bb[1][0])
				bb[1][0] = coordinates[1];
			if (coordinates[1]>bb[1][1])
				bb[1][1] = coordinates[1];
		}
	}

	//create boxes for the additional nodes
	double deltaBox[2];
	int numBoxes1D[2];
	numBoxes1D[0] = (bb[0][1]-bb[0][0])/rMeshSize+1;
	numBoxes1D[1] = (bb[1][1]-bb[1][0])/rMeshSize+1;
	deltaBox[0] = (bb[0][1]-bb[0][0])/numBoxes1D[0];
	deltaBox[1] = (bb[1][1]-bb[1][0])/numBoxes1D[1];

	//std::cout << "bounding box x " << bb[0][0] << " " << bb[1][0] << std::endl;
	//std::cout << "bounding box y " << bb[0][1] << " " << bb[1][1] << std::endl;
	//std::cout << "number of boxes " << numBoxes1D[0] << " " << numBoxes1D[1] << std::endl;
	std::vector<std::vector<NodeBase*> >nodeBoxes(numBoxes1D[0]*numBoxes1D[1]);

	if (rNodeDistanceMerge>deltaBox[0] || rNodeDistanceMerge>deltaBox[1])
	{
		throw MechanicsException("[NuTo::Structure::ElementConvertPlane2D3N] The merge distance is larger than the mesh size, that should not happen.");
	}

	//copy the "old" element pointer, since they are removed in the following loop
	// without the copy, the iterator over the group is no longer valid
	std::vector<std::pair<int,ElementBase*> > oldElements(elementGroup->GetNumMembers());
	unsigned int countElement(0);
	for (Group<ElementBase>::iterator itElement=elementGroup->begin(); itElement!=elementGroup->end();itElement++)
	{
		oldElements[countElement].first = itElement->first;
		oldElements[countElement].second = itElement->second;
		countElement++;
	}

	//add the additional nodes and create the elements
	for (countElement=0; countElement<oldElements.size(); countElement++)
	{
		try
		{
	        ElementBase* theElement = oldElements[countElement].second;
	        //std::cout << "start converting element " << oldElements[countElement].first << std::endl;

			std::vector<NuTo::NodeBase* > newNodeVector(naturalNodeCoordinates.size());
			for (int curNode=0; curNode<naturalNodeCoordinates.size(); curNode++)
			{
				newNodeVector[curNode] = 0;

				//corner nodes already exist
				std::cout << "naturalCoordinates " << naturalNodeCoordinates[curNode][0] << " " << naturalNodeCoordinates[curNode][1] << std::endl;
				if (fabs(naturalNodeCoordinates[curNode][0]<1e-8))
				{
					if (fabs(naturalNodeCoordinates[curNode][1])<1e-8)
					{
						newNodeVector[curNode] = theElement->GetNode(0);
					}
					if (fabs(naturalNodeCoordinates[curNode][1]-1.)<1e-8)
					{
						newNodeVector[curNode] = theElement->GetNode(2);
					}
				}
				else
				{
					if (fabs(naturalNodeCoordinates[curNode][0]-1.)<1e-8 &&
					    fabs(naturalNodeCoordinates[curNode][1])<1e-8)
					{
						newNodeVector[curNode] = theElement->GetNode(1);
					}
				}

				if (newNodeVector[curNode] != 0)
				{
					double newNodeCoordinatesGlobal[2]; //this is a global routine for 3D, though z should always be zero
					newNodeVector[curNode]->GetCoordinates2D(newNodeCoordinatesGlobal);
					std::cout << "coordinates of old/new " << curNode << ". " << newNodeCoordinatesGlobal[0] << " " << newNodeCoordinatesGlobal[1] << std::endl;
				}
				else
				{
					//calculate global coordinates
					double newNodeCoordinatesGlobal[3]; //this is a global routine for 3D, though z should always be zero
					theElement->InterpolateCoordinatesFrom2D(&(naturalNodeCoordinates[curNode][0]),newNodeCoordinatesGlobal);

					//check if node already exists
					//calculate corresponding box
					int theBoxX = ((int)((newNodeCoordinatesGlobal[0]-bb[0][0])/deltaBox[0]))%(numBoxes1D[0]-1);
					int theBoxY = ((int)((newNodeCoordinatesGlobal[1]-bb[1][0])/deltaBox[1]))%(numBoxes1D[1]-1);
					int theBox = theBoxX + theBoxY*numBoxes1D[0];
					//std::cout << "theBox "<< theBox << " x "<< theBoxX << " y " << theBoxY << std::endl;

					//check distance to all nodes in this box and the neighboring boxes
					//(assuming that the merge distance is less than the box dimensions
					newNodeVector[curNode] = 0;
					for (int countX=-1; countX<2 && newNodeVector[curNode] == 0; countX++)
					{
						int curBoxX = theBoxX+countX;
						if (curBoxX<0 || curBoxX==numBoxes1D[0])
							continue;

						for (int countY=-1; countY<2 && newNodeVector[curNode] == 0; countY++)
						{
							int curBoxY = theBoxY+countY;
							if (curBoxY<0 || curBoxY==numBoxes1D[1])
								continue;

							int curBox = curBoxX+ curBoxY*numBoxes1D[0];
							//std::cout << "curBox "<< curBox << " x "<< curBoxX << " y " << curBoxY << std::endl;
							//loop over all nodes in the box and get the distance of those nodes to the new node position
							for (unsigned int countNode=0; countNode<nodeBoxes[curBox].size() && newNodeVector[curNode] == 0;countNode++)
							{
								double coordinates[2];
								nodeBoxes[curBox][countNode]->GetCoordinates2D(coordinates);
								double deltaX = coordinates[0]-newNodeCoordinatesGlobal[0];
								double deltaY = coordinates[1]-newNodeCoordinatesGlobal[1];
								double dist2 = deltaX*deltaX+deltaY*deltaY;
								if (dist2<nodeDistanceMerge2)
								{
									//node already exists, use the existing node
									assert(newNodeVector[curNode]==0);
									//std::cout << "merge new node (local) " << curNode << " with global node " <<  NodeGetId(nodeBoxes[curBox][countNode]) << std::endl;
									newNodeVector[curNode] = nodeBoxes[curBox][countNode];
								}
							}
						}
					}

					//if node does not exist, create new node by copying the first node of the element
					// using the same dofs etc.
					if (newNodeVector[curNode] == 0)
					{
						newNodeVector[curNode] = theElement->GetNode(0)->Clone();

						//set coordinates
						newNodeVector[curNode]->SetCoordinates2D(newNodeCoordinatesGlobal);
						std::cout << "coordinates of new " << curNode << "  " << newNodeCoordinatesGlobal[0] << " " << newNodeCoordinatesGlobal[1] << std::endl;

						//find unused integer id
						int id(mNodeMap.size());
						boost::ptr_map<int,NodeBase>::iterator it = mNodeMap.find(id);
						while (it!=mNodeMap.end())
						{
							id++;
							it = mNodeMap.find(id);
						}

						//insert into map
						this->mNodeMap.insert(id, newNodeVector[curNode]);

						//std::cout << "create new node (local) " << curNode << " with global id " <<  NodeGetId(newNodeVector[curNode]);
						//std::cout << " with coordinates " << newNodeCoordinatesGlobal[0] << " " << newNodeCoordinatesGlobal[1] << std::endl;


						nodeBoxes[theBox].push_back(newNodeVector[curNode]);
					}
				}
			}
			//delete old element, but store the information related to id, element data and ipdata first
			int id = oldElements[countElement].first;
			ElementData::eElementDataType elementDataType(theElement->GetElementDataType());
			IpData::eIpDataType ipDataType(theElement->GetIpDataType(0));

			this->ElementDeleteInternal(id);

			//create new element with same id
			//std::cout << "create element " << id << std::endl;
			this->ElementCreate (id, elementType,newNodeVector , elementDataType, ipDataType);
		}
		catch(NuTo::MechanicsException &e)
		{
			e.AddMessage("[NuTo::Structure::ElementConvertPlane2D3N] Error converting element.");
			throw e;
		}
		catch(...)
		{
			throw NuTo::MechanicsException
			   ("[NuTo::Structure::ElementConvertPlane2D3N] error converting element.");
		}
	}

#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::ElementConvertPlane2D3N] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}
