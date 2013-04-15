// $Id$

#include <sstream>

#include <boost/foreach.hpp>

#include "nuto/base/Debug.h"

#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/cracks/CrackExplicit2D.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/MechanicsException.h"

#include "nuto/mechanics/elements/ElementDataConstitutiveIpCrack.h"
#include "nuto/mechanics/integrationtypes/IntegrationPointBase.h"

//! @brief returns the number of cracks
//! @return number of cracks
unsigned int NuTo::Structure::GetNumCracks() const
{
    return mCrackMap.size();
}

//! @brief returns a reference to a crack
//! @param identifier
//! @return reference to a crack
NuTo::CrackBase* NuTo::Structure::CrackGetCrackPtr(int rIdent)
{
    crackMap_t::iterator it = mCrackMap.find(rIdent);
    if (it!=mCrackMap.end())
        return it->second;
    else
    {
    	std::stringstream message;
    	message << "[NuTo::Structure::CrackGetCrackPtr] Crack with identifier " << rIdent << " does not exist." << std::endl;
    	throw MechanicsException(message.str());
    }
}

//! @brief returns a const reference to a crack
//! @param identifier
//! @return const reference to a crack
const NuTo::CrackBase* NuTo::Structure::CrackGetCrackPtr(int rIdent)const
{
    crackMap_t::const_iterator it = mCrackMap.find(rIdent);
    if (it!=mCrackMap.end())
        return it->second;
    else
    {
    	std::stringstream message;
    	message << "[NuTo::Structure::CrackGetCrackPtr] Crack with identifier " << rIdent << " does not exist." << std::endl;
    	throw MechanicsException(message.str());
    }
}

//! @brief gives the identifier of a crack
//! @param pointer to a crack
//! @return identifier
int NuTo::Structure::CrackGetId(const CrackBase* rCrack)const
{
    for (crackMap_t::const_iterator
            it = mCrackMap.begin(); it!= mCrackMap.end(); it++)
    {
        if (it->second==rCrack)
            return it->first;
    }
    throw MechanicsException("[NuTo::Structure::CrackGetId] Crack does not exist.");
}

//! @brief ... Info routine that prints general information about the cracks
//! @param ... rVerboseLevel describes how detailed the information is
void NuTo::Structure::CrackInfo(int rVerboseLevel)const
{
    std::cout << "number of cracks  : " << mCrackMap.size() << std::endl;
    if (rVerboseLevel>2)
    {
        for (crackMap_t::const_iterator it= mCrackMap.begin(); it!=mCrackMap.end(); it++)
        {
            std::cout << "  Crack " << it->first << std::endl;
            it->second->Info(rVerboseLevel);
            std::cout <<  std::endl;
        }
    }

}

//! @brief ... delete an existing crack
//! @param rIdent ... crack identifier
void NuTo::Structure::CrackDelete(int rIdent)
{
    // find section identifier in map
    crackMap_t::iterator it = this->mCrackMap.find(rIdent);
    if (it == this->mCrackMap.end())
    {
        throw NuTo::MechanicsException("[" + std::string(__PRETTY_FUNCTION__) + "] Crack does not exist.");
    }
    else
    {
        this->mCrackMap.erase(it);
    }
}

//! @brief ... create a new crack
//! @param rIdent ... crack identifier
int NuTo::Structure::CrackCreate()
{
	NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> nodes(0,0);

	return NuTo::Structure::CrackCreate(nodes);
}

//! @brief ... create a new crack with given node-Id's
//! @param rNodes ... vector of node-Id's
int NuTo::Structure::CrackCreate(NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic>& rNodes)
{

	//find unused integer id
    int id(0);
    crackMap_t::iterator it = mCrackMap.find(id);
    while (it!=mCrackMap.end())
    {
        id++;
        it = mCrackMap.find(id);
    }

	// add crack to map

    //! switch dimension of structure
    switch (this->GetDimension())
    {
    case 1:
    	throw NuTo::MechanicsException("[" + std::string(__PRETTY_FUNCTION__) + "] Cracks are not implemented for 1D yet!");
    case 2:
    	this->mCrackMap.insert(id, new NuTo::CrackExplicit2D );
    	break;
    case 3:
    	throw NuTo::MechanicsException("[" + std::string(__PRETTY_FUNCTION__) + "] Cracks are not implemented for 3D yet!");
    default:
    	throw NuTo::MechanicsException("[" + std::string(__PRETTY_FUNCTION__) + "] Wrong dimension of structure found!!!");
    }

    size_t numNodes=rNodes.GetNumRows();
    for(size_t thisNodeCount=0;thisNodeCount<numNodes;++thisNodeCount)
    	CrackPushBack(id,this->NodeGetNodePtr(rNodes(thisNodeCount,0)));

	return id;

}

//! @brief ... extends an existing crack
//! @param rIdent (Input) ... crack identifier
//! @param rNode (Input) ... node Id to be attended to the crack
void NuTo::Structure::CrackPushBack(const int rIdent, const int rNodeNumber)
{
	nodeBasePtr_t nodePtr(this->NodeGetNodePtr(rNodeNumber));
	CrackPushBack(rIdent,nodePtr);
}


//! @brief ... extends an existing crack
//! @param rIdent (Input) ... crack identifier
//! @param rNode (Input) ... node Id to be attended to the crack
void NuTo::Structure::CrackPushFront(const int rIdent, const int rNodeNumber)
{
	nodeBasePtr_t nodePtr(this->NodeGetNodePtr(rNodeNumber));
	CrackPushFront(rIdent,nodePtr);
}

//! @brief ... extends an existing crack
//! @param rIdent ... crack identifier
//! @param rNode ... pointer to the node to be attended to the crack
void NuTo::Structure::CrackPushBack(int rIdent, NuTo::NodeBase* rNode)
{
    switch (this->GetDimension())
    {
    case 2:
		{
			/// find section identifier in map
			boost::ptr_map<int,CrackExplicit2D>::iterator it = this->mCrackMap.find(rIdent);
			if (it == this->mCrackMap.end())
				throw NuTo::MechanicsException("[" + std::string(__PRETTY_FUNCTION__) + "] Crack does not exist.");
			else
				(it->second)->PushBack(rNode);
		}
    	break;
    default:
    	throw NuTo::MechanicsException("[" + std::string(__PRETTY_FUNCTION__) + "] Not implemented for this dimension of structure!!!");
    }
}

//! @brief ... extends an existing crack
//! @param rIdent ... crack identifier
//! @param rNode ... pointer to the node to be attended to the crack
void NuTo::Structure::CrackPushFront(int rIdent, NuTo::NodeBase* rNode)
{
    switch (this->GetDimension())
    {
    case 2:
		{
			/// find section identifier in map
			boost::ptr_map<int,CrackExplicit2D>::iterator it = this->mCrackMap.find(rIdent);
			if (it == this->mCrackMap.end())
				throw NuTo::MechanicsException("[" + std::string(__PRETTY_FUNCTION__) + "] Crack does not exist.");
			else
				(it->second)->PushFront(rNode);
		}
    	break;
    default:
    	throw NuTo::MechanicsException("[" + std::string(__PRETTY_FUNCTION__) + "] Not implemented for this dimension of structure!!!");
    }
}


//! @brief ... shortens an existing crack
//! @param rIdent ... crack identifier
void NuTo::Structure::CrackPopBack(int rIdent)
{
    switch (this->GetDimension())
    {
    case 2:
		{
			/// find section identifier in map
			boost::ptr_map<int,CrackExplicit2D>::iterator it = this->mCrackMap.find(rIdent);
			if (it == this->mCrackMap.end())
				throw NuTo::MechanicsException("[" + std::string(__PRETTY_FUNCTION__) + "] Crack does not exist.");
			else
				(it->second)->PopBack();
		}
    	break;
    default:
    	throw NuTo::MechanicsException("[" + std::string(__PRETTY_FUNCTION__) + "] Not implemented for this dimension of structure!!!");
    }
}

//! @brief ... shortens an existing crack
//! @param rIdent ... crack identifier
void NuTo::Structure::CrackPopFront(int rIdent)
{
    switch (this->GetDimension())
    {
    case 2:
		{
			/// find crack identifier in map
			boost::ptr_map<int,CrackExplicit2D>::iterator it = this->mCrackMap.find(rIdent);
			if (it == this->mCrackMap.end())
				throw NuTo::MechanicsException("[" + std::string(__PRETTY_FUNCTION__) + "] Crack does not exist.");
			else
				(it->second)->PopFront();
		}
    	break;
    default:
    	throw NuTo::MechanicsException("[" + std::string(__PRETTY_FUNCTION__) + "] Not implemented for this dimension of structure!!!");
    }
}


//! @brief ... store all cracks of a structure in a vector
//! @param rElements ... vector of const crack pointer
void NuTo::Structure::GetCracksTotal(std::vector<const CrackBase*>& rCracks) const
{
	rCracks.reserve(mCrackMap.size());
	crackMap_t::const_iterator CrackIter = this->mCrackMap.begin();
    while (CrackIter != this->mCrackMap.end())
    {
    	rCracks.push_back(CrackIter->second);
    	CrackIter++;
    }
}

//! @brief ... store all cracks of a structure in a vector
//! @param rElements ... vector of crack pointer
void NuTo::Structure::GetCracksTotal(std::vector<CrackBase*>& rCracks)
{
	rCracks.reserve(mCrackMap.size());
	crackMap_t::iterator CrackIter = this->mCrackMap.begin();
    while (CrackIter != this->mCrackMap.end())
    {
    	rCracks.push_back(CrackIter->second);
    	CrackIter++;
    }
}

//! @brief ... merge all cracks to the existing structure
void NuTo::Structure::InitiateCracks()
{
	elementBasePtrSet_t crackedElems;
	NuTo::Structure::InitiateCracks(crackedElems);
}

//! @brief ... merge all cracks to the existing structure
//! @param rCrackedElems (Output) ... vector of cracked elements
void NuTo::Structure::InitiateCracks(elementBasePtrSet_t & rCrackedElems)
{
	// go through all cracks
	BOOST_FOREACH( crackMap_t::value_type thisCrack, mCrackMap )
		InitiateCrack(thisCrack.first, rCrackedElems);
}

//! @brief ... merge specified crack to the existing structure
//! @param rIdent ... crack identifier
//! @param rCrackedElems (Output) ... vector of cracked elements
void NuTo::Structure::InitiateCrack(const int rIdent)
{
	elementBasePtrSet_t crackedElems;
	NuTo::Structure::InitiateCrack(rIdent, crackedElems);
}

//! @brief ... merge specified crack to the existing structure
//! @param rIdent ... crack identifier
//! @param rCrackedElems (Output) ... vector of cracked elements
void NuTo::Structure::InitiateCrack(const int rIdent, elementBasePtrSet_t & rCrackedElems)
{
	//! find crack identifier in map
	NuTo::CrackBase* thisCrack = this->CrackGetCrackPtr(rIdent);

	//! initiate crack into given elements
	elementBasePtrVec_t allElems;
	this->GetElementsTotal ( allElems );

	//! now the crack searching the cracked elements
	elementBasePtrVec_t rThisCrackedElems;
	thisCrack->Initiate(allElems, rThisCrackedElems);

	BOOST_FOREACH( elementBasePtr_t thisElPtr, rThisCrackedElems )
	{
		if(thisElPtr->GetElementDataType()!=NuTo::ElementData::CONSTITUTIVELAWIPCRACK)
	    	throw NuTo::MechanicsException("[NuTo::Structure::InitiateCrack] Not implemented for this type of ElementDataType!!!");

		ElementDataBase* ptrElementData=thisElPtr->GetDataPtr();
		/*unsigned int thisCrackId=*/ptrElementData->AddCrack(thisCrack);

		//! check if rCrackedElems contains rThisCrackedElems already, if not: append
		rCrackedElems.insert(thisElPtr);

		/*std::cout << "Element to be cracked: " << ElementGetId(thisElPtr)  << "; mElementData=" << ptrElementData
				<< "; mIntegrationType=" << ptrElementData->GetIntegrationType()
				<< "; mElementDataType=" << thisElPtr->GetElementDataType()
				<< "; thisCrackId=" << thisCrackId
		<< std::endl;
		*/
	}
}

//! @brief ... initiate PhantomNodeMethod
/**
 * This function takes all cracks and merge it to the structure.
 */
void NuTo::Structure::InitiatePhantomNodeMethod()
{
    elementBasePtrSet_t crackedElems;
	this->InitiateCracks(crackedElems);
	this->InitiatePhantomNodeMethod(crackedElems);
}

//! @brief ... take cracked elements and initiate PhantomNodeMethod
//! @param rCrackedElems (Input) ... vector of cracked elements
void NuTo::Structure::InitiatePhantomNodeMethod(elementBasePtrSet_t & rCrackedElems)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
	nodeBasePtrMap_t crackedNodes;

    //! crack the elements
    BOOST_FOREACH(NuTo::Structure::elementBasePtr_t thisCrackedElem, rCrackedElems)
    {
		//! at first only with one Crack!!
		NuTo::Structure::crackBasePtrVec_t crackPtrVec(thisCrackedElem->GetCracks());
		if(crackPtrVec.size()!=1)
			throw NuTo::MechanicsException("[NuTo::Structure::InitiatePhantomNodeMethod] Only 1 crack per element possible!!");
		NuTo::Structure::crackBasePtr_t thisCrackPtr=crackPtrVec[0];
		//! now find the cracked edge and the position of the crack
		//! assumption: linear elements (incremental numeration leads to the edges)
		const unsigned short numElemNodes=thisCrackedElem->GetNumNodes ();

		/** @brief the vector for the node pointers
		 * build up the nodepointer vectors
		 * initialize with the original node pointers
		 */
		nodeBasePtrVec_t oldNodes(numElemNodes), newNodes;
		for(unsigned short i=numElemNodes; i--;)
		{
			oldNodes[i]=thisCrackedElem->GetNode(i);
			newNodes=oldNodes;
		}
		unsigned short locNodeId1stEdge=0, locNodeId2ndEdge=0;
		for(locNodeId1stEdge=0; locNodeId1stEdge<numElemNodes; locNodeId1stEdge++)
		{
			// create a temporary point
			int tmpNodeId=this->NodeCreate("Coordinates");
			double relCoor=0.0;
			size_t segment;
			bool edgeIntersected=thisCrackPtr->Intersect(
					thisCrackedElem->GetNode(locNodeId1stEdge),
					thisCrackedElem->GetNode((locNodeId1stEdge+1)%numElemNodes),
					this->NodeGetNodePtr(tmpNodeId), relCoor, segment);
			if(edgeIntersected)
			{
				/**
				 * first build up new integration types
				 * 1) building up crack segment in element
				 * 2) create new nodes
				 * 3) create the two new elements
				 * 4) get number and position of integration points and create two new modifiable integration types
				 * 5) impose the integration scheme
				 * 6) delete the original element
				 */

				//! building up crack segment in element
				//! build up integration cells for the two elements
				double xA[2], xB[2];
				this->NodeGetNodePtr(tmpNodeId)->GetCoordinates2D(xA);
				//! check if the crack running thru the element
				unsigned short numIntersectedEdges=1; //!< one edge already intersected
				for(locNodeId2ndEdge=locNodeId1stEdge+1; locNodeId2ndEdge<numElemNodes; locNodeId2ndEdge++)
				{
					if(thisCrackPtr->Intersect(
							thisCrackedElem->GetNode(locNodeId2ndEdge),
							thisCrackedElem->GetNode(((locNodeId2ndEdge+1)%numElemNodes)),
							this->NodeGetNodePtr(tmpNodeId), relCoor, segment))
					{
						numIntersectedEdges++;
						break; //!< now we have two intersections --> it's enough
					}
				}
				//! check intersection type
				/**
				 * there are two cases for the cracked elements:
				 * 1) Elements which are torn in two pieces
				 * 2) Elements containing the cracktip
				 *
				 * Therefore we need 4 new nodes for case 1).
				 * For the 2nd case we need only two new nodes at the cracked edge.
				 */

				//! create the first two new nodes
				//! Node A
				nodeBasePtr_t oldNodePtr=thisCrackedElem->GetNode(locNodeId1stEdge);
				nodeBasePtr_t newNodePtr=NULL;
				nodeBasePtrMap_t::iterator crackedNodePairIt = crackedNodes.find(oldNodePtr);
				if(crackedNodePairIt==crackedNodes.end())
				{
					const int newNode=this->NodeCreate("displacements");
					newNodePtr=this->NodeGetNodePtr(newNode);
					double val2D[2];
					oldNodePtr->GetCoordinates2D(val2D);
					newNodePtr->SetCoordinates2D(val2D);
					oldNodePtr->GetDisplacements2D(val2D);
					newNodePtr->SetDisplacements2D(val2D);
					crackedNodes.insert(std::make_pair(oldNodePtr,newNodePtr));
				}else{
					newNodePtr=crackedNodePairIt->second;
				}
				newNodes[locNodeId1stEdge]=newNodePtr;
				//! Node B
				oldNodePtr=thisCrackedElem->GetNode(((locNodeId1stEdge+1)%numElemNodes));
				crackedNodePairIt = crackedNodes.find(oldNodePtr);
				if(crackedNodePairIt==crackedNodes.end())
				{
					const int newNode=this->NodeCreate("displacements");
					newNodePtr=this->NodeGetNodePtr(newNode);
					double val2D[2];
					oldNodePtr->GetCoordinates2D(val2D);
					newNodePtr->SetCoordinates2D(val2D);
					oldNodePtr->GetDisplacements2D(val2D);
					newNodePtr->SetDisplacements2D(val2D);
					crackedNodes.insert(std::make_pair(oldNodePtr,newNodePtr));
				}else{
					newNodePtr=crackedNodePairIt->second;
				}
				newNodes[(locNodeId1stEdge+1)%numElemNodes]=newNodePtr;

				if(numIntersectedEdges==2){ //!< if  numIntersectedEdges==2 -> crack running thru
					this->NodeGetNodePtr(tmpNodeId)->GetCoordinates2D(xB);
					//! for the torn element we have to copy all other element's nodes
					for(unsigned short i=numElemNodes; i--;)
					{
						if(i==locNodeId1stEdge || i==(locNodeId1stEdge+1)%numElemNodes) continue;
						oldNodePtr=thisCrackedElem->GetNode(i);
						crackedNodePairIt = crackedNodes.find(oldNodePtr);
						if(crackedNodePairIt==crackedNodes.end())
						{
							const int newNode=this->NodeCreate("displacements");
							newNodePtr=this->NodeGetNodePtr(newNode);
							double val2D[2];
							oldNodePtr->GetCoordinates2D(val2D);
							newNodePtr->SetCoordinates2D(val2D);
							oldNodePtr->GetDisplacements2D(val2D);
							newNodePtr->SetDisplacements2D(val2D);
							crackedNodes.insert(std::make_pair(oldNodePtr,newNodePtr));
						}else{
							newNodePtr=crackedNodePairIt->second;
						}
						newNodes[i]=newNodePtr;
					}
				}else if(numIntersectedEdges==1){ //!< if  numIntersectedEdges==1 -> crack ends inside the element
					bool edgeCracked=false;
					//! find the projection of the crack to the element corners.
					for(locNodeId2ndEdge=0; locNodeId2ndEdge<numElemNodes; locNodeId2ndEdge++)
					{
						if(locNodeId2ndEdge==locNodeId1stEdge) continue; //!< this edge is already cracked
						int newCrackEnd=this->NodeCreate("Coordinates");
						//! @todo check if the intersecting crack segment is at the end or beginning of the crack
						const unsigned short rayIntersect=thisCrackPtr->ExtendEnd(
								thisCrackedElem->GetNode(locNodeId2ndEdge),
								thisCrackedElem->GetNode(((locNodeId2ndEdge+1)%numElemNodes)),
								this->NodeGetNodePtr(newCrackEnd), relCoor);
						if(rayIntersect)
						{
							//! check the side of the crack: line segment AB is intersected at (0=no intersection, 1=front, 2=end)
							if(rayIntersect==2)
							{
								thisCrackPtr->PushBack(this->NodeGetNodePtr(newCrackEnd));
								//~ this->NodeGetNodePtr(newCrackEnd)->GetCoordinates2D(xB);
							}else if(rayIntersect==1)
							{
								//! in this case we add a new Node to the front of the crack
								//! and for consistency we change the direction of the crack xA->xB
								thisCrackPtr->PushFront(this->NodeGetNodePtr(newCrackEnd));
								//~ unsigned short i=0;
								//~ BOOST_FOREACH(double d, xA) xB[i++]=d;
								//~ this->NodeGetNodePtr(newCrackEnd)->GetCoordinates2D(xA);
							}else
								throw NuTo::MechanicsException("[NuTo::Structure::InitiatePhantomNodeMethod] Got a wrong end of the crack!!!");
							this->NodeGetNodePtr(newCrackEnd)->GetCoordinates2D(xB);
							edgeCracked=true;
							break; //!< now we have two intersections --> it's enough
						}else{
							this->NodeDelete(newCrackEnd,false);
						}
					}
					if(!edgeCracked)
						throw NuTo::MechanicsException("[NuTo::Structure::InitiatePhantomNodeMethod] Edge is not parted by end ray of the crack");
				}else{
					throw NuTo::MechanicsException("[NuTo::Structure::InitiatePhantomNodeMethod] Unhandeled case!!");
				}

				//! @brief 3) create the two new elements
				/**
				 * now we have two intersections of the element's corners.
				 * With this intersection nodes we build up the integration cells,
				 * therefore we check which integration points are outside this area.
				 */
				//! calculate the normal to the cracksegment
				//! @todo get this from a geometry class
				/*!
					In the two-dimensional case, the normal $\boldsymbol{n}$ of the crack is defined as
					\f[\boldsymbol{n} =
						 \dfrac{ \left( \boldsymbol{x}_B - \boldsymbol{x}_A \right) \times \boldsymbol{e}_z}
							   { \left| \left( \boldsymbol{x}_B - \boldsymbol{x}_A \right) \times \boldsymbol{e}_z \right|}
					\f]
				*/
				double normalVec[2] = { xB[1] - xA[1] , xA[0] - xB[0]};
				double norm = sqrt(normalVec[0]*normalVec[0]+normalVec[1]*normalVec[1]);
				normalVec[0] /= norm;
				normalVec[1] /= norm;

				//! @todo Remove this dirty hack --> make it general!!! (perhaps element copy constructor)
				std::vector<nodeBasePtr_t> rNodeVectorA(oldNodes);
				std::vector<nodeBasePtr_t> rNodeVectorB(newNodes);
				//! check the position of the nodes
				for(unsigned short i=oldNodes.size(); i--;)
				{
					double val2D[2];
					oldNodes[i]->GetCoordinates2D(val2D);
					//! calculate the distance of the node to the cracksegment
					/**
						In the two-dimensional case, with the normalized normal $\boldsymbol{n}$ of the crack the distance $d$ is defined as
						\f[ d = \left( \boldsymbol{x}_B - \boldsymbol{x}_A \right) \cdot \boldsymbol{n}
						\f]
					*/
					double dist= (val2D[0] - xA[0])*normalVec[0] +(val2D[1] - xA[1])*normalVec[1];
					/**
					 *  for this first element all IP on the left side are required
					 *  --> all IPs with positive distances have to be deleted in this element
					 *  	--> all IPs without signbit
					 *  --> all other IPs have to be deleted for the second element
					 *  first collect all not needed IPs (don't delete it now, otherwise you will get problems with the iterators inside the integration type)
					 */
					if(!std::signbit(dist))
					{
						rNodeVectorA[i]=newNodes[i];
						rNodeVectorB[i]=oldNodes[i];
					}
				}
				//! create the new element A
				const size_t numElA = this->ElementCreate(thisCrackedElem->GetEnumType(),rNodeVectorA,thisCrackedElem->GetElementDataType(), NuTo::IpData::NOIPDATA);
				elementBasePtr_t newElPtrA = this->ElementGetElementPtr(numElA);
				//! create the new element B
				const size_t numElB = this->ElementCreate(thisCrackedElem->GetEnumType(),rNodeVectorB,thisCrackedElem->GetElementDataType(), NuTo::IpData::NOIPDATA);
				elementBasePtr_t newElPtrB = this->ElementGetElementPtr(numElB);

				//! @brief 4) get numer and position of integration points and create two new modifiable integration types
				double ipCoor[3];
				//! get general information and informations about integration type for this element
				size_t numIp=thisCrackedElem->GetIntegrationType()->GetNumIntegrationPoints();
				std::stringstream  intTypeStr;
				switch(thisCrackedElem->GetEnumType())
				{
				case NuTo::Element::PLANE2D3N:
					intTypeStr << "2D3NGAUSS3IP";
					break;
				case NuTo::Element::PLANE2D4N:
					intTypeStr << "2D4NMOD" << numIp << "IP";
					break;
				default:
					throw NuTo::MechanicsException("[" + std::string(__PRETTY_FUNCTION__) + "] Wrong Elementtype!");
				}
				std::vector<size_t> outsideIpsA(0), outsideIpsB(0);
				//! get the global coordinates of each IP
				for(size_t ip=0; ip<numIp; ++ip){
					thisCrackedElem->GetGlobalIntegrationPointCoordinates(ip,ipCoor);
					//! calculate the distance of the node to the cracksegment
					/**
						In the two-dimensional case, with the normalized normal $\boldsymbol{n}$ of the crack the distance $d$ is defined as
						\f[ d = \left( \boldsymbol{x}_B - \boldsymbol{x}_A \right) \cdot \boldsymbol{n}
						\f]
					*/
					double dist= (ipCoor[0] - xA[0])*normalVec[0] +(ipCoor[1] - xA[1])*normalVec[1];
					/**
					 *  for this first element all IP on the left side are required
					 *  --> all IPs with positive distances have to be deleted in this element
					 *  	--> all IPs without signbit
					 *  --> all other IPs have to be deleted for the second element
					 *  first collect all not needed IPs (don't delete it now, otherwise you will get problems with the iterators inside the integration type)
					 */
					if(!std::signbit(dist))	outsideIpsA.push_back(ip);
						else 				outsideIpsB.push_back(ip);
				}
				//! build up new integration types
				std::stringstream  intTypeStrA, intTypeStrB;
				intTypeStrA << intTypeStr.str() << numElA;
				intTypeStrB << intTypeStr.str() << numElB;
				NuTo::IntegrationTypeBase* intTypePtrA(this->GetPtrIntegrationType(intTypeStrA.str()));
				NuTo::IntegrationTypeBase* intTypePtrB(this->GetPtrIntegrationType(intTypeStrB.str()));
				//! delete non-needed IPs of the two integration types
				BOOST_FOREACH(size_t ip, outsideIpsA)
					intTypePtrA->DeleteIntegrationPoint(ip);
				BOOST_FOREACH(size_t ip, outsideIpsB)
					intTypePtrB->DeleteIntegrationPoint(ip);
				//! @brief 5) impose the integration schemes
				newElPtrA->SetIntegrationType(intTypePtrA,NuTo::IpData::NOIPDATA);
				newElPtrB->SetIntegrationType(intTypePtrB,NuTo::IpData::NOIPDATA);
				//! @brief 6) delete the original element
				this->ElementDelete(thisCrackedElem->ElementGetId());
				//! if edge is cracked, the new element is introduced and the integration cell updated: Go to next cracked element
				this->NodeDelete(tmpNodeId,false);

//~ this->ElementInfo(newElPtrA,5);
//~ this->ElementInfo(newElPtrB,5);
				break; //!< break outer for-loop (element is already cracked)

			}
			// delete the temporary point
			this->NodeDelete(tmpNodeId,false);
		}
    }
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::Structure::InitiatePhantomNodeMethod] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}


//! @brief ... take cracked elements and initiate PhantomNodeMethod
//! @param rNumIp (Input) ... number of integration points for the new (cracked) elements
//! @return  ... id vector of cracked elements
NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> NuTo::Structure::InitiatePhantomNodeMethod(int rNumIp)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    elementBasePtrSet_t crackedElems, elemsToCrack;
	nodeBasePtrMap_t crackedNodes;

	this->InitiateCracks(elemsToCrack);

    //! crack the elements
    BOOST_FOREACH(NuTo::Structure::elementBasePtr_t thisCrackedElem, elemsToCrack)
    {
		//! at first only with one Crack!!
		NuTo::Structure::crackBasePtrVec_t crackPtrVec(thisCrackedElem->GetCracks());
		if(crackPtrVec.size()!=1)
			throw NuTo::MechanicsException("[NuTo::Structure::InitiatePhantomNodeMethod] Only 1 crack per element possible!!");
		NuTo::Structure::crackBasePtr_t thisCrackPtr=crackPtrVec[0];
		//! now find the cracked edge and the position of the crack
		//! assumption: linear elements (incremental numeration leads to the edges)
		const unsigned short numElemNodes=thisCrackedElem->GetNumNodes ();

		/** @brief the vector for the node pointers
		 * build up the nodepointer vectors
		 * initialize with the original node pointers
		 */
		nodeBasePtrVec_t oldNodes(numElemNodes), newNodes;
		for(unsigned short i=numElemNodes; i--;)
		{
			oldNodes[i]=thisCrackedElem->GetNode(i);
			newNodes=oldNodes;
		}
		unsigned short locNodeId1stEdge=0, locNodeId2ndEdge=0;
		for(locNodeId1stEdge=0; locNodeId1stEdge<numElemNodes; locNodeId1stEdge++)
		{
			// create a temporary point
			int tmpNodeId=this->NodeCreate("Coordinates");
			double relCoor=0.0;
			size_t segment;
			bool edgeIntersected=thisCrackPtr->Intersect(
					thisCrackedElem->GetNode(locNodeId1stEdge),
					thisCrackedElem->GetNode((locNodeId1stEdge+1)%numElemNodes),
					this->NodeGetNodePtr(tmpNodeId), relCoor, segment);
			if(edgeIntersected)
			{
				/**
				 * first build up new integration types
				 * 1) building up crack segment in element
				 * 2) create new nodes
				 * 3) create the two new elements
				 * 4) get number and position of integration points and create two new modifiable integration types
				 * 5) impose the integration scheme
				 * 6) delete the original element
				 */

				//! building up crack segment in element
				//! build up integration cells for the two elements
				double xA[2], xB[2];
				this->NodeGetNodePtr(tmpNodeId)->GetCoordinates2D(xA);
				//! check if the crack running thru the element
				unsigned short numIntersectedEdges=1; //!< one edge already intersected
				for(locNodeId2ndEdge=locNodeId1stEdge+1; locNodeId2ndEdge<numElemNodes; locNodeId2ndEdge++)
				{
					if(thisCrackPtr->Intersect(
							thisCrackedElem->GetNode(locNodeId2ndEdge),
							thisCrackedElem->GetNode(((locNodeId2ndEdge+1)%numElemNodes)),
							this->NodeGetNodePtr(tmpNodeId), relCoor, segment))
					{
						numIntersectedEdges++;
						break; //!< now we have two intersections --> it's enough
					}
				}
				//! check intersection type
				/**
				 * there are two cases for the cracked elements:
				 * 1) Elements which are torn in two pieces
				 * 2) Elements containing the cracktip
				 *
				 * Therefor we need 4 new nodes for case 1).
				 * For the 2nd case we need only two new nodes at the cracked edge.
				 */

				//! create the first two new nodes
				//! Node A
				nodeBasePtr_t oldNodePtr=thisCrackedElem->GetNode(locNodeId1stEdge);
				nodeBasePtr_t newNodePtr=NULL;
				nodeBasePtrMap_t::iterator crackedNodePairIt = crackedNodes.find(oldNodePtr);
				if(crackedNodePairIt==crackedNodes.end())
				{
					const int newNode=this->NodeCreate("displacements");
					newNodePtr=this->NodeGetNodePtr(newNode);
					double val2D[2];
					oldNodePtr->GetCoordinates2D(val2D);
					newNodePtr->SetCoordinates2D(val2D);
					oldNodePtr->GetDisplacements2D(val2D);
					newNodePtr->SetDisplacements2D(val2D);
					crackedNodes.insert(std::make_pair(oldNodePtr,newNodePtr));
				}else{
					newNodePtr=crackedNodePairIt->second;
				}
				newNodes[locNodeId1stEdge]=newNodePtr;
				//! Node B
				oldNodePtr=thisCrackedElem->GetNode(((locNodeId1stEdge+1)%numElemNodes));
				crackedNodePairIt = crackedNodes.find(oldNodePtr);
				if(crackedNodePairIt==crackedNodes.end())
				{
					const int newNode=this->NodeCreate("displacements");
					newNodePtr=this->NodeGetNodePtr(newNode);
					double val2D[2];
					oldNodePtr->GetCoordinates2D(val2D);
					newNodePtr->SetCoordinates2D(val2D);
					oldNodePtr->GetDisplacements2D(val2D);
					newNodePtr->SetDisplacements2D(val2D);
					crackedNodes.insert(std::make_pair(oldNodePtr,newNodePtr));
				}else{
					newNodePtr=crackedNodePairIt->second;
				}
				newNodes[(locNodeId1stEdge+1)%numElemNodes]=newNodePtr;

				if(numIntersectedEdges==2){ //!< if  numIntersectedEdges==2 -> crack running thru
					this->NodeGetNodePtr(tmpNodeId)->GetCoordinates2D(xB);
					//! for the torn element we have to copy all other element's nodes
					for(unsigned short i=numElemNodes; i--;)
					{
						if(i==locNodeId1stEdge || i==(locNodeId1stEdge+1)%numElemNodes) continue;
						oldNodePtr=thisCrackedElem->GetNode(i);
						crackedNodePairIt = crackedNodes.find(oldNodePtr);
						if(crackedNodePairIt==crackedNodes.end())
						{
							const int newNode=this->NodeCreate("displacements");
							newNodePtr=this->NodeGetNodePtr(newNode);
							double val2D[2];
							oldNodePtr->GetCoordinates2D(val2D);
							newNodePtr->SetCoordinates2D(val2D);
							oldNodePtr->GetDisplacements2D(val2D);
							newNodePtr->SetDisplacements2D(val2D);
							crackedNodes.insert(std::make_pair(oldNodePtr,newNodePtr));
						}else{
							newNodePtr=crackedNodePairIt->second;
						}
						newNodes[i]=newNodePtr;
					}
				}else if(numIntersectedEdges==1){ //!< if  numIntersectedEdges==1 -> crack ends inside the element
					bool edgeCracked=false;
					//! find the projection of the crack to the element corners.
					for(locNodeId2ndEdge=0; locNodeId2ndEdge<numElemNodes; locNodeId2ndEdge++)
					{
						if(locNodeId2ndEdge==locNodeId1stEdge) continue; //!< this edge is already cracked
						int newCrackEnd=this->NodeCreate("Coordinates");
						//! @todo check if the intersecting crack segment is at the end or beginning of the crack
						const unsigned short rayIntersect=thisCrackPtr->ExtendEnd(
								thisCrackedElem->GetNode(locNodeId2ndEdge),
								thisCrackedElem->GetNode(((locNodeId2ndEdge+1)%numElemNodes)),
								this->NodeGetNodePtr(newCrackEnd), relCoor);
						if(rayIntersect)
						{
							//! check the side of the crack: line segment AB is intersected at (0=no intersection, 1=front, 2=end)
							if(rayIntersect==2)
							{
								thisCrackPtr->PushBack(this->NodeGetNodePtr(newCrackEnd));
								//~ this->NodeGetNodePtr(newCrackEnd)->GetCoordinates2D(xB);
							}else if(rayIntersect==1)
							{
								//! in this case we add a new Node to the front of the crack
								//! and for consistency we change the direction of the crack xA->xB
								thisCrackPtr->PushFront(this->NodeGetNodePtr(newCrackEnd));
								//~ unsigned short i=0;
								//~ BOOST_FOREACH(double d, xA) xB[i++]=d;
								//~ this->NodeGetNodePtr(newCrackEnd)->GetCoordinates2D(xA);
							}else
								throw NuTo::MechanicsException("[NuTo::Structure::InitiatePhantomNodeMethod] Got a wrong end of the crack!!!");
							this->NodeGetNodePtr(newCrackEnd)->GetCoordinates2D(xB);
							edgeCracked=true;
							break; //!< now we have two intersections --> it's enough
						}else{
							this->NodeDelete(newCrackEnd,false);
						}
					}
					if(!edgeCracked)
						throw NuTo::MechanicsException("[NuTo::Structure::InitiatePhantomNodeMethod] Edge is not parted by end ray of the crack");
				}else{
					throw NuTo::MechanicsException("[NuTo::Structure::InitiatePhantomNodeMethod] Unhandeled case!!");
				}

				//! @brief 3) create the two new elements
				/**
				 * now we have two intersections of the element's corners.
				 * With this intersection nodes we build up the integration cells,
				 * therefore we check which integration points are outside this area.
				 */
				//! calculate the normal to the cracksegment
				//! @todo get this from a geometry class
				/*!
					In the two-dimensional case, the normal $\boldsymbol{n}$ of the crack is defined as
					\f[\boldsymbol{n} =
						 \dfrac{ \left( \boldsymbol{x}_B - \boldsymbol{x}_A \right) \times \boldsymbol{e}_z}
							   { \left| \left( \boldsymbol{x}_B - \boldsymbol{x}_A \right) \times \boldsymbol{e}_z \right|}
					\f]
				*/
				double normalVec[2] = { xB[1] - xA[1] , xA[0] - xB[0]};
				double norm = sqrt(normalVec[0]*normalVec[0]+normalVec[1]*normalVec[1]);
				normalVec[0] /= norm;
				normalVec[1] /= norm;

				//! @todo Remove this dirty hack --> make it general!!! (perhaps element copy constructor)
				std::vector<nodeBasePtr_t> rNodeVectorA(oldNodes);
				std::vector<nodeBasePtr_t> rNodeVectorB(newNodes);
				//! check the position of the nodes
				for(unsigned short i=oldNodes.size(); i--;)
				{
					double val2D[2];
					oldNodes[i]->GetCoordinates2D(val2D);
					//! calculate the distance of the node to the cracksegment
					/**
						In the two-dimensional case, with the normalized normal $\boldsymbol{n}$ of the crack the distance $d$ is defined as
						\f[ d = \left( \boldsymbol{x}_B - \boldsymbol{x}_A \right) \cdot \boldsymbol{n}
						\f]
					*/
					double dist= (val2D[0] - xA[0])*normalVec[0] +(val2D[1] - xA[1])*normalVec[1];
					/**
					 *  for this first element all IP on the left side are required
					 *  --> all IPs with positive distances have to be deleted in this element
					 *  	--> all IPs without signbit
					 *  --> all other IPs have to be deleted for the second element
					 *  first collect all not needed IPs (don't delete it now, otherwise you will get problems with the iterators inside the integration type)
					 */
					if(!std::signbit(dist))
					{
						rNodeVectorA[i]=newNodes[i];
						rNodeVectorB[i]=oldNodes[i];
					}
				}
				//! create the new element A
				const size_t numElA = this->ElementCreate(thisCrackedElem->GetEnumType(),rNodeVectorA,thisCrackedElem->GetElementDataType(), NuTo::IpData::NOIPDATA);
				elementBasePtr_t newElPtrA = this->ElementGetElementPtr(numElA);
				crackedElems.insert(newElPtrA);
				//! create the new element B
				const size_t numElB = this->ElementCreate(thisCrackedElem->GetEnumType(),rNodeVectorB,thisCrackedElem->GetElementDataType(), NuTo::IpData::NOIPDATA);
				elementBasePtr_t newElPtrB = this->ElementGetElementPtr(numElB);
				crackedElems.insert(newElPtrB);

				//! @brief 4) get numer and position of integration points and create two new modifiable integration types
				double ipCoor[3];
				std::stringstream  intTypeStr, intTypeStrA, intTypeStrB;
				//! set integration type for the new elements
				switch(thisCrackedElem->GetEnumType())
				{
				case NuTo::Element::PLANE2D3N:
					intTypeStr << "2D3NGAUSS3IP";
					break;
				case NuTo::Element::PLANE2D4N:
					intTypeStr << "2D4NMOD" << rNumIp << "IP";
					break;
				default:
					throw NuTo::MechanicsException("[" + std::string(__PRETTY_FUNCTION__) + "] Wrong Elementtype!");
				}
				//! @brief 6) delete the original element
				this->ElementDelete(thisCrackedElem->ElementGetId());

				//! build up new integration types
				intTypeStrA << intTypeStr.str() << numElA;
				intTypeStrB << intTypeStr.str() << numElB;
				NuTo::IntegrationTypeBase* intTypePtrA(this->GetPtrIntegrationType(intTypeStrA.str()));
				NuTo::IntegrationTypeBase* intTypePtrB(this->GetPtrIntegrationType(intTypeStrB.str()));
				//! @brief 5) impose the integration schemes
				newElPtrA->SetIntegrationType(intTypePtrA,NuTo::IpData::NOIPDATA);
				newElPtrB->SetIntegrationType(intTypePtrB,NuTo::IpData::NOIPDATA);

				//! separate IP's into two integration types
				std::vector<size_t> outsideIpsA(0), outsideIpsB(0);
				//! get the global coordinates of each IP
				size_t numIp=newElPtrA->GetNumIntegrationPoints();
				for(size_t ip=0; ip<numIp; ++ip){
					newElPtrA->GetGlobalIntegrationPointCoordinates(ip,ipCoor);
					//! calculate the distance of the node to the cracksegment
					/**
						In the two-dimensional case, with the normalized normal $\boldsymbol{n}$ of the crack the distance $d$ is defined as
						\f[ d = \left( \boldsymbol{x}_B - \boldsymbol{x}_A \right) \cdot \boldsymbol{n}
						\f]
					*/
					double dist= (ipCoor[0] - xA[0])*normalVec[0] +(ipCoor[1] - xA[1])*normalVec[1];
					/**
					 *  for this first element all IP on the left side are required
					 *  --> all IPs with positive distances have to be deleted in this element
					 *  	--> all IPs without signbit
					 *  --> all other IPs have to be deleted for the second element
					 *  first collect all not needed IPs (don't delete it now, otherwise you will get problems with the iterators inside the integration type)
					 */
					if(!std::signbit(dist))	outsideIpsA.push_back(ip);
						else 				outsideIpsB.push_back(ip);
				}
				//! delete non-needed IPs of the two integration types
				BOOST_FOREACH(size_t ip, outsideIpsA)
					intTypePtrA->DeleteIntegrationPoint(ip);
				BOOST_FOREACH(size_t ip, outsideIpsB)
					intTypePtrB->DeleteIntegrationPoint(ip);
				//! if edge is cracked, the new element is introduced and the integration cell updated: Go to next cracked element
				this->NodeDelete(tmpNodeId,false);

//~ this->ElementInfo(newElPtrA,5);
//~ this->ElementInfo(newElPtrB,5);
				break; //!< break outer for-loop (element is already cracked)

			}
			// delete the temporary point
			this->NodeDelete(tmpNodeId,false);
		}
    }
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::Structure::InitiatePhantomNodeMethod] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif

	std::vector<int> returnVec(0);
	returnVec.reserve(crackedElems.size());
    BOOST_FOREACH(elementBasePtr_t thisElem, crackedElems)
		returnVec.push_back(thisElem->ElementGetId());
    return NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic>(returnVec);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief ... take cracked elements and initiate PhantomNodeMethod
//! @param rNumIp (Input) ... number of integration points for the new (cracked) elements
//! @return  ... id vector of cracked elements
NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> NuTo::Structure::InitiatePhantomNodeMethodTriangle(int rNumIp)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    elementBasePtrSet_t crackedElems, elemsToCrack;
	nodeBasePtrMap_t crackedNodes;

	this->InitiateCracks(elemsToCrack);
	
	//! prepare an integration type for the tip elements
	std::vector<double> coords(2);
	std::vector< std::vector<double> > localCoordinatesArea;
	coords={-1,-1};
	localCoordinatesArea.push_back(coords);
	coords={1,-1};
	localCoordinatesArea.push_back(coords);
	coords={1,1};
	localCoordinatesArea.push_back(coords);
	coords={-1,1};
	localCoordinatesArea.push_back(coords);
	NuTo::IntegrationTypeBase* intTypePtr2D4NTRI(this->GetPtrIntegrationType("2D4NTRI"));
	intTypePtr2D4NTRI->AddIntegrationPoints(localCoordinatesArea,rNumIp);


    //! crack the elements
    BOOST_FOREACH(NuTo::Structure::elementBasePtr_t thisCrackedElem, elemsToCrack)
    {
		//! at first only with one Crack!!
		NuTo::Structure::crackBasePtrVec_t crackPtrVec(thisCrackedElem->GetCracks());
		if(crackPtrVec.size()!=1)
			throw NuTo::MechanicsException("[NuTo::Structure::InitiatePhantomNodeMethodTriangle] Only 1 crack per element possible!!");
		NuTo::Structure::crackBasePtr_t thisCrackPtr=crackPtrVec[0];
		//! now find the cracked edge and the position of the crack
		//! assumption: linear elements (incremental numeration leads to the edges)
		const unsigned short numElemNodes=thisCrackedElem->GetNumNodes ();

		/** @brief the vector for the node pointers
		 * build up the nodepointer vectors
		 * initialize with the original node pointers
		 */
		nodeBasePtrVec_t oldNodes(numElemNodes), newNodes;
		for(unsigned short i=numElemNodes; i--;)
		{
			oldNodes[i]=thisCrackedElem->GetNode(i);
			newNodes=oldNodes;
		}
		unsigned short locNodeId1stEdge=0, locNodeId2ndEdge=0;
		for(locNodeId1stEdge=0; locNodeId1stEdge<numElemNodes; locNodeId1stEdge++)
		{
			// create a temporary point
			int tmpNodeId=this->NodeCreate("Coordinates");
			double relCoor=0.0;
			size_t segment;
			bool edgeIntersected=thisCrackPtr->Intersect(
					thisCrackedElem->GetNode(locNodeId1stEdge),
					thisCrackedElem->GetNode((locNodeId1stEdge+1)%numElemNodes),
					this->NodeGetNodePtr(tmpNodeId), relCoor, segment);
			if(edgeIntersected)
			{
				/**
				 * first build up new integration types
				 * 1) building up crack segment in element
				 * 2) create new nodes
				 * 3) create the two new elements
				 * 4) get number and position of integration points and create two new modifiable integration types
				 * 5) impose the integration scheme
				 * 6) delete the original element
				 */

				//! building up crack segment in element
				//! build up integration cells for the two elements
				double xA[2], xB[2];
				this->NodeGetNodePtr(tmpNodeId)->GetCoordinates2D(xA);
				//! check if the crack running thru the element
				unsigned short numIntersectedEdges=1; //!< one edge already intersected
				for(locNodeId2ndEdge=locNodeId1stEdge+1; locNodeId2ndEdge<numElemNodes; locNodeId2ndEdge++)
				{
					if(thisCrackPtr->Intersect(
							thisCrackedElem->GetNode(locNodeId2ndEdge),
							thisCrackedElem->GetNode(((locNodeId2ndEdge+1)%numElemNodes)),
							this->NodeGetNodePtr(tmpNodeId), relCoor, segment))
					{
						numIntersectedEdges++;
						break; //!< now we have two intersections --> it's enough
					}
				}
				//! check intersection type
				/**
				 * there are two cases for the cracked elements:
				 * 1) Elements which are torn in two pieces
				 * 2) Elements containing the cracktip
				 *
				 * Therefor we need 4 new nodes for case 1).
				 * For the 2nd case we need only two new nodes at the cracked edge.
				 */
				//! create the first two new nodes
				//! Node A
				nodeBasePtr_t oldNodePtr=thisCrackedElem->GetNode(locNodeId1stEdge);
				nodeBasePtr_t newNodePtr=NULL;
				nodeBasePtrMap_t::iterator crackedNodePairIt = crackedNodes.find(oldNodePtr);
				if(crackedNodePairIt==crackedNodes.end())
				{
					const int newNode=this->NodeCreate("displacements");
					newNodePtr=this->NodeGetNodePtr(newNode);
					double val2D[2];
					oldNodePtr->GetCoordinates2D(val2D);
					newNodePtr->SetCoordinates2D(val2D);
					oldNodePtr->GetDisplacements2D(val2D);
					newNodePtr->SetDisplacements2D(val2D);
					crackedNodes.insert(std::make_pair(oldNodePtr,newNodePtr));
				}else{
					newNodePtr=crackedNodePairIt->second;
				}
				newNodes[locNodeId1stEdge]=newNodePtr;
				//! Node B
				oldNodePtr=thisCrackedElem->GetNode(((locNodeId1stEdge+1)%numElemNodes));
				crackedNodePairIt = crackedNodes.find(oldNodePtr);
				if(crackedNodePairIt==crackedNodes.end())
				{
					const int newNode=this->NodeCreate("displacements");
					newNodePtr=this->NodeGetNodePtr(newNode);
					double val2D[2];
					oldNodePtr->GetCoordinates2D(val2D);
					newNodePtr->SetCoordinates2D(val2D);
					oldNodePtr->GetDisplacements2D(val2D);
					newNodePtr->SetDisplacements2D(val2D);
					crackedNodes.insert(std::make_pair(oldNodePtr,newNodePtr));
				}else{
					newNodePtr=crackedNodePairIt->second;
				}
				newNodes[(locNodeId1stEdge+1)%numElemNodes]=newNodePtr;

				if(numIntersectedEdges==2){ //!< if  numIntersectedEdges==2 -> crack running thru
					this->NodeGetNodePtr(tmpNodeId)->GetCoordinates2D(xB);
					//! for the torn element we have to copy all other element's nodes
					for(unsigned short i=numElemNodes; i--;)
					{
						if(i==locNodeId1stEdge || i==(locNodeId1stEdge+1)%numElemNodes) continue;
						oldNodePtr=thisCrackedElem->GetNode(i);
						crackedNodePairIt = crackedNodes.find(oldNodePtr);
						if(crackedNodePairIt==crackedNodes.end())
						{
							const int newNode=this->NodeCreate("displacements");
							newNodePtr=this->NodeGetNodePtr(newNode);
							double val2D[2];
							oldNodePtr->GetCoordinates2D(val2D);
							newNodePtr->SetCoordinates2D(val2D);
							oldNodePtr->GetDisplacements2D(val2D);
							newNodePtr->SetDisplacements2D(val2D);
							crackedNodes.insert(std::make_pair(oldNodePtr,newNodePtr));
						}else{
							newNodePtr=crackedNodePairIt->second;
						}
						newNodes[i]=newNodePtr;
					}
				}else if(numIntersectedEdges==1){ //!< if  numIntersectedEdges==1 -> crack ends inside the element
					bool edgeCracked=false;
					nodeBasePtr_t crackTipNodePtr[2];
					//! find the projection of the crack to the element corners.
					for(locNodeId2ndEdge=0; locNodeId2ndEdge<numElemNodes; locNodeId2ndEdge++)
					{
						if(locNodeId2ndEdge==locNodeId1stEdge) continue; //!< this edge is already cracked
						int newCrackEnd=this->NodeCreate("Coordinates");
						//! store the two line nodes to find find the connected elements
						crackTipNodePtr[0]=thisCrackedElem->GetNode(locNodeId2ndEdge);
						crackTipNodePtr[1]=thisCrackedElem->GetNode(((locNodeId2ndEdge+1)%numElemNodes));
						//! @todo check if the intersecting crack segment is at the end or beginning of the crack
						const unsigned short rayIntersect=thisCrackPtr->ExtendEnd(
								crackTipNodePtr[0],crackTipNodePtr[1],
								this->NodeGetNodePtr(newCrackEnd), relCoor);
						if(rayIntersect)
						{
							//! check the side of the crack: line segment AB is intersected at (0=no intersection, 1=front, 2=end)
							if(rayIntersect==2)
							{
								thisCrackPtr->PushBack(this->NodeGetNodePtr(newCrackEnd));
								//~ this->NodeGetNodePtr(newCrackEnd)->GetCoordinates2D(xB);
							}else if(rayIntersect==1)
							{
								//! in this case we add a new Node to the front of the crack
								//! and for consistency we change the direction of the crack xA->xB
								thisCrackPtr->PushFront(this->NodeGetNodePtr(newCrackEnd));
								//~ unsigned short i=0;
								//~ BOOST_FOREACH(double d, xA) xB[i++]=d;
								//~ this->NodeGetNodePtr(newCrackEnd)->GetCoordinates2D(xA);
							}else
								throw NuTo::MechanicsException("[NuTo::Structure::InitiatePhantomNodeMethodTriangle] Got a wrong end of the crack!!!");
							this->NodeGetNodePtr(newCrackEnd)->GetCoordinates2D(xB);
							edgeCracked=true;
							//! change the integration scheme of the elements around the crack tip
							std::vector< elementBasePtr_t > elementVec;
							std::set< elementBasePtr_t > tipElements;
							for(size_t iTipNode=0;iTipNode<2;++iTipNode){
								this->NodeGetElements (crackTipNodePtr[iTipNode], elementVec);
								for(std::vector< elementBasePtr_t >::iterator it=elementVec.begin();it!=elementVec.end();++it){
									tipElements.insert(*it);
									DBG_PRINT_VAL(NodeGetId(crackTipNodePtr[iTipNode]));
									DBG_PRINT_VAL(ElementGetId((*it)));
								}
							}
							for(std::set< elementBasePtr_t >::iterator it=tipElements.begin();it!=tipElements.end();++it){
								//! set integration type for the new elements
								switch((*it)->GetEnumType())
								{
								case NuTo::Element::PLANE2D4N:
									(*it)->SetIntegrationType(intTypePtr2D4NTRI,NuTo::IpData::NOIPDATA);
									break;
								default:
									throw NuTo::MechanicsException("[" + std::string(__PRETTY_FUNCTION__) + "] Wrong Elementtype!");
								}
							}
							break; //!< now we have two intersections --> it's enough
						}else{
							this->NodeDelete(newCrackEnd,false);
						}
					}
					if(!edgeCracked)
						throw NuTo::MechanicsException("[NuTo::Structure::InitiatePhantomNodeMethodTriangle] Edge is not parted by end ray of the crack");
				}else{
					throw NuTo::MechanicsException("[NuTo::Structure::InitiatePhantomNodeMethodTriangle] Unhandeled case!!");
				}
				//! @brief 3) create the two new elements
				/**
				 * now we have two intersections of the element's corners.
				 * With this intersection nodes we build up the integration cells,
				 * therefore we check which integration points are outside this area.
				 */
				//! calculate the normal to the cracksegment
				//! @todo get this from a geometry class
				/*!
					In the two-dimensional case, the normal $\boldsymbol{n}$ of the crack is defined as
					\f[\boldsymbol{n} =
						 \dfrac{ \left( \boldsymbol{x}_B - \boldsymbol{x}_A \right) \times \boldsymbol{e}_z}
							   { \left| \left( \boldsymbol{x}_B - \boldsymbol{x}_A \right) \times \boldsymbol{e}_z \right|}
					\f]
				*/
				double normalVec[2] = { xB[1] - xA[1] , xA[0] - xB[0]};
				double norm = sqrt(normalVec[0]*normalVec[0]+normalVec[1]*normalVec[1]);
				normalVec[0] /= norm;
				normalVec[1] /= norm;

				//! @todo Remove this dirty hack --> make it general!!! (perhaps element copy constructor)
				std::vector<nodeBasePtr_t> rNodeVectorA(oldNodes);
				std::vector<nodeBasePtr_t> rNodeVectorB(newNodes);
				std::vector<double> nodeVec(2);
				std::vector< std::vector<double> > integrationAreaA, localCoordinatesAreaA, integrationAreaB, localCoordinatesAreaB;
				nodeVec.assign(xA,xA+2);
				integrationAreaB.push_back(nodeVec);
				nodeVec.assign(xB,xB+2);
				integrationAreaB.push_back(nodeVec);
				//! check the position of the nodes
				for(unsigned short i=oldNodes.size(); i--;)
				{
					double val2D[2];
					oldNodes[i]->GetCoordinates2D(val2D);
					//! calculate the distance of the node to the cracksegment
					/**
						In the two-dimensional case, with the normalized normal $\boldsymbol{n}$ of the crack the distance $d$ is defined as
						\f[ d = \left( \boldsymbol{x}_B - \boldsymbol{x}_A \right) \cdot \boldsymbol{n}
						\f]
					*/
					double dist= (val2D[0] - xA[0])*normalVec[0] +(val2D[1] - xA[1])*normalVec[1];
					/**
					 *  for this first element all IP on the left side are required
					 *  --> all IPs with positive distances have to be deleted in this element
					 *  	--> all IPs without signbit
					 *  --> all other IPs have to be deleted for the second element
					 *  first collect all not needed IPs (don't delete it now, otherwise you will get problems with the iterators inside the integration type)
					 */
					if(!std::signbit(dist))
					{
						rNodeVectorA[i]=newNodes[i];
						rNodeVectorB[i]=oldNodes[i];
						nodeVec.assign(val2D,val2D+2);
						integrationAreaB.push_back(nodeVec);
					}else{ //< left area corresponds to element A
						nodeVec.assign(val2D,val2D+2);
						integrationAreaA.push_back(nodeVec);
					}
				}
				nodeVec.assign(xB,xB+2);
				integrationAreaA.push_back(nodeVec);
				nodeVec.assign(xA,xA+2);
				integrationAreaA.push_back(nodeVec);
				
				//! create the new element A
				const size_t numElA = this->ElementCreate(thisCrackedElem->GetEnumType(),rNodeVectorA,thisCrackedElem->GetElementDataType(), NuTo::IpData::NOIPDATA);
				elementBasePtr_t newElPtrA = this->ElementGetElementPtr(numElA);
				crackedElems.insert(newElPtrA);
				//! create the new element B
				const size_t numElB = this->ElementCreate(thisCrackedElem->GetEnumType(),rNodeVectorB,thisCrackedElem->GetElementDataType(), NuTo::IpData::NOIPDATA);
				elementBasePtr_t newElPtrB = this->ElementGetElementPtr(numElB);
				crackedElems.insert(newElPtrB);

				//! @brief 4) get numer and position of integration points and create two new modifiable integration types
				std::stringstream  intTypeStr, intTypeStrA, intTypeStrB;
				//! set integration type for the new elements
				switch(thisCrackedElem->GetEnumType())
				{
				case NuTo::Element::PLANE2D3N:
					intTypeStr << "2D3NGAUSS3IP";
					break;
				case NuTo::Element::PLANE2D4N:
					//~ intTypeStr << "2D4NTRI" << rNumIp << "IP";
					intTypeStr << "2D4NTRI";
					break;
				default:
					throw NuTo::MechanicsException("[" + std::string(__PRETTY_FUNCTION__) + "] Wrong Elementtype!");
				}
				//! @brief 6) delete the original element
				this->ElementDelete(thisCrackedElem->ElementGetId());

				//! build up new integration types
				std::vector<double> ipCoords;
				//! a) domain A: get local coordinates of the bounding polygon and give to integration scheme subroutine
				intTypeStrA << intTypeStr.str() << numElA;
				NuTo::IntegrationTypeBase* intTypePtrA(this->GetPtrIntegrationType(intTypeStrA.str()));
				BOOST_FOREACH(std::vector<double> pt, integrationAreaA){
					double locCoord[2];
					double globCoord[2]={pt[0],pt[1]};
					if(newElPtrA->GetLocalPointCoordinates(globCoord,locCoord)){
						ipCoords.assign(locCoord,locCoord+2);
						localCoordinatesAreaA.push_back(ipCoords);
					}else
						throw NuTo::MechanicsException("[" + std::string(__PRETTY_FUNCTION__) + "] Integration point not within element!");
				}
				intTypePtrA->AddIntegrationPoints(localCoordinatesAreaA,rNumIp);
				//! b) domain B: get local coordinates of the bounding polygon and give to integration scheme subroutine
				intTypeStrB << intTypeStr.str() << numElB;
				NuTo::IntegrationTypeBase* intTypePtrB(this->GetPtrIntegrationType(intTypeStrB.str()));
				BOOST_FOREACH(std::vector<double> pt, integrationAreaB){
					double locCoord[2];
					double globCoord[2]={pt[0],pt[1]};
					if(newElPtrB->GetLocalPointCoordinates(globCoord,locCoord)){
						ipCoords.assign(locCoord,locCoord+2);
						localCoordinatesAreaB.push_back(ipCoords);
					}else
						throw NuTo::MechanicsException("[" + std::string(__PRETTY_FUNCTION__) + "] Integration point not within element!");
				}
				intTypePtrB->AddIntegrationPoints(localCoordinatesAreaB,rNumIp);

				//! @brief 5) impose the integration schemes
				newElPtrA->SetIntegrationType(intTypePtrA,NuTo::IpData::NOIPDATA);
				newElPtrB->SetIntegrationType(intTypePtrB,NuTo::IpData::NOIPDATA);

				//! if edge is cracked, the new element is introduced and the integration cell updated: Go to next cracked element
				this->NodeDelete(tmpNodeId,false);

				break; //!< break outer for-loop (element is already cracked)

			}
			// delete the temporary point
			this->NodeDelete(tmpNodeId,false);
		}
    }
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::Structure::InitiatePhantomNodeMethodTriangle] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif

	std::vector<int> returnVec(0);
	returnVec.reserve(crackedElems.size());
    BOOST_FOREACH(elementBasePtr_t thisElem, crackedElems)
		returnVec.push_back(thisElem->ElementGetId());
    return NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic>(returnVec);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

