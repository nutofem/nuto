// $Id$

#include <sstream>

#include <boost/foreach.hpp>

#include "nuto/base/Debug.h"

#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/cracks/CrackExplicit2D.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/MechanicsException.h"

#include "nuto/mechanics/elements/ElementDataConstitutiveIpCrack.h"

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

	return id;
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
		unsigned int thisCrackId=ptrElementData->AddCrack(thisCrack);

		//! check if rCrackedElems contains rThisCrackedElems already, if not: append
		rCrackedElems.insert(thisElPtr);

		std::cout << "Element to be cracked: " << ElementGetId(thisElPtr)  << "; mElementData=" << ptrElementData
				<< "; mIntegrationType=" << ptrElementData->GetIntegrationType()
				<< "; mElementDataType=" << thisElPtr->GetElementDataType()
				<< "; thisCrackId=" << thisCrackId
		<< std::endl;

	}
}

//! @brief ... take cracked elements and initiate PhantomNodeMethod
//! @param rCrackedElems (Input) ... vector of cracked elements
void NuTo::Structure::InitiatePhantomNodeMethod(elementBasePtrSet_t & rCrackedElems)
{
    nodeBasePtrMap_t crackedNodes;
	nodeBasePtr_t oldNodePtrA, oldNodePtrB, oldNodePtrC, oldNodePtrD;
	int newNodeA(0), newNodeB(0), newNodeC(0), newNodeD(0);
	nodeBasePtr_t newNodePtrA=NULL, newNodePtrB=NULL, newNodePtrC=NULL, newNodePtrD=NULL;

	//! crack the elements
    BOOST_FOREACH(NuTo::Structure::elementBasePtr_t thisCrackedElem, rCrackedElems)
    {
		//! at first only with one Crack!!
		NuTo::Structure::crackBasePtrVec_t crackPtrVec(thisCrackedElem->GetCracks());
		NuTo::Structure::crackBasePtr_t thisCrackPtr=crackPtrVec[0];
		thisCrackPtr->Info(5);
		//! now find the cracked edge and the position of the crack
		//! assumption: linear elements (incremental numeration leads to the edges)
		const unsigned short numElemNodes=thisCrackedElem->GetNumNodes ();
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

				//! create the first two new nodes
				//! Node A
				oldNodePtrA=thisCrackedElem->GetNode(locNodeId1stEdge);
				nodeBasePtrMap_t::iterator crackedNodePairIt = crackedNodes.find(oldNodePtrA);
				if(crackedNodePairIt==crackedNodes.end())
				{
					newNodeA=this->NodeCreate("displacements");
					newNodePtrA=this->NodeGetNodePtr(newNodeA);
					double val2D[2];
					oldNodePtrA->GetCoordinates2D(val2D);
					newNodePtrA->SetCoordinates2D(val2D);
					oldNodePtrA->GetDisplacements2D(val2D);
					newNodePtrA->SetDisplacements2D(val2D);
					crackedNodes.insert(std::make_pair(oldNodePtrA,newNodePtrA));
				}else{
					newNodePtrA=crackedNodePairIt->second;
					newNodeA=this->NodeGetId(newNodePtrA);
				}

				//! Node B
				oldNodePtrB=thisCrackedElem->GetNode(((locNodeId1stEdge+1)%numElemNodes));
				crackedNodePairIt = crackedNodes.find(oldNodePtrB);
				if(crackedNodePairIt==crackedNodes.end())
				{
					newNodeB=this->NodeCreate("displacements");
					newNodePtrB=this->NodeGetNodePtr(newNodeB);
					double val2D[2];
					oldNodePtrB->GetCoordinates2D(val2D);
					newNodePtrB->SetCoordinates2D(val2D);
					oldNodePtrB->GetDisplacements2D(val2D);
					newNodePtrB->SetDisplacements2D(val2D);
					crackedNodes.insert(std::make_pair(oldNodePtrB,newNodePtrB));
				}else{
					newNodePtrB=crackedNodePairIt->second;
					newNodeB=this->NodeGetId(newNodePtrB);
				}


				/**
				 * first build up new integration types
				 * 1) building up crack segment in element
				 * 2) create two new nodes
				 * 3) create the two new elements
				 * 4) get numer and position of integration points and create two new modifiable integration types
				 * 5) impose the integration scheme
				 * 6) delete the original element
				 */
				//~ /**
				 //~ * @brief find neighbors and create new nodes, if necessary
				 //~ * find the elements containing nodes i and i+1
				 //~ */
				//~ elementBasePtrVec_t elements1(0), elements2(0);
				//~ elementBasePtrSet_t connectedElements;
				//~ this->NodeGetElements(thisCrackedElem->GetNode(i),elements1);
				//~ this->NodeGetElements(thisCrackedElem->GetNode(i),elements2);
				//~ BOOST_FOREACH(elementBasePtr_t thisElem1, elements1)
				//~ {
					//~ if(thisElem1==thisCrackedElem) continue;
					//~ BOOST_FOREACH(elementBasePtr_t thisElem2, elements2)
					//~ if(thisElem1==thisElem2) connectedElements.insert(thisElem1);
				//~ }
				//~
				//~ BOOST_FOREACH(elementBasePtr_t thisElem, connectedElements)
				//~ {
					//~ //! check if the neighbor element is cracked
					//~ DBG_PRINT_VAL(thisElem->ElementGetId())
					//~ DBG_PRINT_VAL(thisElem->IsCracked())
					//~ //! if element is not cracked yet: introduce new nodes
				//~ }

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
						//! create the second two new nodes
						//! Node C
						oldNodePtrC=thisCrackedElem->GetNode(locNodeId2ndEdge);
						crackedNodePairIt = crackedNodes.find(oldNodePtrC);
						if(crackedNodePairIt==crackedNodes.end())
						{
							newNodeC=this->NodeCreate("displacements");
							newNodePtrC=this->NodeGetNodePtr(newNodeC);
							double val2D[2];
							oldNodePtrC->GetCoordinates2D(val2D);
							newNodePtrC->SetCoordinates2D(val2D);
							oldNodePtrC->GetDisplacements2D(val2D);
							newNodePtrC->SetDisplacements2D(val2D);
							crackedNodes.insert(std::make_pair(oldNodePtrC,newNodePtrC));
						}else{
							newNodePtrC=crackedNodePairIt->second;
							newNodeC=this->NodeGetId(newNodePtrC);
						}

						//! Node D
						oldNodePtrD=thisCrackedElem->GetNode(((locNodeId2ndEdge+1)%numElemNodes));
						crackedNodePairIt = crackedNodes.find(oldNodePtrD);
						if(crackedNodePairIt==crackedNodes.end())
						{
							newNodeD=this->NodeCreate("displacements");
							newNodePtrD=this->NodeGetNodePtr(newNodeD);
							double val2D[2];
							oldNodePtrD->GetCoordinates2D(val2D);
							newNodePtrD->SetCoordinates2D(val2D);
							oldNodePtrD->GetDisplacements2D(val2D);
							newNodePtrD->SetDisplacements2D(val2D);
							crackedNodes.insert(std::make_pair(oldNodePtrD,newNodePtrD));
						}else{
							newNodePtrD=crackedNodePairIt->second;
							newNodeD=this->NodeGetId(newNodePtrD);
						}
						break; //!< now we have two intersections --> it's enough
					}
				}
				//! check intersection type
				if(numIntersectedEdges==2){ //< if  numIntersectedEdges==2 -> crack running thru
					this->NodeGetNodePtr(tmpNodeId)->GetCoordinates2D(xB);
				}else if(numIntersectedEdges==1){ //< if  numIntersectedEdges==1 -> crack ends inside the element
					//! find the projection of the crack to the element corners.
					for(locNodeId2ndEdge=0; locNodeId2ndEdge<numElemNodes; locNodeId2ndEdge++)
					{
						if(locNodeId2ndEdge==locNodeId1stEdge) continue; //< this edge is already cracked
						int newCrackEnd=this->NodeCreate("Coordinates");
						//! @todo check if the intersecting crack segment is at the end or beginning of the crack
						if(thisCrackPtr->ExtendEnd(
								thisCrackedElem->GetNode(locNodeId2ndEdge),
								thisCrackedElem->GetNode(((locNodeId2ndEdge+1)%numElemNodes)),
								this->NodeGetNodePtr(newCrackEnd), relCoor))
						{
							thisCrackPtr->PushBack(this->NodeGetNodePtr(newCrackEnd));
							this->NodeGetNodePtr(newCrackEnd)->GetCoordinates2D(xB);
							oldNodePtrC=thisCrackedElem->GetNode(locNodeId2ndEdge);
							newNodePtrC=oldNodePtrC;
							oldNodePtrD=thisCrackedElem->GetNode(((locNodeId2ndEdge+1)%numElemNodes));
							newNodePtrD=oldNodePtrD;
							break; //!< now we have two intersections --> it's enough
						}else{
							this->NodeDelete(newCrackEnd);
						}
					}
				}else{
					throw NuTo::MechanicsException("[NuTo::Structure::InitiatePhantomNodeMethod] Unhandeled case!!");
				}

				//! @brief 3) create the two new elements
				//! @todo Remove this dirty hack --> make it general!!! (perhaps element copy constructor)
				std::vector<nodeBasePtr_t> rNodeVector;
				//! 3a) elementA
DBG_PRINT_VAL(locNodeId1stEdge)
DBG_PRINT_VAL(locNodeId1stEdge+1)
DBG_PRINT_VAL(locNodeId2ndEdge)
DBG_PRINT_VAL(locNodeId2ndEdge+1)
				rNodeVector=std::vector<nodeBasePtr_t>(0);
				for(int vertex=0; vertex<numElemNodes; vertex++)
				{
DBG_PRINT_VAL(vertex)
					if(vertex==(locNodeId1stEdge)%numElemNodes)
						rNodeVector.push_back(oldNodePtrA);
					else if(vertex==(locNodeId1stEdge+1)%numElemNodes)
						rNodeVector.push_back(newNodePtrB);
					else if(vertex==(locNodeId2ndEdge)%numElemNodes)
						rNodeVector.push_back(newNodePtrC);
					else if(vertex==(locNodeId2ndEdge+1)%numElemNodes)
						rNodeVector.push_back(oldNodePtrD);
				}

				BOOST_FOREACH(nodeBasePtr_t thisNode, rNodeVector)
					DBG_PRINT_VAL(NodeGetId(thisNode))

				//! create the new element A
				const size_t numElA = this->ElementCreate(thisCrackedElem->GetEnumType(),rNodeVector,thisCrackedElem->GetElementDataType(), NuTo::IpData::NOIPDATA);
				elementBasePtr_t newElPtrA = this->ElementGetElementPtr(numElA);

				//! 3b) elementB: all nodes but the ith
				rNodeVector=std::vector<nodeBasePtr_t>(0);
				for(int vertex=0; vertex<numElemNodes; vertex++)
				{
					if(vertex==(locNodeId1stEdge)%numElemNodes)
						rNodeVector.push_back(newNodePtrA);
					else if(vertex==(locNodeId1stEdge+1)%numElemNodes)
						rNodeVector.push_back(oldNodePtrB);
					else if(vertex==(locNodeId2ndEdge)%numElemNodes)
						rNodeVector.push_back(oldNodePtrC);
					else if(vertex==(locNodeId2ndEdge+1)%numElemNodes)
						rNodeVector.push_back(newNodePtrD);
				}
				//! create the new element B
				const size_t numElB = this->ElementCreate(thisCrackedElem->GetEnumType(),rNodeVector,thisCrackedElem->GetElementDataType(), NuTo::IpData::NOIPDATA);
				elementBasePtr_t newElPtrB = this->ElementGetElementPtr(numElB);

				//! @brief 4) get numer and position of integration points and create two new modifiable integration types
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
					/// calculate the distance of the node to the cracksegment
					/*!
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
				this->NodeDelete(tmpNodeId);


				break; //!< break upper for-loop (element is already cracked)

			}
			// delete the temporary point
			this->NodeDelete(tmpNodeId);
		}
    }

}
