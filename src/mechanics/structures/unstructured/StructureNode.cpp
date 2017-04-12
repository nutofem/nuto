#include <boost/tokenizer.hpp>
#include <sstream>

#include "base/Timer.h"


#include "mechanics/elements/ElementBase.h"
#include "mechanics/structures/StructureOutputBlockVector.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "math/SparseMatrixCSRGeneral.h"
#include "mechanics/groups/Group.h"
#include "mechanics/groups/GroupEnum.h"
#include "mechanics/nodes/NodeDof.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/structures/Assembler.h"

int NuTo::Structure::GetNumNodes() const
{
    return mNodeMap.size();
}

NuTo::NodeBase* NuTo::Structure::NodeGetNodePtr(int rIdent)
{
    boost::ptr_map<int,NodeBase>::iterator it = mNodeMap.find(rIdent);
    if (it!=mNodeMap.end())
        return it->second;
    else
    {
    	std::stringstream out;
    	out << rIdent;
    	throw MechanicsException("[NuTo::Structure::NodeGetNodePtr] Node with id " + out.str() +" does not exist.");
    }
}

const NuTo::NodeBase* NuTo::Structure::NodeGetNodePtr(int rIdent)const
{
    boost::ptr_map<int,NodeBase>::const_iterator it = mNodeMap.find(rIdent);
    if (it!=mNodeMap.end())
        return it->second;
    else
    {
    	std::stringstream out;
    	out << rIdent;
    	throw MechanicsException("[NuTo::Structure::NodeGetNodePtr] Node with id " + out.str() +" does not exist.");
    }
}

int NuTo::Structure::NodeGetId(const NodeBase* rNode)const
{
    for (boost::ptr_map<int,NodeBase>::const_iterator
            it = mNodeMap.begin(); it!= mNodeMap.end(); it++)
    {
        if (it->second==rNode)
            return it->first;
    }
    throw MechanicsException("[NuTo::Structure::GetNodeId] Node does not exist.");
}

const boost::ptr_map<int, NuTo::NodeBase>& NuTo::Structure::NodeGetNodeMap() const
{
    return mNodeMap;
}


void NuTo::Structure::NodeGetElements(const int rNodeId, std::vector<int>& rElementNumbers)
{
	const NuTo::NodeBase* nodePtr(NuTo::Structure::NodeGetNodePtr(rNodeId));
	std::vector<NuTo::ElementBase*> elementPtrs;
	this->NodeGetElements(nodePtr,elementPtrs);
	rElementNumbers.resize(elementPtrs.size());
	size_t i=0;
	BOOST_FOREACH(NuTo::ElementBase* thisElPtr,elementPtrs)
		rElementNumbers[++i] = ElementGetId(thisElPtr);
}

void NuTo::Structure::NodeGetElements(const NuTo::NodeBase* rNodePtr, std::vector<NuTo::ElementBase*>& rElements)
{
	rElements.resize(0);
	boost::ptr_map<int,NuTo::ElementBase>::iterator ElementIter = this->mElementMap.begin();
	while (ElementIter != this->mElementMap.end())
	{
		for(size_t thisNode=ElementIter->second->GetNumNodes();thisNode--;){
			NuTo::NodeBase* thisNodePtr=ElementIter->second->GetNode(thisNode);
			if(thisNodePtr==rNodePtr)
			{
				rElements.push_back(ElementIter->second);
				break;
			}
		}
		ElementIter++;
	}
}


void NuTo::Structure::NodeInfo(int rVerboseLevel)const
{
    mLogger <<"number of nodes   : " << mNodeMap.size() <<"\n";
    if (rVerboseLevel>3)
    {
        mLogger << "\t\tnodes :" <<"\n";
        for (boost::ptr_map<int,NodeBase>::const_iterator it = mNodeMap.begin(); it!= mNodeMap.end(); it++)
        {
            mLogger << "\t\t" << it->first;
            if (rVerboseLevel>4)
            {
                const NodeBase& node = *(it->second);
                auto dofTypes = node.GetDofTypes();

                for (Node::eDof dofType: dofTypes)
                {
                    mLogger << " \t \t" << Node::DofToString(dofType);
                    int numDofs = node.GetNum(dofType);
                    for (unsigned short iDof = 0; iDof < numDofs; ++iDof)
                    {
                        mLogger << " \t \t" << it->second->Get(dofType)[iDof];
                        if (it->second->IsDof(dofType))
                            mLogger << "("<< it->second->GetDof(dofType, iDof)<< ")" ;
                    }
                }
            }
            mLogger << "\n";
        }

    }
}

int NuTo::Structure::NodeCreate()
{
	Eigen::VectorXd coordinates(this->GetDimension());
	coordinates.setZero();

	//return int identifier of the new node
    return NodeCreate(coordinates);
}

int NuTo::Structure::NodeCreate(Eigen::VectorXd rCoordinates)
{
    int id = GetUnusedId(mNodeMap);
    this->NodeCreate(id, rCoordinates);
    return id;
}

int NuTo::Structure::NodeCreate(Eigen::VectorXd rCoordinates, std::set<NuTo::Node::eDof> rDofs)
{
    int id = GetUnusedId(mNodeMap);

    NodeBase* nodePtr = NodePtrCreate(rDofs, rCoordinates);

    // add node to map
    this->mNodeMap.insert(id, nodePtr);

    GetAssembler().SetNodeVectorChanged();

    //return int identifier of the new node
    return id;
}

int NuTo::Structure::GetDofDimension(Node::eDof rDof)
{
    switch (rDof)
    {
    // **************************************
    // dofs with dimension = global dimension
    // **************************************
    case Node::eDof::COORDINATES:
    case Node::eDof::DISPLACEMENTS:
        return GetDimension();

    // **************************************
    // scalars:
    // **************************************
    case Node::eDof::TEMPERATURE:
    case Node::eDof::NONLOCALEQSTRAIN:
    case Node::eDof::WATERVOLUMEFRACTION:
    case Node::eDof::RELATIVEHUMIDITY:
    case Node::eDof::CRACKPHASEFIELD:
    case Node::eDof::ELECTRICPOTENTIAL:
        return 1;

    // **************************************
    // others:
    // **************************************
    case Node::eDof::ROTATIONS:
    {
        if (GetDimension() == 2)
            return 1;

        if (GetDimension() == 3)
            return 3;

        throw MechanicsException(__PRETTY_FUNCTION__, "Rotations are only defined for structural dimension 2 and 3");
    }
    case Node::eDof::NONLOCALEQPLASTICSTRAIN:
    case Node::eDof::NONLOCALTOTALSTRAIN:
    {
        if (GetDimension() == 1)
            return 1;

        if (GetDimension() == 2)
            return 3;

        if (GetDimension() == 3)
            return 6;

        break;
    }
    default:
    break; // throw below gets called...
    }
    throw MechanicsException(__PRETTY_FUNCTION__, "Dimensions of the required DOF " + Node::DofToString(rDof) + " not defined.");
}

NuTo::NodeBase* NuTo::Structure::NodePtrCreate(std::set<Node::eDof> rDOFs, Eigen::VectorXd rCoordinates)
{

    if (rCoordinates.rows() != mDimension)
        throw MechanicsException(__PRETTY_FUNCTION__, "Dimension of the coordinate vector does not fit the dimension of the structure");

    NodeBase* nodePtr = nullptr;

    std::map<Node::eDof, NodeDofInfo> dofInfos;

    // somehow always add coordinates
    rDOFs.insert(Node::eDof::COORDINATES);

    for (Node::eDof dof : rDOFs)
    {
        NodeDofInfo& dofInfo = dofInfos[dof];

        dofInfo.mDimension = GetDofDimension(dof);
        dofInfo.mNumTimeDerivatives = GetNumTimeDerivatives();
        dofInfo.mIsDof = true;

        if (dof == Node::eDof::COORDINATES)
        {
            dofInfo.mNumTimeDerivatives = 0;
            dofInfo.mIsDof = false;
        }
    }

    nodePtr = new NodeDof(dofInfos);

    nodePtr->Set(Node::eDof::COORDINATES, rCoordinates);
    return nodePtr;
}

void NuTo::Structure::NodeCreate(int rNodeNumber, Eigen::VectorXd rCoordinates)
{
	// check node number
	boost::ptr_map<int,NodeBase>::iterator it = this->mNodeMap.find(rNodeNumber);
	if(it != this->mNodeMap.end())
	{
		throw MechanicsException("[NuTo::Structure::NodeCreate] Node already exists.");
	}

	std::set<Node::eDof> dofs;
    dofs.insert(Node::eDof::COORDINATES);

    NodeBase* nodePtr = NodePtrCreate(dofs, rCoordinates);

    // add node to map
    this->mNodeMap.insert(rNodeNumber, nodePtr);

    GetAssembler().SetNodeVectorChanged();
}

int NuTo::Structure::NodeCreateDOFs(std::string rDOFs)
{
    Eigen::VectorXd coordinates(this->GetDimension());
    coordinates.setZero();

    //return int identifier of the new node
    return NodeCreateDOFs(rDOFs, coordinates);
}

int NuTo::Structure::NodeCreateDOFs(std::string rDOFs, Eigen::VectorXd rCoordinates)
{
    int id = GetUnusedId(mNodeMap);
    NodeCreateDOFs(id, rDOFs, rCoordinates);
    return id;
}

void NuTo::Structure::NodeCreateDOFs(int rNodeNumber, std::string rDOFs, Eigen::VectorXd rCoordinates)
{

    // check node number
    boost::ptr_map<int,NodeBase>::iterator it = this->mNodeMap.find(rNodeNumber);
    if(it != this->mNodeMap.end())
    {
        throw MechanicsException("[NuTo::Structure::NodeCreateDOFs] Node already exists.");
    }

    std::set<Node::eDof> dofs;
    dofs.insert(Node::eDof::COORDINATES);

    // transform string to uppercase
    std::transform(rDOFs.begin(), rDOFs.end(), rDOFs.begin(), toupper);
    boost::char_separator<char> sep(" ");
    boost::tokenizer<boost::char_separator<char> > tok(rDOFs, sep);
    for (boost::tokenizer<boost::char_separator<char> >::iterator beg = tok.begin(); beg != tok.end(); ++beg)
    {
        try
        {
            std::string dofString = *beg;
            Node::eDof dofEnum = Node::DofToEnum(dofString);
            dofs.insert(dofEnum);
        }
        catch (NuTo::MechanicsException& e)
        {
            e.AddMessage("[NuTo::Structure::NodeCreate] invalid dof type: " + *beg + ".");
            throw;
        }
    }

    NodeBase* nodePtr = NodePtrCreate(dofs, rCoordinates);

    // add node to map
    this->mNodeMap.insert(rNodeNumber, nodePtr);

    GetAssembler().SetNodeVectorChanged();
}


int NuTo::Structure::NodeCreateDOFs(std::set<NuTo::Node::eDof> rDOFs)
{
    Eigen::VectorXd coordinates(this->GetDimension());
    coordinates.setZero();

    //return int identifier of the new node
    return NodeCreateDOFs(rDOFs, coordinates);
}


int NuTo::Structure::NodeCreateDOFs(std::set<NuTo::Node::eDof> rDOFs, Eigen::VectorXd rCoordinates)
{
    int id = GetUnusedId(mNodeMap);
    NodeCreateDOFs(id, rDOFs, rCoordinates);
    return id;
}


void NuTo::Structure::NodeCreateDOFs(int rNodeNumber, std::set<NuTo::Node::eDof> rDOFs, Eigen::VectorXd rCoordinates)
{
    // check node number
    boost::ptr_map<int,NodeBase>::iterator it = this->mNodeMap.find(rNodeNumber);
    if(it != this->mNodeMap.end())
    {
        throw MechanicsException("[NuTo::Structure::NodeCreateDOFs] Node already exists.");
    }

    rDOFs.insert(Node::eDof::COORDINATES);

    NodeBase* nodePtr = NodePtrCreate(rDOFs, rCoordinates);

    // add node to map
    this->mNodeMap.insert(rNodeNumber, nodePtr);

    GetAssembler().SetNodeVectorChanged();
}


std::vector<int> NuTo::Structure::NodesCreate(Eigen::MatrixXd& rCoordinates)
{
	std::vector<int> idVec;
	/// go through the nodes
	for(size_t i=0 ; i<(size_t)rCoordinates.cols(); ++i)
	{
		Eigen::VectorXd coordinate = rCoordinates.col(i);
		idVec.push_back(NodeCreate(coordinate));
	}

    return idVec;
}

void NuTo::Structure::NodeDelete(int rNodeNumber)
{
	NodeDelete(rNodeNumber,true);
}

void NuTo::Structure::NodeDelete(int rNodeNumber, bool checkElements)
{
    boost::ptr_map<int,NodeBase>::iterator itNode = mNodeMap.find(rNodeNumber);
    if (itNode==this->mNodeMap.end())
    {
        throw MechanicsException("[NuTo::Structure::NodeDelete] Node with the given id does not exist.");
    }
    else
    {
    	if (checkElements)
    	{
    		NodeBase* nodePtr = itNode->second;
    	 	/// Search for node in elements: using a loop over all elements
			for(boost::ptr_map<int,ElementBase>::const_iterator elemIt=mElementMap.begin(); elemIt!=mElementMap.end(); ++elemIt)
			{
				// loop over all element nodes
				for(unsigned short iNode=0; iNode<elemIt->second->GetNumNodes(); ++iNode)
				{
					// if the id of the node to be deleted is found, throw an error
					if(nodePtr == elemIt->second->GetNode(iNode) )
					{
						std::stringstream outNode, outElement;
						outNode << rNodeNumber;
						outElement << (elemIt->first);
						throw MechanicsException( "[NuTo::Structure::NodeDelete] Node " + outNode.str() + " is used by element " + outElement.str() + ", delete element first");
					}
				}
			}
    	}

    	                        for(boost::ptr_map<int,GroupBase>::iterator groupIt=mGroupMap.begin();groupIt!=mGroupMap.end(); ++groupIt){
        	if(groupIt->second->GetType()==NuTo::eGroupId::Nodes){
        		if(groupIt->second->Contain(itNode->first))
        		{
        			groupIt->second->RemoveMember(itNode->first);
        		}
        	}
        }

        // delete element from map
        this->mNodeMap.erase(itNode);

        GetAssembler().SetNodeVectorChanged();
    }
}


void NuTo::Structure::NodeBuildGlobalDofs(std::string rCallerName)
{
    if (not GetAssembler().RenumberingRequired())
        return;
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    try
    {
        UpdateDofStatus();
        std::vector<NodeBase*> nodes;
        this->GetNodesTotal(nodes);
        GetAssembler().BuildGlobalDofs(nodes); 
        UpdateDofStatus();
    }
    catch (MathException& e)
    {
        e.AddMessage(std::string("[") + __PRETTY_FUNCTION__ + "] Error building global dof numbering. \n");
        if (not rCallerName.empty())
            e.AddMessage(std::string(" -- was called by [") + rCallerName + "]");

        throw;
    }
    catch (MechanicsException& e)
    {
        e.AddMessage(std::string("[") + __PRETTY_FUNCTION__ + "] Error building global dof numbering. \n");
        if (not rCallerName.empty())
            e.AddMessage(std::string(" -- was called by [") + rCallerName + "]");

        throw;
    }
}




NuTo::StructureOutputBlockVector NuTo::Structure::NodeExtractDofValues(int rTimeDerivative) const
{
    GetAssembler().ThrowIfRenumberingRequred();

    StructureOutputBlockVector dofValues(GetDofStatus(), true); // with resize

    for (auto dofType : DofTypesGet())
    {
        auto& actDofValues = dofValues.J[dofType];
        auto& depDofValues = dofValues.K[dofType];

        int numActiveDofs = GetNumActiveDofs(dofType);

        for (auto it : mNodeMap)
        {
            const NodeBase& node = *it->second;
            if (not node.IsDof(dofType))
                continue;

            const auto& values = node.Get(dofType, rTimeDerivative);
            for (int i = 0; i < values.rows(); ++i)
            {
                int dofNumber = node.GetDof(dofType, i); 
                if (dofNumber < numActiveDofs)
                    actDofValues[dofNumber] = values[i]; 
                else
                    depDofValues[dofNumber - numActiveDofs] = values[i];
            }
        }
    }
    return dofValues;
}

void NuTo::Structure::NodeMergeDofValues(int rTimeDerivative, const NuTo::BlockFullVector<double>& rActiveDofValues, const NuTo::BlockFullVector<double>& rDependentDofValues)
{
    GetAssembler().ThrowIfRenumberingRequred();

    for (auto dofType : DofTypesGetActive())
    {
        if (rActiveDofValues[dofType].rows() != GetNumActiveDofs(dofType))
            throw MechanicsException(__PRETTY_FUNCTION__,  "invalid dimension of active dof vector for " + Node::DofToString(dofType));

        if (rDependentDofValues[dofType].rows() != GetNumDependentDofs(dofType))
            throw MechanicsException(__PRETTY_FUNCTION__, "invalid dimension of dependent dof vector for " + Node::DofToString(dofType));

        auto& actDofValues = rActiveDofValues[dofType];
        auto& depDofValues = rDependentDofValues[dofType];
        int numActiveDofs = GetNumActiveDofs(dofType);

        for (auto it : mNodeMap)
        {
            NodeBase& node = *it->second;
            if (not node.IsDof(dofType))
                continue;

            int numValues = node.GetNum(dofType);
            Eigen::VectorXd values(numValues);
            for (int i = 0; i < numValues; ++i)
            {
                int dofNumber = node.GetDof(dofType, i); 
                if (dofNumber < numActiveDofs)
                    values[i] = actDofValues[dofNumber]; 
                else
                    values[i] = depDofValues[dofNumber - numActiveDofs];
            }
            node.Set(dofType, rTimeDerivative, values);
        }

    }
    this->mUpdateTmpStaticDataRequired=true;
}




NuTo::BlockFullVector<double> NuTo::Structure::NodeCalculateDependentDofValues(const NuTo::BlockFullVector<double>& rActiveDofValues) const
{
    GetAssembler().ThrowIfRenumberingRequred();

    return this->GetAssembler().GetConstraintRhs() - GetAssembler().GetConstraintMatrix() * rActiveDofValues;
}


// store all nodes of a structure in a vector
void NuTo::Structure::GetNodesTotal(std::vector<const NodeBase*>& rNodes) const
{
	rNodes.reserve(mNodeMap.size());
    rNodes.resize(0);
	boost::ptr_map<int,NodeBase>::const_iterator NodeIter = this->mNodeMap.begin();
    while (NodeIter != this->mNodeMap.end())
    {
    	rNodes.push_back(NodeIter->second);
        NodeIter++;
    }
}

// store all nodes of a structure in a vector
void NuTo::Structure::GetNodesTotal(std::vector<std::pair<int, const NodeBase*> >& rNodes) const
{
	rNodes.reserve(mNodeMap.size());
    rNodes.resize(0);
	boost::ptr_map<int,NodeBase>::const_iterator NodeIter = this->mNodeMap.begin();
    while (NodeIter != this->mNodeMap.end())
    {
    	rNodes.push_back(std::pair<int, const NodeBase*>(NodeIter->first,NodeIter->second));
        NodeIter++;
    }
}

// store all nodes of a structure in a vector
void NuTo::Structure::GetNodesTotal(std::vector<NodeBase*>& rNodes)
{
	rNodes.reserve(mNodeMap.size());
	rNodes.resize(0);
    boost::ptr_map<int,NodeBase>::iterator NodeIter = this->mNodeMap.begin();
    while (NodeIter != this->mNodeMap.end())
    {
    	rNodes.push_back(NodeIter->second);
        NodeIter++;
    }
}

// store all nodes of a structure in a vector
void NuTo::Structure::GetNodesTotal(std::vector<std::pair<int,NodeBase*> >& rNodes)
{
	rNodes.reserve(mNodeMap.size());
	rNodes.resize(0);
    boost::ptr_map<int,NodeBase>::iterator NodeIter = this->mNodeMap.begin();
    while (NodeIter != this->mNodeMap.end())
    {
    	rNodes.push_back(std::pair<int, NodeBase*>(NodeIter->first,NodeIter->second));
        NodeIter++;
    }
}

void NuTo::Structure::NodeExchangePtr(int rId, NuTo::NodeBase* rOldPtr, NuTo::NodeBase* rNewPtr, std::vector<ElementBase*> rElements)
{
    //in node map
    //find it
    if (mNodeMap.erase(rId)!=1)
    {
        throw MechanicsException("[NuTo::Structure::NodeExchangePtr] Pointer to node (to exchange) does not exist.");
    }

    mNodeMap.insert(rId,rNewPtr);

    if (rElements.empty())
    {
        // in all elements
        for (boost::ptr_map<int,ElementBase>::iterator itElement = mElementMap.begin(); itElement!= mElementMap.end(); itElement++)
        {
            itElement->second->ExchangeNodePtr(rOldPtr, rNewPtr);
        }

    }
    else
    {
        // in specific elements:
        for (ElementBase* element : rElements)
        {
            element->ExchangeNodePtr(rOldPtr, rNewPtr);
        }
    }

    //in groups
    for(boost::ptr_map<int,GroupBase>::iterator groupIt=mGroupMap.begin();groupIt!=mGroupMap.end(); ++groupIt)
    {
        if(groupIt->second->GetType()==NuTo::eGroupId::Nodes)
        {
        	if (groupIt->second->Contain(rId))
                groupIt->second->ExchangePtr(rId, rOldPtr, rNewPtr);
        }
    }

    //in constraints
    GetAssembler().GetConstraints().ExchangeNodePtr(*rOldPtr, *rNewPtr);
}
