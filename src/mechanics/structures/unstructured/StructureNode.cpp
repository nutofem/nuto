#include <boost/tokenizer.hpp>
#include <sstream>

#include "base/Timer.h"


#include "mechanics/constraints/ConstraintBase.h"
#include "mechanics/elements/ElementBase.h"
#include "mechanics/structures/StructureOutputBlockVector.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "math/SparseMatrixCSRGeneral.h"
#include "mechanics/groups/Group.h"
#include "mechanics/groups/GroupEnum.h"
#include "mechanics/nodes/NodeDof.h"
#include "mechanics/nodes/NodeEnum.h"

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

    //renumbering of dofs for global matrices required
    this->mNodeNumberingRequired  = true;

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

    //renumbering of dofs for global matrices required
    this->mNodeNumberingRequired  = true;
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

    //renumbering of dofs for global matrices required
    this->mNodeNumberingRequired  = true;

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

    //renumbering of dofs for global matrices required
    this->mNodeNumberingRequired  = true;
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

        this->mNodeNumberingRequired = true;
    }
}


void NuTo::Structure::NodeBuildGlobalDofs(std::string rCallerName)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    try
    {
        std::map<Node::eDof, int> numDofsMap;

        // build initial node numbering

        //
        // number Lagrange multipliers in constraint equations defined in StructureBase
        // currently removed
        // ConstraintNumberGlobalDofs(this->mNumDofs);
        //

        UpdateDofStatus();

        this->mNodeNumberingRequired = false;

        BOOST_FOREACH(auto it, mNodeMap)
            it.second->SetGlobalDofsNumbers(numDofsMap);


        mConstraintMatrix.AllocateSubmatrices();
        mConstraintMappingRHS.AllocateSubmatrices();
        mConstraintRHS.AllocateSubvectors();


        for (auto dof : DofTypesGet())
        {
            mDofStatus.SetNumDependentDofs(dof, ConstraintGetNumLinearConstraints(dof));
            mDofStatus.SetNumActiveDofs(dof, numDofsMap[dof] - GetNumDependentDofs(dof));
        }



        // build constraint matrix for all dofs
        mConstraintMatrix = ConstraintGetConstraintMatrixBeforeGaussElimination();

        for (auto dof : DofTypesGet())
        {
            auto& constraintMatrix = mConstraintMatrix(dof, dof);

            const int numActiveDofs = GetNumActiveDofs(dof);
            const int numDependentDofs = GetNumDependentDofs(dof);
            const int numDofs = GetNumDofs(dof);

            //init RHSMatrix as a diagonal identity matrix
            auto& constraintMappingRHS = mConstraintMappingRHS(dof, dof);

            constraintMappingRHS.Resize(numDependentDofs, numDependentDofs);
            for (int i = 0; i < numDependentDofs ; ++i)
                constraintMappingRHS.AddValue(i,i,1.);

            // perform gauss algorithm
            std::vector<int> mappingInitialToNewOrdering;
            std::vector<int> mappingNewToInitialOrdering;

            constraintMatrix.Gauss(constraintMappingRHS, mappingNewToInitialOrdering, mappingInitialToNewOrdering);

            // move dependent dofs at the end
            // Warning!!! after this loop mappingNewToInitialOrdering is no longer valid !!!
            std::vector<int> tmpMapping;
            for (int dependentDofCount = 0; dependentDofCount < numDependentDofs; dependentDofCount++)
            {
                tmpMapping.push_back(numActiveDofs + dependentDofCount);
                mappingInitialToNewOrdering[mappingNewToInitialOrdering[dependentDofCount]] += numActiveDofs;
            }
            for (int activeDofCount = numDependentDofs; activeDofCount < numDofs; activeDofCount++)
            {
                tmpMapping.push_back(activeDofCount - numDependentDofs);
                mappingInitialToNewOrdering[mappingNewToInitialOrdering[activeDofCount]] -= numDependentDofs;
            }
            mappingNewToInitialOrdering.clear();

            // reorder columns
            constraintMatrix.ReorderColumns(tmpMapping);


            // remove columns of dependent dofs
            // check if the submatrix which is removed is a diagonal matrix
            const auto& columns = constraintMatrix.GetColumns();

            for (unsigned int iRow = 0; iRow < columns.size(); iRow++)
                for (unsigned int iPos = 0; iPos < columns[iRow].size(); iPos++)
                {
                    int column = columns[iRow][iPos];
                    if (column > numActiveDofs)
                        if (column - numActiveDofs != (int)iRow)
                            throw MechanicsException("[NuTo::Structure::NodeBuildGlobalDofs] invalid matrix structure.");
                }


            constraintMatrix.RemoveLastColumns(numDependentDofs);

            // renumber dofs

            BOOST_FOREACH(auto it, mNodeMap)
            {
                it->second->RenumberGlobalDofs(dof, mappingInitialToNewOrdering);
            }
        }

        // since only the diagonals were set, the off-diagonal submatrices have to be resized
        // to guarantee the right dimensions in arithmetic operations
        mConstraintMatrix.FixOffDiagonalDimensions();
        mConstraintMappingRHS.FixOffDiagonalDimensions();

        mConstraintMatrix.CheckDimensions();
        mConstraintMappingRHS.CheckDimensions();


        // number Lagrange multipliers in constraint equations defined in StructureBase
        // currently removed
        //renumber DOFS in constraints (Lagrange multiplier)
        //        ConstraintRenumberGlobalDofs(mappingInitialToNewOrdering);

        //calculate current rhs matrix
        this->ConstraintUpdateRHSAfterGaussElimination();

        mNodeNumberingRequired = false;
        UpdateDofStatus(); // again

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
    if (this->mNodeNumberingRequired)
        throw MechanicsException(std::string("[") + __PRETTY_FUNCTION__ +"] a valid dof numbering was not found (build dof numbering using NodeBuildGlobalDofs).");

    StructureOutputBlockVector dofValues(GetDofStatus(), true); // with resize

    for (auto dofType : DofTypesGet())
    {
        auto& actDofValues = dofValues.J[dofType];
        auto& depDofValues = dofValues.K[dofType];

        // extract dof values from nodes
        BOOST_FOREACH(auto it, mNodeMap)
            it->second->GetGlobalDofValues(rTimeDerivative, dofType, actDofValues, depDofValues);

        //extract dof values of additional DOFs
        // NodeExtractAdditionalGlobalDofValues(rTimeDerivative, rActiveDofValues[dof],rDependentDofValues[dof]);
        // I do not get this method. But obviously, it is not implemented at the moment. - TT
    }
    return dofValues;
}

void NuTo::Structure::NodeMergeDofValues(int rTimeDerivative, const NuTo::BlockFullVector<double>& rActiveDofValues, const NuTo::BlockFullVector<double>& rDependentDofValues)
{
    if (this->mNodeNumberingRequired)
        throw MechanicsException(std::string("[") + __PRETTY_FUNCTION__ +"] a valid dof numbering was not found (build dof numbering using NodeBuildGlobalDofs).");

    for (auto dofType : DofTypesGetActive())
    {
        if (rActiveDofValues[dofType].rows() != GetNumActiveDofs(dofType))
            throw MechanicsException(std::string("[") + __PRETTY_FUNCTION__ +"] invalid dimension of active dof vector for " + Node::DofToString(dofType));

        if (rDependentDofValues[dofType].rows() != GetNumDependentDofs(dofType))
            throw MechanicsException(std::string("[") + __PRETTY_FUNCTION__ +"] invalid dimension of dependent dof vector for " + Node::DofToString(dofType));

        // write dof values to the nodes
        BOOST_FOREACH(auto it, mNodeMap)
            it->second->SetGlobalDofValues(rTimeDerivative, dofType, rActiveDofValues[dofType], rDependentDofValues[dofType]);


        //write the dof values of the Lagrange multipliers
//        ConstraintMergeGlobalDofValues(rActiveDofValues, rDependentDofValues);

        //write dof values of additional DOFs
//        NodeMergeAdditionalGlobalDofValues(rTimeDerivative,rActiveDofValues,rDependentDofValues);
        // I do not get this method. But obviously, it is not implemented at the moment. - TT
    }
    this->mUpdateTmpStaticDataRequired=true;
}




NuTo::BlockFullVector<double> NuTo::Structure::NodeCalculateDependentDofValues(const NuTo::BlockFullVector<double>& rActiveDofValues) const
{
    if (this->mNodeNumberingRequired)
        throw MechanicsException(std::string("[") + __PRETTY_FUNCTION__ +"] a valid dof numbering was not found (build dof numbering using NodeBuildGlobalDofs).");

    return this->mConstraintRHS - this->mConstraintMatrix * rActiveDofValues;
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
    for(boost::ptr_map<int,ConstraintBase>::iterator constraintIt=mConstraintMap.begin();constraintIt!=mConstraintMap.end(); ++constraintIt)
    {
        constraintIt->second->ExchangeNodePtr(rOldPtr, rNewPtr);
    }
}
