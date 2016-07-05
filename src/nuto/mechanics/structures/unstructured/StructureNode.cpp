#include <boost/tokenizer.hpp>
#include <sstream>

#include "nuto/base/Timer.h"

#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/mechanics/groups/Group.h"
#include "nuto/mechanics/nodes/NodeDof.h"

//! @brief returns the number of nodes
//! @return number of nodes
int NuTo::Structure::GetNumNodes() const
{
    return mNodeMap.size();
}

//! @brief a reference to a node
//! @param identifier
//! @return reference to a node
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

//! @brief a reference to a node
//! @param identifier
//! @return reference to a node
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

//! @brief gives the identifier of a node
//! @param reference to a node
//! @return identifier
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


//! @brief ... return the global dof number of the displacement component of a node
//! @param rNodeId (Input) 			... node id
//! @param rDispDof 	... local disp dof (0,1 or 2 for x,y or z)
//! @returnrglobal dof number
int NuTo::Structure::NodeGetDofDisplacement(int rNodeId, int rDispDof)
{
	NodeBase* nodePtr = NodeGetNodePtr(rNodeId);
	if (nodePtr->GetNumDisplacements()>rDispDof)
	{
		return nodePtr->GetDofDisplacement(rDispDof);
	}
	else
	   throw MechanicsException("[NuTo::Structure::NodeGetDisplacementDofs] Node does have sufficient disp dofs.");
}

//! @brief ... store all element ids connected to this node in a vector
//! @param rNode (Input) 			... node id
//! @param rElementNumbers (Output) ... vector of element ids
void NuTo::Structure::NodeGetElements(const int rNodeId, NuTo::FullVector<int,Eigen::Dynamic>& rElementNumbers)
{
	const NuTo::NodeBase* nodePtr(NuTo::Structure::NodeGetNodePtr(rNodeId));
	std::vector<NuTo::ElementBase*> elementPtrs;
	this->NodeGetElements(nodePtr,elementPtrs);
	rElementNumbers.Resize(elementPtrs.size());
	size_t i=0;
	BOOST_FOREACH(NuTo::ElementBase* thisElPtr,elementPtrs)
		rElementNumbers(++i,1)=thisElPtr->ElementGetId();
}

//! @brief ... store all elements connected to this node in a vector
//! @param rNode (Input) 		... node pointer
//! @param rElements (Output) 	... vector of element pointers
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


//! @brief info about the elements in the Structure
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
				mLogger << "\t coords:";
				for(unsigned short iDof=0; iDof<it->second->GetNumCoordinates(); ++iDof)
					mLogger << "\t" << it->second->GetCoordinate(iDof);
                if (it->second->GetNumDisplacements()>0)
                {
                	mLogger << "\t disp:";
                    for(unsigned short iDof=0; iDof<it->second->GetNumDisplacements(); ++iDof)
                    {
                    	mLogger << "\t" << it->second->GetDisplacement(iDof) ;
                       	mLogger << "("<< it->second->GetDofDisplacement(iDof)<< ")" ;
                    }
                }
                if (it->second->GetNumNonlocalEqStrain()>0)
                {
                	mLogger << "\t nonlocal eq strain:";
                    for(unsigned short iDof=0; iDof<it->second->GetNumNonlocalEqStrain(); ++iDof)
                    {
                    	mLogger << "\t" << it->second->GetNonlocalEqStrain(iDof) ;
                       	mLogger << "("<< it->second->GetDofNonlocalEqStrain()<< ")" ;
                    }
                }
                if (it->second->GetNumRotations()>0 )
                {
                	mLogger << "\t rotations:";
                    for(unsigned short iDof=0; iDof<it->second->GetNumRotations(); ++iDof)
                    {
                    	mLogger << "\t" << it->second->GetRotation(iDof) ;
                        mLogger << "("<< it->second->GetDofRotation(iDof)<< ")" ;
                    }
                }
                if (it->second->GetNumTemperature()>0)
                {
                    mLogger << "\t Temperature:";
                    for(unsigned short iDof=0; iDof < it->second->GetNumTemperature(); ++iDof)
                    {
                    	mLogger << "\t" << it->second->GetTemperature();
                        mLogger << "("<< it->second->GetDofTemperature() << ")";
                    }
                }
                if (it->second->GetNumNonlocalTotalStrain()>0 )
                {
                	mLogger << "\t ns:";
                    for(unsigned short iDof=0; iDof<it->second->GetNumNonlocalTotalStrain(); ++iDof)
                    {
                    	mLogger << "\t" << it->second->GetNonlocalTotalStrain(iDof) ;
                        mLogger << "("<< it->second->GetDofNonlocalTotalStrain(iDof)<< ")" ;
                    }
                }
                if (it->second->GetNumWaterVolumeFraction()>0 )
                {
                    mLogger << "\t water volume fraction:";
                    mLogger << "\t" << it->second->GetWaterVolumeFraction() ;
                    mLogger << "("<< it->second->GetDofWaterVolumeFraction()<< ")" ;
                }
                if (it->second->GetNumRelativeHumidity()>0 )
                {
                    mLogger << "\t relative humidity:";
                    mLogger << "\t" << it->second->GetRelativeHumidity() ;
                    mLogger << "("<< it->second->GetDofRelativeHumidity()<< ")" ;
                }
                if (it->second->GetNumDamage()>0 )
                {
                    mLogger << "\t damage:";
                    mLogger << "\t" << it->second->GetDamage() ;
                    mLogger << "("<< it->second->GetDofDamage()<< ")" ;
                }
			}
			mLogger << "\n";
    	}

    }
}

//! creates a node at coordinate's origin
int NuTo::Structure::NodeCreate()
{
	NuTo::FullVector<double,Eigen::Dynamic> coordinates(this->GetDimension());
	coordinates.setZero();

	//return int identifier of the new node
    return NodeCreate(coordinates);
}

//! creates a node
int NuTo::Structure::NodeCreate(NuTo::FullVector<double,Eigen::Dynamic> rCoordinates)
{

    //find unused integer id
    int id(mNodeMap.size());
    boost::ptr_map<int,NodeBase>::iterator it = mNodeMap.find(id);
    while (it!=mNodeMap.end())
    {
        id++;
        it = mNodeMap.find(id);
    }
    this->NodeCreate(id, rCoordinates);

    //return int identifier of the new node
    return id;
}

//! creates a node with rDofs degrees of freedom
int NuTo::Structure::NodeCreate(NuTo::FullVector<double,Eigen::Dynamic> rCoordinates, std::set<NuTo::Node::eDof> rDofs)
{

    //find unused integer id
    int id(mNodeMap.size());
    boost::ptr_map<int,NodeBase>::iterator it = mNodeMap.find(id);
    while (it!=mNodeMap.end())
    {
        id++;
        it = mNodeMap.find(id);
    }

    NodeBase* nodePtr = NodePtrCreate(rDofs, rCoordinates);

    // add node to map
    this->mNodeMap.insert(id, nodePtr);

    //renumbering of dofs for global matrices required
    this->mNodeNumberingRequired  = true;

    //return int identifier of the new node
    return id;
}




NuTo::NodeBase* NuTo::Structure::NodePtrCreate(std::set<Node::eDof> rDOFs, NuTo::FullVector<double, Eigen::Dynamic> rCoordinates)
{

    if (rCoordinates.GetNumRows() != mDimension || rCoordinates.GetNumColumns() != 1)
        throw MechanicsException("[NuTo::Structure::NodeCreate] Dimension of the coordinate vector does not fit the dimension of the structure");

    int attributes = 1 << Node::COORDINATES;
    for (Node::eDof dof : rDOFs)
    {
        attributes = attributes | 1 << dof;
    }

    NodeBase* nodePtr;
    switch (attributes)
    {
    // the << shifts the 1 bitwise to the left, so 1<<n = 2^n
    // it actually sets the n-th bit (from the right) to 1, and all the other to zero
    case (1 << Node::COORDINATES):
        // reference node having only coordinates
        switch (mDimension)
        {
        case 1:
            nodePtr = new NuTo::NodeDof<1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0>();
            break;
        case 2:
            nodePtr = new NuTo::NodeDof<2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0>();
            break;
        case 3:
            nodePtr = new NuTo::NodeDof<3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0>();
            break;
        default:
            throw MechanicsException("[NuTo::Structure::NodeCreate] Dimension of the structure is not valid.");
        }
        break;
    case (1 << Node::COORDINATES) | (1 << Node::DISPLACEMENTS):
        // coordinates and displacements
        switch (mNumTimeDerivatives)
        {
        case 0:
            switch (mDimension)
            {
            case 1:
                nodePtr = new NuTo::NodeDof<1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0>();
                break;
            case 2:
                nodePtr = new NuTo::NodeDof<2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0>();
                break;
            case 3:
                nodePtr = new NuTo::NodeDof<3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0>();
                break;
            default:
                throw MechanicsException("[NuTo::Structure::NodeCreate] Dimension of the structure is not valid.");
            }
            break;
        case 1:
            switch (mDimension)
            {
            case 1:
                nodePtr = new NuTo::NodeDof<1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0>();
                break;
            case 2:
                nodePtr = new NuTo::NodeDof<2, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0>();
                break;
            case 3:
                nodePtr = new NuTo::NodeDof<3, 1, 3, 0, 0, 0, 0, 0, 0, 0, 0>();
                break;
            default:
                throw MechanicsException("[NuTo::Structure::NodeCreate] Dimension of the structure is not valid.");
            }
            break;
        case 2:
            switch (mDimension)
            {
            case 1:
                nodePtr = new NuTo::NodeDof<1, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0>();
                break;
            case 2:
                nodePtr = new NuTo::NodeDof<2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0>();
                break;
            case 3:
                nodePtr = new NuTo::NodeDof<3, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0>();
                break;
            default:
                throw MechanicsException("[NuTo::Structure::NodeCreate] Dimension of the structure is not valid.");
            }
            break;
        default:
            throw MechanicsException("[NuTo::Structure::NodeCreate] Coordinates and Displacements only implemented for 0,1 and 2 time derivatives.");
        }
        break;
    case (1 << Node::COORDINATES) | (1 << Node::DISPLACEMENTS) | (1 << Node::ROTATIONS):
        // coordinates and displacements and rotations
        switch (mNumTimeDerivatives)
        {
        case 0:
            switch (mDimension)
            {
            case 1:
                nodePtr = new NuTo::NodeDof<1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0>();
                break;
            case 2:
                nodePtr = new NuTo::NodeDof<2, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0>();
                break;
            case 3:
                nodePtr = new NuTo::NodeDof<3, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0>();
                break;
            default:
                throw MechanicsException("[NuTo::Structure::NodeCreate] Dimension of the structure is not valid.");
            }
            break;
        case 2:
            switch (mDimension)
            {
            case 1:
                nodePtr = new NuTo::NodeDof<1, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0>();
                break;
            case 2:
                nodePtr = new NuTo::NodeDof<2, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0>();
                break;
            case 3:
                nodePtr = new NuTo::NodeDof<3, 2, 3, 3, 0, 0, 0, 0, 0, 0, 0>();
                break;
            default:
                throw MechanicsException("[NuTo::Structure::NodeCreate] Dimension of the structure is not valid.");
            }
            break;
        default:
            throw MechanicsException("[NuTo::Structure::NodeCreate] Coordinates, Displacements and Rotations only implemented for 0 and 2 time derivatives.");
        }
        break;
    case (1 << Node::COORDINATES) | (1 << Node::TEMPERATURE):
        // coordinates and temperatures
        switch (mNumTimeDerivatives)
        {
        case 0:
            switch (mDimension)
            {
            case 1:
                nodePtr = new NuTo::NodeDof<1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0>();
                break;
            case 2:
                nodePtr = new NuTo::NodeDof<2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0>();
                break;
            case 3:
                nodePtr = new NuTo::NodeDof<3, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0>();
                break;
            default:
                throw MechanicsException("[NuTo::Structure::NodeCreate] Dimension of the structure is not valid.");
            }
            break;
        case 1:
            switch (mDimension)
            {
            case 1:
                nodePtr = new NuTo::NodeDof<1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0>();
                break;
            case 2:
                nodePtr = new NuTo::NodeDof<2, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0>();
                break;
            case 3:
                nodePtr = new NuTo::NodeDof<3, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0>();
                break;
            default:
                throw MechanicsException("[NuTo::Structure::NodeCreate] Dimension of the structure is not valid.");
            }
            break;
        case 2:
            switch (mDimension)
            {
            case 1:
                nodePtr = new NuTo::NodeDof<1, 2, 0, 0, 1, 0, 0, 0, 0, 0, 0>();
                break;
            case 2:
                nodePtr = new NuTo::NodeDof<2, 2, 0, 0, 1, 0, 0, 0, 0, 0, 0>();
                break;
            case 3:
                nodePtr = new NuTo::NodeDof<3, 2, 0, 0, 1, 0, 0, 0, 0, 0, 0>();
                break;
            default:
                throw MechanicsException("[NuTo::Structure::NodeCreate] Dimension of the structure is not valid.");
            }
            break;
        default:
            throw MechanicsException("[NuTo::Structure::NodeCreate] Temperatures only implemented for 0, 1 and 2 time derivatives.");
        }
        break;
    case (1 << Node::COORDINATES) | (1 << Node::TEMPERATURE) | (1 << Node::DISPLACEMENTS):
        // coordinates and temperatures
        switch (mNumTimeDerivatives)
        {
        case 0:
            switch (mDimension)
            {
            case 1:
                nodePtr = new NuTo::NodeDof<1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0>();
                break;
            case 2:
                nodePtr = new NuTo::NodeDof<2, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0>();
                break;
            case 3:
                nodePtr = new NuTo::NodeDof<3, 0, 3, 0, 1, 0, 0, 0, 0, 0, 0>();
                break;
            default:
                throw MechanicsException("[NuTo::Structure::NodeCreate] Dimension of the structure is not valid.");
            }
            break;
        case 1:
            switch (mDimension)
            {
            case 1:
                nodePtr = new NuTo::NodeDof<1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0>();
                break;
            case 2:
                nodePtr = new NuTo::NodeDof<2, 1, 2, 0, 1, 0, 0, 0, 0, 0, 0>();
                break;
            case 3:
                nodePtr = new NuTo::NodeDof<3, 1, 3, 0, 1, 0, 0, 0, 0, 0, 0>();
                break;
            default:
                throw MechanicsException("[NuTo::Structure::NodeCreate] Dimension of the structure is not valid.");
            }
            break;
        case 2:
            switch (mDimension)
            {
            case 1:
                nodePtr = new NuTo::NodeDof<1, 2, 1, 0, 1, 0, 0, 0, 0, 0, 0>();
                break;
            case 2:
                nodePtr = new NuTo::NodeDof<2, 2, 2, 0, 1, 0, 0, 0, 0, 0, 0>();
                break;
            case 3:
                nodePtr = new NuTo::NodeDof<3, 2, 3, 0, 1, 0, 0, 0, 0, 0, 0>();
                break;
            default:
                throw MechanicsException("[NuTo::Structure::NodeCreate] Dimension of the structure is not valid.");
            }
            break;
        default:
            throw MechanicsException("[NuTo::Structure::NodeCreate] Temperatures only implemented for 0, 1 and 2 time derivatives.");
        }
        break;
    case (1 << Node::COORDINATES) | (1 << Node::DISPLACEMENTS) | (1 << Node::NONLOCALDATA):
        // coordinates, displacements and nonlocal data
        switch (mDimension)
        {
        case 1:
            throw MechanicsException("[NuTo::Structure::NodeCreate] A 1D node with nonlocal data is not implemented");
            break;
        case 2:
            throw MechanicsException("[NuTo::Structure::NodeCreate] nonlocal not yet implemented.");
            //nodePtr = new NuTo::NodeCoordinatesDisplacementsNonlocalData2D();
            break;
            //case 3:
            //	nodePtr = new NuTo::NodeCoordinatesDisplacementsNonlocalData3D();
            //	break;
        default:
            throw MechanicsException("[NuTo::Structure::NodeCreate] Dimension of the structure is not valid.");
        }
        break;
    case (1 << Node::COORDINATES) | (1 << Node::DISPLACEMENTS) | (1 << Node::NONLOCALEQPLASTICSTRAIN):
        // coordinates and displacements and nonlocal eq. plastic strains
        switch (mNumTimeDerivatives)
        {
        case 0:
            switch (mDimension)
            {
            case 1:
                nodePtr = new NuTo::NodeDof<1, 0, 1, 0, 0, 2, 0, 0, 0, 0, 0>();
                break;
            case 2:
                nodePtr = new NuTo::NodeDof<2, 0, 2, 0, 0, 2, 0, 0, 0, 0, 0>();
                break;
            case 3:
                nodePtr = new NuTo::NodeDof<3, 0, 3, 0, 0, 2, 0, 0, 0, 0, 0>();
                break;
            default:
                throw MechanicsException("[NuTo::Structure::NodeCreate] Dimension of the structure is not valid.");
            }
            break;
        case 2:
            switch (mDimension)
            {
            case 1:
                nodePtr = new NuTo::NodeDof<1, 2, 1, 0, 0, 2, 0, 0, 0, 0, 0>();
                break;
            case 2:
                nodePtr = new NuTo::NodeDof<2, 2, 2, 0, 0, 2, 0, 0, 0, 0, 0>();
                break;
            case 3:
                nodePtr = new NuTo::NodeDof<3, 2, 3, 0, 0, 2, 0, 0, 0, 0, 0>();
                break;
            default:
                throw MechanicsException("[NuTo::Structure::NodeCreate] Dimension of the structure is not valid.");
            }
            break;
        default:
            throw MechanicsException("[NuTo::Structure::NodeCreate] Coordinates, Displacements and nonlocal eq. plastic strains only implemented for 0 and 2 time derivatives.");
        }
        break;
    case (1 << Node::COORDINATES) | (1 << Node::DISPLACEMENTS) | (1 << Node::NONLOCALTOTALSTRAIN):
        // coordinates and displacements and nonlocal total strains
        switch (mNumTimeDerivatives)
        {
        case 0:
            switch (mDimension)
            {
            case 1:
                nodePtr = new NuTo::NodeDof<1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0>();
                break;
            case 2:
                nodePtr = new NuTo::NodeDof<2, 0, 2, 0, 0, 0, 3, 0, 0, 0, 0>();
                break;
            case 3:
                nodePtr = new NuTo::NodeDof<3, 0, 3, 0, 0, 0, 6, 0, 0, 0, 0>();
                break;
            default:
                throw MechanicsException("[NuTo::Structure::NodeCreate] Dimension of the structure is not valid.");
            }
            break;
        case 2:
            switch (mDimension)
            {
            case 1:
                nodePtr = new NuTo::NodeDof<1, 2, 1, 0, 0, 0, 1, 0, 0, 0, 0>();
                break;
            case 2:
                nodePtr = new NuTo::NodeDof<2, 2, 2, 0, 0, 0, 3, 0, 0, 0, 0>();
                break;
            case 3:
                nodePtr = new NuTo::NodeDof<3, 2, 3, 0, 0, 0, 6, 0, 0, 0, 0>();
                break;
            default:
                throw MechanicsException("[NuTo::Structure::NodeCreate] Dimension of the structure is not valid.");
            }
            break;
        default:
            throw MechanicsException("[NuTo::Structure::NodeCreate] Coordinates, Displacements and onlocal total strains only implemented for 0 and 2 time derivatives.");
        }
        break;
    case (1 << Node::COORDINATES) | (1 << Node::DISPLACEMENTS) | (1 << Node::NONLOCALEQSTRAIN):
        // coordinates and displacements and nonlocal eq. strains
        switch (mNumTimeDerivatives)
        {
        case 0:
            switch (mDimension)
            {
            case 1:
                nodePtr = new NuTo::NodeDof<1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0>();
                break;
            case 2:
                nodePtr = new NuTo::NodeDof<2, 0, 2, 0, 0, 0, 0, 1, 0, 0, 0>();
                break;
            case 3:
                nodePtr = new NuTo::NodeDof<3, 0, 3, 0, 0, 0, 0, 1, 0, 0, 0>();
                break;
            default:
                throw MechanicsException("[NuTo::Structure::NodeCreate] Dimension of the structure is not valid.");
            }
            break;
        case 2:

            switch (mDimension)
            {
            case 1:
                nodePtr = new NuTo::NodeDof<1, 2, 1, 0, 0, 0, 0, 1, 0, 0, 0>();
                break;
            case 2:
                nodePtr = new NuTo::NodeDof<2, 2, 2, 0, 0, 0, 0, 1, 0, 0, 0>();
                break;
            case 3:
                nodePtr = new NuTo::NodeDof<3, 2, 3, 0, 0, 0, 0, 1, 0, 0, 0>();
                break;
            default:
                throw MechanicsException("[NuTo::Structure::NodeCreate] Dimension of the structure is not valid.");
            }
            break;
        default:
            throw MechanicsException("[NuTo::Structure::NodeCreate] Coordinates, Displacements and nonlocal eq strains only implemented for 0 time derivatives.");
        }
        break;
    case (1 << Node::COORDINATES) | (1 << Node::NONLOCALEQSTRAIN):
        // coordinates and displacements and nonlocal eq. strains
        switch (mNumTimeDerivatives)
        {
        case 0:
            switch (mDimension)
            {
            case 1:
                nodePtr = new NuTo::NodeDof<1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0>();
                break;
            case 2:
                nodePtr = new NuTo::NodeDof<2, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0>();
                break;
            case 3:
                nodePtr = new NuTo::NodeDof<3, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0>();
                break;
            default:
                throw MechanicsException("[NuTo::Structure::NodeCreate] Dimension of the structure is not valid.");
            }
            break;
            //          case 2:
            //
            //              switch (mDimension)
            //              {
            //              case 1:
            //                  nodePtr = new NuTo::NodeDof<1,2,0,0,0,0,0,1,0,0, 0>();
            //                  break;
            //              case 2:
            //                  nodePtr = new NuTo::NodeDof<2,2,0,0,0,0,0,1,0,0, 0>();
            //                  break;
            //              case 3:
            //                  nodePtr = new NuTo::NodeDof<3,2,0,0,0,0,0,1,0,0, 0>();
            //                  break;
            //              default:
            //                  throw MechanicsException("[NuTo::Structure::NodeCreate] Dimension of the structure is not valid.");
            //              }
            //          break;
        default:
            throw MechanicsException("[NuTo::Structure::NodeCreate] Coordinates, Displacements and nonlocal eq strains only implemented for 0 time derivatives.");
        }
        break;
        // Moisture Transport --- Beginn
    case (1 << Node::COORDINATES) | (1 << Node::DISPLACEMENTS) | (1 << Node::WATERVOLUMEFRACTION) | (1 << Node::RELATIVEHUMIDITY):
        // coordinates, displacements, water volume fraction and relative humidity
        switch (mNumTimeDerivatives)
        {
        case 0:
            switch (mDimension)
            {
            case 1:
                nodePtr = new NuTo::NodeDof<1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0>();
                break;
            case 2:
                nodePtr = new NuTo::NodeDof<2, 0, 2, 0, 0, 0, 0, 0, 1, 1, 0>();
                break;
            case 3:
                nodePtr = new NuTo::NodeDof<3, 0, 3, 0, 0, 0, 0, 0, 1, 1, 0>();
                break;
            default:
                throw MechanicsException("[NuTo::Structure::NodeCreate] Dimension of the structure is not valid.");
            }
            break;
        case 1:
            switch (mDimension)
            {
            case 1:
                nodePtr = new NuTo::NodeDof<1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0>();
                break;
            case 2:
                nodePtr = new NuTo::NodeDof<2, 1, 2, 0, 0, 0, 0, 0, 1, 1, 0>();
                break;
            case 3:
                nodePtr = new NuTo::NodeDof<3, 1, 3, 0, 0, 0, 0, 0, 1, 1, 0>();
                break;
            default:
                throw MechanicsException("[NuTo::Structure::NodeCreate] Dimension of the structure is not valid.");
            }
            break;
        case 2:
            switch (mDimension)
            {
            case 1:
                nodePtr = new NuTo::NodeDof<1, 2, 1, 0, 0, 0, 0, 0, 1, 1, 0>();
                break;
            case 2:
                nodePtr = new NuTo::NodeDof<2, 2, 2, 0, 0, 0, 0, 0, 1, 1, 0>();
                break;
            case 3:
                nodePtr = new NuTo::NodeDof<3, 2, 3, 0, 0, 0, 0, 0, 1, 1, 0>();
                break;
            default:
                throw MechanicsException("[NuTo::Structure::NodeCreate] Dimension of the structure is not valid.");
            }
            break;
        default:
            throw MechanicsException("[NuTo::Structure::NodeCreate] Coordinates, Displacements, Relative Humidity and Water Volume Fraction only implemented for 0 and 2 time derivatives.");
        }
        break;
    case (1 << Node::COORDINATES) | (1 << Node::WATERVOLUMEFRACTION) | (1 << Node::RELATIVEHUMIDITY):
        // coordinates, water volume fraction and relative humidity
        switch (mNumTimeDerivatives)
        {
        case 0:
            switch (mDimension)
            {
            case 1:
                nodePtr = new NuTo::NodeDof<1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0>();
                break;
            case 2:
                nodePtr = new NuTo::NodeDof<2, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0>();
                break;
            case 3:
                nodePtr = new NuTo::NodeDof<3, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0>();
                break;
            default:
                throw MechanicsException("[NuTo::Structure::NodeCreate] Dimension of the structure is not valid.");
            }
            break;
        case 1:
            switch (mDimension)
            {
            case 1:
                nodePtr = new NuTo::NodeDof<1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0>();
                break;
            case 2:
                nodePtr = new NuTo::NodeDof<2, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0>();
                break;
            case 3:
                nodePtr = new NuTo::NodeDof<3, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0>();
                break;
            default:
                throw MechanicsException("[NuTo::Structure::NodeCreate] Dimension of the structure is not valid.");
            }
            break;
        case 2:
            switch (mDimension)
            {
            case 1:
                nodePtr = new NuTo::NodeDof<1, 2, 0, 0, 0, 0, 0, 0, 1, 1, 0>();
                break;
            case 2:
                nodePtr = new NuTo::NodeDof<2, 2, 0, 0, 0, 0, 0, 0, 1, 1, 0>();
                break;
            case 3:
                nodePtr = new NuTo::NodeDof<3, 2, 0, 0, 0, 0, 0, 0, 1, 1, 0>();
                break;
            default:
                throw MechanicsException("[NuTo::Structure::NodeCreate] Dimension of the structure is not valid.");
            }
            break;
        default:
            throw MechanicsException("[NuTo::Structure::NodeCreate] Coordinates, Relative Humidity and Water Volume Fraction only implemented for 0 and 2 time derivatives.");
        }
        break;
    case (1 << Node::COORDINATES) | (1 << Node::RELATIVEHUMIDITY):
        // coordinates and relative humidity
        switch (mNumTimeDerivatives)
        {
        case 0:
            switch (mDimension)
            {
            case 1:
                nodePtr = new NuTo::NodeDof<1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0>();
                break;
            case 2:
                nodePtr = new NuTo::NodeDof<2, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0>();
                break;
            case 3:
                nodePtr = new NuTo::NodeDof<3, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0>();
                break;
            default:
                throw MechanicsException("[NuTo::Structure::NodeCreate] Dimension of the structure is not valid.");
            }
            break;
        case 1:
            switch (mDimension)
            {
            case 1:
                nodePtr = new NuTo::NodeDof<1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0>();
                break;
            case 2:
                nodePtr = new NuTo::NodeDof<2, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0>();
                break;
            case 3:
                nodePtr = new NuTo::NodeDof<3, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0>();
                break;
            default:
                throw MechanicsException("[NuTo::Structure::NodeCreate] Dimension of the structure is not valid.");
            }
            break;
        case 2:
            switch (mDimension)
            {
            case 1:
                nodePtr = new NuTo::NodeDof<1, 2, 0, 0, 0, 0, 0, 0, 0, 1, 0>();
                break;
            case 2:
                nodePtr = new NuTo::NodeDof<2, 2, 0, 0, 0, 0, 0, 0, 0, 1, 0>();
                break;
            case 3:
                nodePtr = new NuTo::NodeDof<3, 2, 0, 0, 0, 0, 0, 0, 0, 1, 0>();
                break;
            default:
                throw MechanicsException("[NuTo::Structure::NodeCreate] Dimension of the structure is not valid.");
            }
            break;
        default:
            throw MechanicsException("[NuTo::Structure::NodeCreate] Coordinates and Relative Humidity only implemented for 0 and 2 time derivatives.");
        }
        break;
    case (1 << Node::COORDINATES) | (1 << Node::WATERVOLUMEFRACTION):
        // coordinates and water volume fraction
        switch (mNumTimeDerivatives)
        {
        case 0:
            switch (mDimension)
            {
            case 1:
                nodePtr = new NuTo::NodeDof<1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0>();
                break;
            case 2:
                nodePtr = new NuTo::NodeDof<2, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0>();
                break;
            case 3:
                nodePtr = new NuTo::NodeDof<3, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0>();
                break;
            default:
                throw MechanicsException("[NuTo::Structure::NodeCreate] Dimension of the structure is not valid.");
            }
            break;
        case 1:
            switch (mDimension)
            {
            case 1:
                nodePtr = new NuTo::NodeDof<1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0>();
                break;
            case 2:
                nodePtr = new NuTo::NodeDof<2, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0>();
                break;
            case 3:
                nodePtr = new NuTo::NodeDof<3, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0>();
                break;
            default:
                throw MechanicsException("[NuTo::Structure::NodeCreate] Dimension of the structure is not valid.");
            }
            break;
        case 2:
            switch (mDimension)
            {
            case 1:
                nodePtr = new NuTo::NodeDof<1, 2, 0, 0, 0, 0, 0, 0, 1, 0, 0>();
                break;
            case 2:
                nodePtr = new NuTo::NodeDof<2, 2, 0, 0, 0, 0, 0, 0, 1, 0, 0>();
                break;
            case 3:
                nodePtr = new NuTo::NodeDof<3, 2, 0, 0, 0, 0, 0, 0, 1, 0, 0>();
                break;
            default:
                throw MechanicsException("[NuTo::Structure::NodeCreate] Dimension of the structure is not valid.");
            }
            break;
        default:
            throw MechanicsException("[NuTo::Structure::NodeCreate] Coordinates and Water Volume Fraction only implemented for 0 and 2 time derivatives.");
        }
        break;
        // Moisture Transport --- Ende
    case (1 << Node::COORDINATES) | (1 << Node::DISPLACEMENTS) | (1 << Node::DAMAGE):
        // coordinates and displacements and damage
        switch (mNumTimeDerivatives)
        {
        case 0:
            switch (mDimension)
            {
            case 1:
                nodePtr = new NuTo::NodeDof<1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1>();
                break;
            case 2:
                nodePtr = new NuTo::NodeDof<2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 1>();
                break;
            case 3:
                nodePtr = new NuTo::NodeDof<3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 1>();
                break;
            default:
                throw MechanicsException("[NuTo::Structure::NodeCreate] Dimension of the structure is not valid.");
            }
            break;
        case 1:

            switch (mDimension)
            {
            case 1:
                nodePtr = new NuTo::NodeDof<1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1>();
                break;
            case 2:
                nodePtr = new NuTo::NodeDof<2, 1, 2, 0, 0, 0, 0, 0, 0, 0, 1>();
                break;
            case 3:
                nodePtr = new NuTo::NodeDof<3, 1, 3, 0, 0, 0, 0, 0, 0, 0, 1>();
                break;
            default:
                throw MechanicsException("[NuTo::Structure::NodeCreate] Dimension of the structure is not valid.");
            }
            break;
        default:
            throw MechanicsException("[NuTo::Structure::NodeCreate] Coordinates, Displacements and damage only implemented for 0,1 time derivatives.");
        }
        break;
    default:
        throw MechanicsException("[NuTo::Structure::NodeCreate] This combination of attributes is not implemented (just add in the source file the relevant combination).");
    }
    //add coordinates
    try
    {
        switch (mDimension)
        {
        case 1:
            nodePtr->SetCoordinates1D(rCoordinates);
            break;
        case 2:
            nodePtr->SetCoordinates2D(rCoordinates);
            break;
        case 3:
            nodePtr->SetCoordinates3D(rCoordinates);
            break;
        case 0:
            throw MechanicsException("[NuTo::StructureBase::NodeCreate] Node has no displacements.");
            break;
        }
    } catch (std::bad_cast& b)
    {
        throw MechanicsException("[NuTo::StructureBase::NodeCreate] Node has no coordinates or its dimension is not equivalent to the dimension of the input matrix.");
    } catch (NuTo::MechanicsException& b)
    {
        b.AddMessage("[NuTo::StructureBase::NodeCreate] Error setting coordinates.");
        throw b;
    } catch (...)
    {
        throw MechanicsException("[NuTo::StructureBase::NodeCreate] Error setting coordinates of node (unspecified exception).");
    }
    return nodePtr;
}

void NuTo::Structure::NodeCreate(int rNodeNumber, NuTo::FullVector<double,Eigen::Dynamic> rCoordinates)
{
	// check node number
	boost::ptr_map<int,NodeBase>::iterator it = this->mNodeMap.find(rNodeNumber);
	if(it != this->mNodeMap.end())
	{
		throw MechanicsException("[NuTo::Structure::NodeCreate] Node already exists.");
	}

	std::set<Node::eDof> dofs;
	dofs.insert(Node::COORDINATES);

    NodeBase* nodePtr = NodePtrCreate(dofs, rCoordinates);

    // add node to map
    this->mNodeMap.insert(rNodeNumber, nodePtr);

    //renumbering of dofs for global matrices required
    this->mNodeNumberingRequired  = true;
}

//! @brief creates a node with specific dofs at coordinate's origin
//! @param rDOFs ... space separated string containing the node dofs (e.g. displacements, rotations, temperatures)
//! @param rCoordinates ...  node coordinates
//! @return node number
int NuTo::Structure::NodeCreateDOFs(std::string rDOFs)
{
    NuTo::FullVector<double,Eigen::Dynamic> coordinates(this->GetDimension());
    coordinates.setZero();

    //return int identifier of the new node
    return NodeCreateDOFs(rDOFs, coordinates);
}

int NuTo::Structure::NodeCreateDOFs(std::string rDOFs, NuTo::FullVector<double,Eigen::Dynamic> rCoordinates)
{
    //find unused integer id
    int id(mNodeMap.size());
    boost::ptr_map<int,NodeBase>::iterator it = mNodeMap.find(id);
    while (it!=mNodeMap.end())
    {
        id++;
        it = mNodeMap.find(id);
    }

    NodeCreateDOFs(id, rDOFs, rCoordinates);

    return id;
}

//! @brief creates a node with specific dofs
//! @param node number
//! @param rDOFs ... space separated string containing the node dofs (e.g. displacements, rotations, temperatures)
//! @param rCoordinates ...  node coordinates
void NuTo::Structure::NodeCreateDOFs(int rNodeNumber, std::string rDOFs, NuTo::FullVector<double,Eigen::Dynamic> rCoordinates)
{

    // check node number
    boost::ptr_map<int,NodeBase>::iterator it = this->mNodeMap.find(rNodeNumber);
    if(it != this->mNodeMap.end())
    {
        throw MechanicsException("[NuTo::Structure::NodeCreateDOFs] Node already exists.");
    }

    std::set<Node::eDof> dofs;
    dofs.insert(Node::COORDINATES);

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
            throw e;
        }
    }

    NodeBase* nodePtr = NodePtrCreate(dofs, rCoordinates);

    // add node to map
    this->mNodeMap.insert(rNodeNumber, nodePtr);

    //renumbering of dofs for global matrices required
    this->mNodeNumberingRequired  = true;

}


//! @brief creates a node with specific dofs at coordinate's origin
//! @param rDOFs ... set containing the node dof enums (e.g. displacements, rotations, temperatures)
//! @param rCoordinates ...  node coordinates
//! @return node number
int NuTo::Structure::NodeCreateDOFs(std::set<NuTo::Node::eDof> rDOFs)
{
    NuTo::FullVector<double,Eigen::Dynamic> coordinates(this->GetDimension());
    coordinates.setZero();

    //return int identifier of the new node
    return NodeCreateDOFs(rDOFs, coordinates);
}


//! @brief creates a node with specific dofs
//! @param rDOFs ... set containing the node dof enums (e.g. displacements, rotations, temperatures)
//! @param rCoordinates ...  node coordinates
//! @return node number
int NuTo::Structure::NodeCreateDOFs(std::set<NuTo::Node::eDof> rDOFs, NuTo::FullVector<double, Eigen::Dynamic> rCoordinates)
{
    //find unused integer id
    int id(mNodeMap.size());
    boost::ptr_map<int,NodeBase>::iterator it = mNodeMap.find(id);
    while (it!=mNodeMap.end())
    {
        id++;
        it = mNodeMap.find(id);
    }

    NodeCreateDOFs(id, rDOFs, rCoordinates);

    return id;
}


//! @brief creates a node with specific dofs
//! @param node number
//! @param rDOFs ... set containing the node dof enums (e.g. displacements, rotations, temperatures)
//! @param rCoordinates ...  node coordinates
void NuTo::Structure::NodeCreateDOFs(int rNodeNumber, std::set<NuTo::Node::eDof> rDOFs, NuTo::FullVector<double, Eigen::Dynamic> rCoordinates)
{
    // check node number
    boost::ptr_map<int,NodeBase>::iterator it = this->mNodeMap.find(rNodeNumber);
    if(it != this->mNodeMap.end())
    {
        throw MechanicsException("[NuTo::Structure::NodeCreateDOFs] Node already exists.");
    }

    rDOFs.insert(Node::COORDINATES);

    NodeBase* nodePtr = NodePtrCreate(rDOFs, rCoordinates);

    // add node to map
    this->mNodeMap.insert(rNodeNumber, nodePtr);

    //renumbering of dofs for global matrices required
    this->mNodeNumberingRequired  = true;
}


//! @brief create multiple nodes
//! @param reference to a FullMatrix containing the node coordinates (row->coordinate; col->nodes)
//! @return a FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic> with the identifiers
NuTo::FullVector<int,Eigen::Dynamic> NuTo::Structure::NodesCreate(NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rCoordinates)
{
	std::vector<int> idVec;
	/// go through the nodes
	for(size_t i=0 ; i<(size_t)rCoordinates.GetNumColumns(); ++i)
	{
		NuTo::FullVector<double,Eigen::Dynamic> coordinate(rCoordinates.GetColumn(i));
		idVec.push_back(NodeCreate(coordinate));
	}

    //return int identifiers of the new nodes as FullMatrix
	NuTo::FullVector<int,Eigen::Dynamic> ids(idVec);

    return ids;
}

void NuTo::Structure::NodeDelete(int rNodeNumber)
{
	NodeDelete(rNodeNumber,true);
}

//! @brief Deletes a nodes
//! @param rNodeNumber identifier for the node
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

    	//! Search for node in groups
        //! using a loop over all groups
        //! remove the entry from all groups
        for(boost::ptr_map<int,GroupBase>::iterator groupIt=mGroupMap.begin();groupIt!=mGroupMap.end(); ++groupIt){
        	if(groupIt->second->GetType()==NuTo::Groups::Nodes){
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


//! @brief number the dofs in the structure
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

        throw e;
    }
    catch (MechanicsException& e)
    {
        e.AddMessage(std::string("[") + __PRETTY_FUNCTION__ + "] Error building global dof numbering. \n");
        if (not rCallerName.empty())
            e.AddMessage(std::string(" -- was called by [") + rCallerName + "]");

        throw e;
    }
}




//! @brief extract dof values (e.g. displacements, temperatures to the nodes)
//! @param rTimeDerivative time derivative (0 disp 1 vel 2 acc)
//! @return ... StructureBlockVector containing the dofs (J and K)
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

//! @brief write dof values (e.g. displacements, temperatures to the nodes)
//! @param rTimeDerivative time derivative (0 disp 1 vel 2 acc)
//! @param rActiveDofValues ... vector of independent dof values (ordering according to global dofs, size is number of active dofs)
//! @param rDependentDofValues ... vector of dependent  dof values (ordering according to global dofs, size is number of active dofs)
void NuTo::Structure::NodeMergeDofValues(int rTimeDerivative, const NuTo::BlockFullVector<double>& rActiveDofValues, const NuTo::BlockFullVector<double>& rDependentDofValues)
{
    if (this->mNodeNumberingRequired)
        throw MechanicsException(std::string("[") + __PRETTY_FUNCTION__ +"] a valid dof numbering was not found (build dof numbering using NodeBuildGlobalDofs).");

    for (auto dofType : DofTypesGetActive())
    {
        if (rActiveDofValues[dofType].GetNumRows() != GetNumActiveDofs(dofType))
            throw MechanicsException(std::string("[") + __PRETTY_FUNCTION__ +"] invalid dimension of active dof vector for " + Node::DofToString(dofType));

        if (rDependentDofValues[dofType].GetNumRows() != GetNumDependentDofs(dofType))
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




//! @brief calculate dependent dof values (for the zeroth time derivative)
//! @param rActiveDofValues ... vector of independent dof values (ordering according to global dofs, size is number of active dofs)
//! @return  ... vector of dependent  dof values (ordering according to global dofs, size is number of active dofs)
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

//! @brief exchanges the node ptr in the full data set (elements, groups, loads, constraints etc.)
//! this routine is used, if e.g. the data type of a node has changed, but the restraints, elements etc. are still identical
//! @param rId ... old node id
//! @param rOldPtr ... old node ptr
//! @param rNewPtr ... new node ptr
//! @param rElements (optional) ... vector of all elements that contain the node - speedup!
void NuTo::Structure::NodeExchangePtr(int rId, NuTo::NodeBase* rOldPtr, NuTo::NodeBase* rNewPtr, std::vector<ElementBase*> rElements)
{
    //in node map
    //find it
    if (mNodeMap.erase(rId)!=1)
    {
        throw MechanicsException("[NuTo::Structure::NodeExchangePtr] Pointer to node (to exchange) does not exist.");
    }

    mNodeMap.insert(rId,rNewPtr);

    if (rElements.size() == 0)
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
        if(groupIt->second->GetType()==NuTo::Groups::Nodes)
        {
        	if (groupIt->second->Contain(rId))
                groupIt->second->ExchangePtr(rId, rOldPtr, rNewPtr);
        }
    }

    //in constraints
    //in groups
    for(boost::ptr_map<int,ConstraintBase>::iterator constraintIt=mConstraintMap.begin();constraintIt!=mConstraintMap.end(); ++constraintIt)
    {
        constraintIt->second->ExchangeNodePtr(rOldPtr, rNewPtr);
    }


}

