#include "nuto/mechanics/elements/ContinuumElementIGA.h"
#include "nuto/mechanics/nodes/NodeBase.h"

#include "nuto/mechanics/sections/SectionTruss.h"
#include "nuto/mechanics/sections/SectionPlane.h"
#include "nuto/mechanics/elements/ElementOutputBase.h"
#include "nuto/mechanics/elements/ElementOutputIpData.h"
#include "nuto/mechanics/elements/ElementDataBase.h"

#include "nuto/mechanics/interpolationtypes/InterpolationType.h"

#include "nuto/mechanics/elements/EvaluateDataContinuum.h"

#include "nuto/mechanics/dofSubMatrixStorage/BlockFullVector.h"
#include "nuto/mechanics/dofSubMatrixStorage/BlockFullMatrix.h"

#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOBase.h"

template<int TDim>
NuTo::ContinuumElementIGA<TDim>::ContinuumElementIGA(const NuTo::StructureBase*          rStructure,
                                                     const std::vector<NuTo::NodeBase*> &rNodes,
                                                     const Eigen::MatrixXd              &rKnots,
                                                     const Eigen::VectorXi              &rKnotIDs,
                                                     ElementData::eElementDataType       rElementDataType,
                                                     IpData::eIpDataType                 rIpDataType,
                                                     InterpolationType                  *rInterpolationType)
    :
        ContinuumElement<TDim>(rStructure, rNodes, rElementDataType, rIpDataType, rInterpolationType),
        mKnots(rKnots),
        mKnotIDs(rKnotIDs)
{}

template<int TDim>
Eigen::Matrix<double, TDim, TDim> NuTo::ContinuumElementIGA<TDim>::CalculateJacobianParametricSpaceIGA() const
{
    Eigen::Matrix<double, TDim, TDim> jac;
    jac.setZero(TDim, TDim);
    for(int i = 0; i < TDim; i++) jac(i,i) = 0.5*(mKnots(i,1) - mKnots(i,0));

    return jac;
}

template<int TDim>
const Eigen::VectorXd NuTo::ContinuumElementIGA<TDim>::GetIntegrationPointVolume() const
{
    Eigen::MatrixXd nodeCoordinates = this->ExtractNodeValues(0, Node::COORDINATES);

    Eigen::VectorXd volume(this->GetNumIntegrationPoints());
    for (int theIP = 0; theIP < this->GetNumIntegrationPoints(); theIP++)
    {
        Eigen::MatrixXd derivativeShapeFunctionsNatural = this->mInterpolationType->Get(Node::COORDINATES).CalculateDerivativeShapeFunctionsNatural(theIP, mKnotIDs);
        double detJacobian = this->CalculateJacobian(derivativeShapeFunctionsNatural, nodeCoordinates).determinant();
        volume[theIP] = detJacobian * this->mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP);
    }
    return volume;
}

template<int TDim>
void NuTo::ContinuumElementIGA<TDim>::CheckElement()
{
    int numIntegrationPoints = this->GetNumIntegrationPoints();

    if (numIntegrationPoints < 1)
    {
        MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] invalid integration type.");
    }

    int theIP = 0;

    Eigen::MatrixXd derivativeShapeFunctions = this->mInterpolationType->Get(Node::COORDINATES).CalculateDerivativeShapeFunctionsNatural(theIP, mKnotIDs);
    Eigen::MatrixXd nodeCoordinates = this->ExtractNodeValues(0, Node::COORDINATES);
    double detJacobian = this->CalculateJacobian(derivativeShapeFunctions, nodeCoordinates).determinant();
    if (detJacobian < 0)
    {
        this->ReorderNodes();
        // recalculate node coordinates after reordering
        nodeCoordinates = this->ExtractNodeValues(0, Node::COORDINATES);
    }

    double size = 0;
    for (int iIP = 0; iIP < numIntegrationPoints; ++iIP)
    {
        Eigen::MatrixXd derivativeShapeFunctions = this->mInterpolationType->Get(Node::COORDINATES).CalculateDerivativeShapeFunctionsNatural(iIP, mKnotIDs);
        detJacobian = this->CalculateJacobian(derivativeShapeFunctions, nodeCoordinates).determinant();
        if (detJacobian <= 0)
        {
            throw MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] Determinant of the Jacobian <= zero, no inversion possible.");
        }
        size += this->GetIntegrationPointWeight(iIP) * detJacobian;
    }

    assert(std::abs(GetIntegrationPointVolume().sum() / size - 1) < 1.e-6); // just to be sure ...

    // check element length
    if (size < 1e-14)
    {
        MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] element with zero size (check nodes).");
    }

}

template<int TDim>
void NuTo::ContinuumElementIGA<TDim>::CalculateNMatrixBMatrixDetJacobian(EvaluateDataContinuum<TDim> &rData, int rTheIP) const
{
    // calculate Jacobian
    const Eigen::MatrixXd& derivativeShapeFunctionsGeometryNatural = this->mInterpolationType->Get(Node::COORDINATES).CalculateDerivativeShapeFunctionsNatural(rTheIP, mKnotIDs);

    Eigen::Matrix<double, TDim, TDim> jacobian = this->CalculateJacobian(derivativeShapeFunctionsGeometryNatural, rData.mNodalValues[Node::COORDINATES]);

    // there are two mappings in IGA: 1) reference element <=> parametric space (knots) 2) parametric space <=> physical space
    // the B-matrix only the invJac of mapping 2) is needed
    rData.mDetJacobian = jacobian.determinant();
    for(int i = 0; i < TDim; i++) rData.mDetJacobian *= 0.5*(mKnots(i,1) - mKnots(i,0));

    if (rData.mDetJacobian == 0)
        throw MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] Determinant of the Jacobian is zero, no inversion possible.");

    Eigen::Matrix<double, TDim, TDim> invJacobian = jacobian.inverse();

    // calculate shape functions and their derivatives
    for (auto dof : this->mInterpolationType->GetDofs())
    {
        if (dof == Node::COORDINATES)
            continue;
        const InterpolationBase& interpolationType = this->mInterpolationType->Get(dof);
        rData.mNIGA[dof] = interpolationType.CalculateMatrixN(rTheIP, mKnotIDs);

        rData.mB[dof] = this->CalculateMatrixB(dof, interpolationType.CalculateDerivativeShapeFunctionsNatural(rTheIP, mKnotIDs), invJacobian);
    }
}

template<int TDim>
NuTo::NodeBase* NuTo::ContinuumElementIGA<TDim>::GetNode(int rLocalNodeNumber)
{
    assert(rLocalNodeNumber >= 0);
    assert(rLocalNodeNumber < (int )((this->mNodes).size()));
    return (this->mNodes)[rLocalNodeNumber];
}

template<int TDim>
const NuTo::NodeBase* NuTo::ContinuumElementIGA<TDim>::GetNode(int rLocalNodeNumber) const
{
    assert(rLocalNodeNumber >= 0);
    assert(rLocalNodeNumber < (int )((this->mNodes).size()));
    return (this->mNodes)[rLocalNodeNumber];
}

template<int TDim>
void NuTo::ContinuumElementIGA<TDim>::SetNode(int rLocalNodeNumber, NodeBase* rNode)
{
    assert(rLocalNodeNumber >= 0);
    assert(rLocalNodeNumber < (int )((this->mNodes).size()));
    assert(rNode != nullptr);
    (this->mNodes)[rLocalNodeNumber] = rNode;
}

template<int TDim>
Eigen::VectorXd NuTo::ContinuumElementIGA<TDim>::InterpolateDofGlobal(int rTimeDerivative, const Eigen::VectorXd& rNaturalCoordinates, Node::eDof rDofType) const
{

    const InterpolationBase& interpolationType = this->mInterpolationType->Get(rDofType);
    Eigen::MatrixXd nodalValues = this->ExtractNodeValues(rTimeDerivative, rDofType);
    Eigen::MatrixXd matrixN = interpolationType.CalculateMatrixN(rNaturalCoordinates, mKnotIDs);

//    std::cout << matrixN * nodalValues << std::endl;

    return matrixN * nodalValues;
}

namespace NuTo // template specialization in *.cpp somehow requires the definition to be in the namespace...
{

template<>
double NuTo::ContinuumElementIGA<1>::CalculateDetJxWeightIPxSection(double rDetJacobian, int rTheIP) const
{
    Eigen::MatrixXd matrixN = mInterpolationType->Get(Node::COORDINATES).CalculateMatrixN(rTheIP, mKnotIDs);
    Eigen::VectorXd globalIPCoordinate = matrixN * ExtractNodeValues(0, Node::COORDINATES);

    return rDetJacobian * mElementData->GetIntegrationType()->GetIntegrationPointWeight(rTheIP) * mSection->GetArea() * mSection->AsSectionTruss()->GetAreaFactor(globalIPCoordinate.at(0, 0));
}

template<>
double NuTo::ContinuumElementIGA<2>::CalculateDetJxWeightIPxSection(double rDetJacobian, int rTheIP) const
{
    return rDetJacobian * mElementData->GetIntegrationType()->GetIntegrationPointWeight(rTheIP) * mSection->GetThickness();
}

template<>
double NuTo::ContinuumElementIGA<3>::CalculateDetJxWeightIPxSection(double rDetJacobian, int rTheIP) const
{
    return rDetJacobian * mElementData->GetIntegrationType()->GetIntegrationPointWeight(rTheIP);
}

}  // namespace NuTo

template class NuTo::ContinuumElementIGA<1>;
template class NuTo::ContinuumElementIGA<2>;
template class NuTo::ContinuumElementIGA<3>;


#ifdef ENABLE_SERIALIZATION
template void NuTo::ContinuumElementIGA<1>::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ContinuumElementIGA<2>::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ContinuumElementIGA<3>::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ContinuumElementIGA<1>::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ContinuumElementIGA<2>::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ContinuumElementIGA<3>::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ContinuumElementIGA<1>::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ContinuumElementIGA<2>::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ContinuumElementIGA<3>::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template<int TDim>
template<class Archive>
void NuTo::ContinuumElementIGA<TDim>::save(Archive & ar, const unsigned int version)const
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ContinuumElementIGA " << std::endl;
#endif
    ar & boost::serialization::make_nvp("ContinuumElementIGA_ElementBase",boost::serialization::base_object<ElementBase >(*this));
    ar & boost::serialization::make_nvp("mSection", const_cast<SectionBase*&>(mSection));

    const std::uintptr_t* mNodesAddress = reinterpret_cast<const std::uintptr_t*>(mNodes.data());
    int size = mNodes.size();
    ar & boost::serialization::make_nvp("mNodes_size", size);
    ar & boost::serialization::make_nvp("mNodes", boost::serialization::make_array(mNodesAddress, size));
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ContinuumElementIGA" << std::endl;
#endif
}

template void NuTo::ContinuumElementIGA<1>::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ContinuumElementIGA<2>::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ContinuumElementIGA<3>::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ContinuumElementIGA<1>::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ContinuumElementIGA<2>::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ContinuumElementIGA<3>::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ContinuumElementIGA<1>::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template void NuTo::ContinuumElementIGA<2>::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template void NuTo::ContinuumElementIGA<3>::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<int TDim>
template<class Archive>
void NuTo::ContinuumElementIGA<TDim>::load(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start deserialize ContinuumElementIGA " << std::endl;
#endif
    ar & boost::serialization::make_nvp("ContinuumElementIGA_ElementBase",boost::serialization::base_object<ElementBase >(*this));
    ar & boost::serialization::make_nvp("mSection", const_cast<SectionBase*&>(mSection));

    int size = 0;
    ar & boost::serialization::make_nvp("mNodes_size", size);
    std::uintptr_t* mNodesAddress = new std::uintptr_t[size];
    ar & boost::serialization::make_nvp("mNodes", boost::serialization::make_array(mNodesAddress, size));
    mNodes.assign(reinterpret_cast<NodeBase**>(&mNodesAddress[0]), reinterpret_cast<NodeBase**>(&mNodesAddress[size]));
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish deserialize ContinuumElementIGA" << std::endl;
#endif
}

template<int TDim>
void NuTo::ContinuumElementIGA<TDim>::SetNodePtrAfterSerialization(const std::map<std::uintptr_t, std::uintptr_t>& mNodeMapCast)
{
    for(std::vector<NodeBase*>::iterator it = mNodes.begin(); it != mNodes.end(); it++)
    {
        std::map<std::uintptr_t, std::uintptr_t>::const_iterator itCast = mNodeMapCast.find(reinterpret_cast<std::uintptr_t>(*it));
        if (itCast!=mNodeMapCast.end())
        {
            *it = reinterpret_cast<NodeBase*>(itCast->second);
        }
        else
            throw MechanicsException(__PRETTY_FUNCTION__, "The NodeBase-Pointer could not be updated.");
    }
}

BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::ContinuumElementIGA<1>)), "ContinuumElementIGA_1")
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::ContinuumElementIGA<2>)), "ContinuumElementIGA_2")
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::ContinuumElementIGA<3>)), "ContinuumElementIGA_3")
#endif // ENABLE_SERIALIZATION
