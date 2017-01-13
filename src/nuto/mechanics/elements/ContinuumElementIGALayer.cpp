#include "nuto/base/ErrorEnum.h"

#include "nuto/math/FullMatrix.h"

#include "nuto/mechanics/elements/ContinuumElementIGALayer.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/nodes/NodeEnum.h"

#include "nuto/mechanics/sections/SectionTruss.h"
#include "nuto/mechanics/sections/SectionPlane.h"
#include "nuto/mechanics/elements/ElementOutputBase.h"
#include "nuto/mechanics/elements/ElementOutputIpData.h"
#include "nuto/mechanics/elements/ElementDataBase.h"

#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"
#include "nuto/mechanics/interpolationtypes/InterpolationBase.h"
#include "nuto/mechanics/interpolationtypes/InterpolationType.h"

#include "nuto/mechanics/elements/ElementEnum.h"
#include "nuto/mechanics/elements/EvaluateDataContinuum.h"
#include "nuto/mechanics/elements/IpDataEnum.h"

#include "nuto/mechanics/dofSubMatrixStorage/BlockFullVector.h"
#include "nuto/mechanics/dofSubMatrixStorage/BlockFullMatrix.h"

#include "nuto/mechanics/constitutive/ConstitutiveBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOBase.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveScalar.h"
#include "nuto/mechanics/constitutive/inputoutput/EngineeringStrain.h"
#include "nuto/mechanics/constitutive/inputoutput/EngineeringStress.h"

#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/elements/ElementBase.h"

template<int TDim>
NuTo::ContinuumElementIGALayer<TDim>::ContinuumElementIGALayer(const NuTo::StructureBase*          rStructure,
                                                               const std::vector<NuTo::NodeBase*> &rNodes,
                                                               const Eigen::MatrixXd              &rKnots,
                                                               const Eigen::VectorXi              &rKnotIDs,
                                                               ElementData::eElementDataType       rElementDataType,
                                                               IpData::eIpDataType                 rIpDataType,
                                                               InterpolationType                  *rInterpolationType)

    :  ContinuumElementIGA<TDim>(rStructure, rNodes, rKnots, rKnotIDs, rElementDataType, rIpDataType, rInterpolationType)
{}


////! @brief calculates output data for the element
////! @param rInput ... constitutive input map for the constitutive law
////! @param rOutput ...  coefficient matrix 0 1 or 2  (mass, damping and stiffness) and internal force (which includes inertia terms)
template<int TDim>
NuTo::eError NuTo::ContinuumElementIGALayer<TDim>::Evaluate(const ConstitutiveInputMap& rInput, std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>>& rOutput)
{
    // So far the Layer element is used for contact element containing a list of slave and master elements
    for (auto dof : this->mStructure->GetDofStatus().GetActiveDofTypes())
    {
        if (not (this->mInterpolationType->IsDof(dof)) && dof != Node::eDof::DISPLACEMENTS)
        {
//            rGlobalRowDofs[dof].Resize(0);
            continue;
        }

        const InterpolationBase& interpolationType = this->GetInterpolationType()->Get(dof);
        const int numNodes = interpolationType.GetNumNodes();

        FullVector<int, Eigen::Dynamic> dofWiseGlobalRowDofs;

        int numSlaveDofs = interpolationType.GetNumDofs();
        dofWiseGlobalRowDofs.setZero(numSlaveDofs);

        unsigned int numDofsPerType = this->GetNode(interpolationType.GetNodeIndex(0))->GetNum(dof);

        for (int iNodeDof = 0; iNodeDof < numNodes; ++iNodeDof)
        {
            const NodeBase* nodePtr = this->GetNode(interpolationType.GetNodeIndex(iNodeDof));

            for (unsigned iDof = 0; iDof < numDofsPerType; ++iDof)
            {
                dofWiseGlobalRowDofs[numDofsPerType * iNodeDof + iDof] = nodePtr->GetDof(dof, iDof);
            }
        }
    }

    return eError::SUCCESSFUL;
}

template<int TDim>
NuTo::Element::eElementType NuTo::ContinuumElementIGALayer<TDim>::GetEnumType() const
{
    return Element::eElementType::CONTINUUMELEMENTIGA;
}

template<int TDim>
Eigen::MatrixXd NuTo::ContinuumElementIGALayer<TDim>::CalculateJacobian(const Eigen::MatrixXd& rDerivativeShapeFunctions, const Eigen::VectorXd& rNodeCoordinates) const
{
    int numCoordinateNodes = this->GetNumNodes(Node::eDof::COORDINATES);
    assert(rDerivativeShapeFunctions.rows() == numCoordinateNodes);
    assert(rDerivativeShapeFunctions.cols() == TDim);

    // boundary layer or full iga element
    assert(rNodeCoordinates.rows() == (TDim+1)*this->GetNumNodes(Node::eDof::COORDINATES) );

    Eigen::MatrixXd nodeBlockCoordinates(TDim+1, numCoordinateNodes);
    // convert the coordinates to a block structure
    // x0  x1  x1  x2 ...
    // y0  y1  y2  y3 ...
    // z0  z1  z2  z3 ...
    for (int i = 0; i < numCoordinateNodes; ++i)
        nodeBlockCoordinates.col(i) = rNodeCoordinates.block<TDim+1, 1>((TDim+1) * i, 0);

    return nodeBlockCoordinates.lazyProduct(rDerivativeShapeFunctions);
}


template<int TDim>
Eigen::VectorXd NuTo::ContinuumElementIGALayer<TDim>::CalculateJacobianSurface(const Eigen::VectorXd &rParameter, const Eigen::VectorXd &rNodalCoordinates, int rSurfaceId) const
{
    const InterpolationBase&  interpolationTypeCoords =  this->mInterpolationType->Get(Node::eDof::COORDINATES);
    Eigen::MatrixXd derivativeShapeFunctionsNaturalSlave =  interpolationTypeCoords.CalculateDerivativeShapeFunctionsNatural(rParameter);
    const Eigen::MatrixXd jacobianStd = CalculateJacobian(derivativeShapeFunctionsNaturalSlave, rNodalCoordinates);// = [dX / dXi]
    // in case of non IGA - just an identity matrix
    const Eigen::MatrixXd jacobianIGA = CalculateJacobianParametricSpaceIGA();// = [dXi / d\tilde{Xi}]
    Eigen::VectorXd ipCoordsSurface(1);
    Eigen::MatrixXd derivativeNaturalSurfaceCoordinates =  interpolationTypeCoords.CalculateDerivativeNaturalSurfaceCoordinates(ipCoordsSurface, rSurfaceId); // = [dXi / dAlpha]

    return jacobianStd * jacobianIGA * derivativeNaturalSurfaceCoordinates; // = || [dX / dXi] * [dXi / dAlpha] ||

}

template<int TDim>
Eigen::Matrix<double, TDim, TDim> NuTo::ContinuumElementIGALayer<TDim>::CalculateJacobianParametricSpaceIGA() const
{
    Eigen::Matrix<double, TDim, TDim> jac;
    jac.setZero(TDim, TDim);
    for(int i = 0; i < TDim; i++) jac(i,i) = 0.5*(this->mKnots(i,1) - this->mKnots(i,0));

    return jac;
}

template<int TDim>
NuTo::NodeBase* NuTo::ContinuumElementIGALayer<TDim>::GetNode(int rLocalNodeNumber)
{
    assert(rLocalNodeNumber >= 0);
    assert(rLocalNodeNumber < (int )((this->mNodes).size()));
    return (this->mNodes)[rLocalNodeNumber];
}

template<int TDim>
const NuTo::NodeBase* NuTo::ContinuumElementIGALayer<TDim>::GetNode(int rLocalNodeNumber) const
{
    assert(rLocalNodeNumber >= 0);
    assert(rLocalNodeNumber < (int )((this->mNodes).size()));
    return (this->mNodes)[rLocalNodeNumber];
}

template<int TDim>
void NuTo::ContinuumElementIGALayer<TDim>::SetNode(int rLocalNodeNumber, NodeBase* rNode)
{
    assert(rLocalNodeNumber >= 0);
    assert(rLocalNodeNumber < (int )((this->mNodes).size()));
    assert(rNode != nullptr);
    (this->mNodes)[rLocalNodeNumber] = rNode;
}

template<int TDim>
Eigen::VectorXd NuTo::ContinuumElementIGALayer<TDim>::InterpolateDofGlobal(int rTimeDerivative, const Eigen::VectorXd& rNaturalCoordinates, Node::eDof rDofType) const
{
    const InterpolationBase& interpolationType = this->mInterpolationType->Get(rDofType);
    Eigen::MatrixXd nodalValues = this->ExtractNodeValues(rTimeDerivative, rDofType);
    Eigen::MatrixXd matrixN = interpolationType.CalculateMatrixN(rNaturalCoordinates, this->mKnotIDs);

//    std::cout << matrixN * nodalValues << std::endl;

    return matrixN * nodalValues;
}

template<int TDim>
Eigen::VectorXd NuTo::ContinuumElementIGALayer<TDim>::InterpolateDofGlobalCurrentConfiguration(int rTimeDerivative, const Eigen::VectorXd& rNaturalCoordinates, Node::eDof rDofTypeInit, Node::eDof rDofTypeCurrent) const
{
    const InterpolationBase& interpolationTypeInit = this->mInterpolationType->Get(rDofTypeInit);
    Eigen::VectorXd nodalInit    = this->ExtractNodeValues(rTimeDerivative, rDofTypeInit);
    Eigen::VectorXd nodalCurrent = this->ExtractNodeValues(rTimeDerivative, rDofTypeCurrent);

    Eigen::MatrixXd matrixN = interpolationTypeInit.CalculateMatrixN(rNaturalCoordinates, this->mKnotIDs);

    return matrixN * (nodalInit + nodalCurrent);
}

template<int TDim>
Eigen::VectorXd NuTo::ContinuumElementIGALayer<TDim>::InterpolateDofGlobalSurfaceDerivative(int rTimeDerivative, const Eigen::VectorXd& rParameter, int rDerivative, int rDirection) const
{
    Eigen::VectorXd nodalInitial       = this->ExtractNodeValues(rTimeDerivative, Node::eDof::COORDINATES);
    Eigen::VectorXd nodalDisplacements = this->ExtractNodeValues(rTimeDerivative, Node::eDof::DISPLACEMENTS);

    Eigen::MatrixXd matrixNDerivativeCoords = this->mInterpolationType->Get(Node::eDof::COORDINATES).CalculateMatrixNDerivative(rParameter, this->mKnotIDs, rDerivative, rDirection);
    Eigen::MatrixXd matrixNDerivativeDisp = this->mInterpolationType->Get(Node::eDof::DISPLACEMENTS).CalculateMatrixNDerivative(rParameter, this->mKnotIDs, rDerivative, rDirection);

    return matrixNDerivativeCoords * nodalInitial + matrixNDerivativeDisp * nodalDisplacements;
}

namespace NuTo // template specialization in *.cpp somehow requires the definition to be in the namespace...
{

template<>
Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic> NuTo::ContinuumElementIGALayer<1>::InterpolateDofGlobalSurfaceDerivativeTotal(int rTimeDerivative, const Eigen::VectorXd& rParameter, int rDerivative) const
{
    Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic> derivative(1,1);

    derivative(0,0) = InterpolateDofGlobalSurfaceDerivative(rTimeDerivative, rParameter, rDerivative, 0);

    return derivative;
}

template<>
Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic> NuTo::ContinuumElementIGALayer<2>::InterpolateDofGlobalSurfaceDerivativeTotal(int rTimeDerivative, const Eigen::VectorXd& rParameter, int rDerivative) const
{
    // think of it somehow as of a tensor ....
    Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic> derivative;
    switch (rDerivative)
    {
    case 0:
    {
        derivative.resize(1,1);
        derivative(0,0) = InterpolateDofGlobalSurfaceDerivative(rTimeDerivative, rParameter, rDerivative, 0);
        break;
    }
    case 1:
    {
        derivative.resize(1,2);
        derivative(0,0) = InterpolateDofGlobalSurfaceDerivative(rTimeDerivative, rParameter, rDerivative, 0);
        derivative(0,1) = InterpolateDofGlobalSurfaceDerivative(rTimeDerivative, rParameter, rDerivative, 1);
        break;
    }
    case 2:
    {
        derivative.resize(2,2);
        derivative(0,0) = InterpolateDofGlobalSurfaceDerivative(rTimeDerivative, rParameter, rDerivative, 0);
        derivative(0,1) = InterpolateDofGlobalSurfaceDerivative(rTimeDerivative, rParameter, rDerivative, 2);
        derivative(1,0) = InterpolateDofGlobalSurfaceDerivative(rTimeDerivative, rParameter, rDerivative, 1);
        derivative(1,1) = InterpolateDofGlobalSurfaceDerivative(rTimeDerivative, rParameter, rDerivative, 2);
        break;
    }
    default:
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Maximum derivative is of order 2!");
        break;
    }

    return derivative;
}

template<>
const ContinuumElementIGALayer<1>& ContinuumElementIGALayer<1>::AsContinuumElementIGALayer1D() const
{
    return *this;
}

template<>
const ContinuumElementIGALayer<2>& ContinuumElementIGALayer<2>::AsContinuumElementIGALayer2D() const
{
    return *this;
}

template<>
ContinuumElementIGALayer<1>& ContinuumElementIGALayer<1>::AsContinuumElementIGALayer1D()
{
    return *this;
}

template<>
ContinuumElementIGALayer<2>& ContinuumElementIGALayer<2>::AsContinuumElementIGALayer2D()
{
    return *this;
}

}  // namespace NuTo

template class NuTo::ContinuumElementIGALayer<1>;
template class NuTo::ContinuumElementIGALayer<2>;


#ifdef ENABLE_SERIALIZATION
template void NuTo::ContinuumElementIGALayer<1>::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ContinuumElementIGALayer<2>::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ContinuumElementIGALayer<1>::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ContinuumElementIGALayer<2>::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ContinuumElementIGALayer<1>::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ContinuumElementIGALayer<2>::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template<int TDim>
template<class Archive>
void NuTo::ContinuumElementIGALayer<TDim>::save(Archive & ar, const unsigned int version)const
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ContinuumElementIGALayer " << std::endl;
#endif
    ar & boost::serialization::make_nvp("ContinuumElementIGALayer_Base",boost::serialization::base_object<ContinuumElementIGA<TDim> >(*this));
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ContinuumElementIGALayer" << std::endl;
#endif
}

template void NuTo::ContinuumElementIGALayer<1>::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ContinuumElementIGALayer<2>::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ContinuumElementIGALayer<1>::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ContinuumElementIGALayer<2>::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ContinuumElementIGALayer<1>::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template void NuTo::ContinuumElementIGALayer<2>::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<int TDim>
template<class Archive>
void NuTo::ContinuumElementIGALayer<TDim>::load(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start deserialize ContinuumElementIGALayer " << std::endl;
#endif
    ar & boost::serialization::make_nvp("ContinuumElementIGALayer_Base",boost::serialization::base_object<ContinuumElementIGA<TDim> >(*this));
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish deserialize ContinuumElementIGALayer" << std::endl;
#endif
}

BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::ContinuumElementIGALayer<1>)), "ContinuumElementIGALayer_1")
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::ContinuumElementIGALayer<2>)), "ContinuumElementIGALayer_2")
#endif // ENABLE_SERIALIZATION
