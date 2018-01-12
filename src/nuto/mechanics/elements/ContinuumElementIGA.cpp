#include "nuto/base/ErrorEnum.h"

#include "nuto/math/FullMatrix.h"

#include "nuto/mechanics/elements/ContinuumElementIGA.h"
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
NuTo::Element::eElementType NuTo::ContinuumElementIGA<TDim>::GetEnumType() const
{
    return Element::eElementType::CONTINUUMELEMENTIGA;
}

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
    Eigen::MatrixXd nodeCoordinates = this->ExtractNodeValues(0, Node::eDof::COORDINATES);

    Eigen::VectorXd volume(this->GetNumIntegrationPoints());
    for (int theIP = 0; theIP < this->GetNumIntegrationPoints(); theIP++)
    {
        Eigen::MatrixXd derivativeShapeFunctionsNatural = this->mInterpolationType->Get(Node::eDof::COORDINATES).CalculateDerivativeShapeFunctionsNatural(theIP, mKnotIDs);
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

    Eigen::MatrixXd derivativeShapeFunctions = this->mInterpolationType->Get(Node::eDof::COORDINATES).CalculateDerivativeShapeFunctionsNatural(theIP, mKnotIDs);
    Eigen::MatrixXd nodeCoordinates = this->ExtractNodeValues(0, Node::eDof::COORDINATES);
    double detJacobian = this->CalculateJacobian(derivativeShapeFunctions, nodeCoordinates).determinant();
    if (detJacobian < 0)
    {
        this->ReorderNodes();
        // recalculate node coordinates after reordering
        nodeCoordinates = this->ExtractNodeValues(0, Node::eDof::COORDINATES);
    }

    double size = 0;
    for (int iIP = 0; iIP < numIntegrationPoints; ++iIP)
    {
        Eigen::MatrixXd derivativeShapeFunctions = this->mInterpolationType->Get(Node::eDof::COORDINATES).CalculateDerivativeShapeFunctionsNatural(iIP, mKnotIDs);
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
    const Eigen::MatrixXd& derivativeShapeFunctionsGeometryNatural = this->mInterpolationType->Get(Node::eDof::COORDINATES).CalculateDerivativeShapeFunctionsNatural(rTheIP, mKnotIDs);

    Eigen::MatrixXd jacobian = this->CalculateJacobian(derivativeShapeFunctionsGeometryNatural, rData.mNodalValues[Node::eDof::COORDINATES]);

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
        if (dof == Node::eDof::COORDINATES)
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

template<int TDim>
Eigen::VectorXd NuTo::ContinuumElementIGA<TDim>::InterpolateDofGlobalCurrentConfiguration(int rTimeDerivative, const Eigen::VectorXd& rNaturalCoordinates, Node::eDof rDofTypeInit, Node::eDof rDofTypeCurrent) const
{
    const InterpolationBase& interpolationTypeInit = this->mInterpolationType->Get(rDofTypeInit);
    Eigen::VectorXd nodalInit    = this->ExtractNodeValues(rTimeDerivative, rDofTypeInit);
    Eigen::VectorXd nodalCurrent = this->ExtractNodeValues(rTimeDerivative, rDofTypeCurrent);

    Eigen::MatrixXd matrixN = interpolationTypeInit.CalculateMatrixN(rNaturalCoordinates);

    return matrixN * (nodalInit + nodalCurrent);
}

template<int TDim>
Eigen::VectorXd NuTo::ContinuumElementIGA<TDim>::InterpolateDofGlobalSurfaceDerivative(int rTimeDerivative, const Eigen::VectorXd& rParameter, int rDerivative, int rDirection) const
{
    Eigen::VectorXd nodalInitial       = this->ExtractNodeValues(rTimeDerivative, Node::eDof::COORDINATES);
    Eigen::VectorXd nodalDisplacements = this->ExtractNodeValues(rTimeDerivative, Node::eDof::DISPLACEMENTS);

    Eigen::MatrixXd matrixNDerivativeCoordinates = this->mInterpolationType->Get(Node::eDof::COORDINATES).CalculateMatrixNDerivative(rParameter, mKnotIDs, rDerivative, rDirection);
    Eigen::MatrixXd matrixNDerivativeDisplacements = this->mInterpolationType->Get(Node::eDof::DISPLACEMENTS).CalculateMatrixNDerivative(rParameter, mKnotIDs, rDerivative, rDirection);

    return matrixNDerivativeCoordinates * nodalInitial + matrixNDerivativeDisplacements * nodalDisplacements;
}

template<int TDim>
Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic> NuTo::ContinuumElementIGA<TDim>::InterpolateDofGlobalSurfaceDerivativeTotal(int rTimeDerivative, const Eigen::VectorXd& rParameter, int rDerivative, int rSurface) const
{
    Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic> derivative(1,1);

    int direction = -1;
    switch(rSurface)
    {
    case 0: direction = 0;
        break;
    case 1: direction = 1;
        break;
    case 2: direction = 0;
        break;
    case 3: direction = 1;
        break;
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "IGA2D has exactly four surfaces, 0 to 3. You tried to access " + std::to_string(rSurface) + ".");
    }

    derivative(0,0) = InterpolateDofGlobalSurfaceDerivative(rTimeDerivative, rParameter, rDerivative, direction);

    return derivative;
}

namespace NuTo // template specialization in *.cpp somehow requires the definition to be in the namespace...
{

template<>
double ContinuumElementIGA<1>::CalculateDetJxWeightIPxSection(double rDetJacobian, int rTheIP) const
{
    Eigen::MatrixXd matrixN = mInterpolationType->Get(Node::eDof::COORDINATES).CalculateMatrixN(rTheIP, mKnotIDs);
    Eigen::VectorXd globalIPCoordinate = matrixN * this->ExtractNodeValues(0, Node::eDof::COORDINATES);

    return rDetJacobian * mElementData->GetIntegrationType()->GetIntegrationPointWeight(rTheIP) * mSection->GetArea() * mSection->AsSectionTruss()->GetAreaFactor(globalIPCoordinate(0, 0));
}

template<>
double ContinuumElementIGA<2>::CalculateDetJxWeightIPxSection(double rDetJacobian, int rTheIP) const
{
    return rDetJacobian * mElementData->GetIntegrationType()->GetIntegrationPointWeight(rTheIP) * mSection->GetThickness();
}

template<>
Eigen::VectorXd ContinuumElementIGA<1>::InterpolateDofGlobalSurfaceNormal(const Eigen::VectorXd& rParameter) const
{
    throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Normal of a 1D element!");
}

template<>
Eigen::VectorXd ContinuumElementIGA<2>::InterpolateDofGlobalSurfaceNormal(const Eigen::VectorXd& rParameter) const
{
    Eigen::VectorXd tangent = InterpolateDofGlobalSurfaceDerivative(0, rParameter, 1, 0);

    if(tangent.rows() != 2)
        throw MechanicsException(__PRETTY_FUNCTION__, "The normal only available for 2D domains.");

    Eigen::Vector2d normal(-tangent(1), tangent(0));
    normal.normalize();

    return normal;
}

template<>
Eigen::VectorXd ContinuumElementIGA<1>::CalculateJacobianSurface(const Eigen::VectorXd &rParameter, const Eigen::VectorXd &rNodalCoordinates, int rSurfaceId) const
{
    throw MechanicsException(__PRETTY_FUNCTION__, "There is no surface in 1D case.");
}

template<>
Eigen::VectorXd ContinuumElementIGA<2>::CalculateJacobianSurface(const Eigen::VectorXd &rParameter, const Eigen::VectorXd &rNodalCoordinates, int rSurfaceId) const
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

template<>
const ContinuumElementIGA<1>& ContinuumElementIGA<1>::AsContinuumElementIGA1D() const
{
    return *this;
}

template<>
const ContinuumElementIGA<2>& ContinuumElementIGA<2>::AsContinuumElementIGA2D() const
{
    return *this;
}

template<>
ContinuumElementIGA<1>& ContinuumElementIGA<1>::AsContinuumElementIGA1D()
{
    return *this;
}

template<>
ContinuumElementIGA<2>& ContinuumElementIGA<2>::AsContinuumElementIGA2D()
{
    return *this;
}

}  // namespace NuTo

template class NuTo::ContinuumElementIGA<1>;
template class NuTo::ContinuumElementIGA<2>;


#ifdef ENABLE_SERIALIZATION
template void NuTo::ContinuumElementIGA<1>::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ContinuumElementIGA<2>::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ContinuumElementIGA<1>::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ContinuumElementIGA<2>::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ContinuumElementIGA<1>::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ContinuumElementIGA<2>::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template<int TDim>
template<class Archive>
void NuTo::ContinuumElementIGA<TDim>::save(Archive & ar, const unsigned int version)const
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ContinuumElementIGA " << std::endl;
#endif
    ar & boost::serialization::make_nvp("ContinuumElement",boost::serialization::base_object<ContinuumElement<TDim> >(*this));
    ar & boost::serialization::make_nvp("mKnots_size", mKnots);
    ar & boost::serialization::make_nvp("mKnotIDs_size", mKnotIDs);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ContinuumElementIGA" << std::endl;
#endif
}

template void NuTo::ContinuumElementIGA<1>::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ContinuumElementIGA<2>::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ContinuumElementIGA<1>::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ContinuumElementIGA<2>::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ContinuumElementIGA<1>::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template void NuTo::ContinuumElementIGA<2>::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<int TDim>
template<class Archive>
void NuTo::ContinuumElementIGA<TDim>::load(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start deserialize ContinuumElementIGA " << std::endl;
#endif
    ar & boost::serialization::make_nvp("ContinuumElementIGA_ElementBase",boost::serialization::base_object<ContinuumElement<TDim> >(*this));
    ar & BOOST_SERIALIZATION_NVP(mKnots);
    ar & BOOST_SERIALIZATION_NVP(mKnotIDs);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish deserialize ContinuumElementIGA" << std::endl;
#endif
}

BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::ContinuumElementIGA<1>)), "ContinuumElementIGA_1")
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::ContinuumElementIGA<2>)), "ContinuumElementIGA_2")
#endif // ENABLE_SERIALIZATION