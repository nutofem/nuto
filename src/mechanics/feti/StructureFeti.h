
#pragma once

#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/structures/StructureOutputBlockMatrix.h"
#include "mechanics/nodes/NodeEnum.h"

#include <eigen3/Eigen/Sparse>
#include <cstring>
#include <fstream>
#include <set>
#include <mechanics/MechanicsEnums.h>
#include "mechanics/mesh/MeshGenerator.h"


namespace NuTo
{


/**     \author Philip Huschke
 *      \date September 2016
 *      \brief Class for FETI methods
 *
 *
 *
 *      The standard interface problem for FETI methods
 *
 *      \f[
 *      \begin{bmatrix}
 *      \sum_{s=1}^{N_s} \boldsymbol{B}_s \boldsymbol{K}_s^+ \boldsymbol{B}_s^T & \boldsymbol{B}_{1} \boldsymbol{R}_{1}
 * & \dots & \boldsymbol{B}_{{N_s}} \boldsymbol{R}_{{N_s}}
 *      \\
 *      \boldsymbol{R}_{1}^T \boldsymbol{B}_{1}^T  & \boldsymbol{0} & \dots & \boldsymbol{0}
 *      \\
 *      \vdots & \vdots & \ddots & \vdots
 *      \\
 *      \boldsymbol{R}_{{N_s}}^T \boldsymbol{B}_{{N_s}}^T & \boldsymbol{0} & \dots & \boldsymbol{0}
 *      \end{bmatrix}
 *      \begin{bmatrix}
 *      \boldsymbol{\lambda}
 *      \\
 *      \boldsymbol{\alpha}_{1}
 *      \\
 *      \vdots
 *      \\
 *      \boldsymbol{\alpha}_{{N_s}}
 *      \end{bmatrix}
 *      =
 *      \begin{bmatrix}
 *      \sum_{s=1}^{N_s} \boldsymbol{B}_s \boldsymbol{K}_s^+ \boldsymbol{f}_s
 *      \\
 *      \boldsymbol{R}_{{1}}^T \boldsymbol{f}_{{1}}
 *      \\
 *      \vdots
 *      \\
 *      \boldsymbol{R}_{{N_s}}^T \boldsymbol{f}_{{N_s}}
 *      \end{bmatrix}
 *      \label{eq_interface_problem}
 *      \f]
 *
 *
 *      \f$ \boldsymbol{K}_s \f$ stiffness matrix of subdomain s (possibly singular)
 *
 *
 *      \f$ \boldsymbol{B}_s \f$ StructureFETI::mConnectivityMatrix of subdomain s
 *
 *
 *      \f$ \boldsymbol{R}_s \f$ StructureFETI::mRigidBodyModes of subdomain s
 *
 *
 *      \f$ \boldsymbol{\lambda} \f$ Lagrange multiplier
 *
 *
 *      \f$ \boldsymbol{\alpha}_s \f$ Linear combination of rigid body modes of subdomain s
 *
 *
 *      \f$ N_s \f$ Number of subdomains StructureFETI::mNumProcesses
 *
 *
 *
 *
 */
class StructureFeti : public Structure
{
private:
    ///
    /// \brief The Node struct
    ///
    struct Node
    {
        ///
        /// \brief mCoordinates
        ///
        Eigen::Vector3d mCoordinates;
        int mId;
    };

    struct Element
    {
        int mId;
        std::vector<int> mNodeIds;
    };

    struct Interface
    {
        std::map<int, int> mNodeIdsMap;
        int mValue;
    };

    using NodeList = std::vector<Node>;
    using ElementList = std::vector<Element>;
    using InterfaceList = std::vector<Interface>;
    using Matrix = Eigen::MatrixXd;
    using SparseMatrix = Eigen::SparseMatrix<double>;
    using SparseMatrixMap = std::map<NuTo::Node::eDof, SparseMatrix>;


public:
    //! @brief Constructor
    //! @param rDimension   Structural dimension (1,2 or 3)
    StructureFeti(int rDimension);

    ///
    /// \brief mRank
    ///
    /// MPI_Comm_rank(..)
    const int mRank;

    ///
    /// \brief mNumProcesses
    ///
    /// MPI_Comm_size(...)
    const int mNumProcesses;

    ///
    /// \brief FindKeywordInFile
    /// \param file
    /// \param keyword
    ///
    void FindKeywordInFile(std::ifstream& file, std::string keyword);

    ///
    /// \brief GetRigidBodyModes
    /// \return mRigidBodyModes
    ///
    Matrix& GetRigidBodyModes()
    {
        return mRigidBodyModes;
    }

    const int GetNumRigidBodyModes()
    {
        return mNumRigidBodyModes;
    }
    const int GetNumRigidBodyModesTotal()
    {
        return mNumRigidBodyModesTotal;
    }

    ///
    /// \brief GetInterfaceRigidBodyModes
    /// \return
    ///
    Matrix& GetInterfaceRigidBodyModes()
    {
        return mInterfaceRigidBodyModes;
    }


    Matrix& GetG()
    {
        return mG;
    }

    ///
    /// \brief GetConnectivityMatrix
    /// \return mConnectivityMatrix
    ///
    SparseMatrix& GetConnectivityMatrix()
    {
        return mConnectivityMatrix;
    }


    Matrix& GetProjectionMatrix()
    {
        return mProjectionMatrix;
    }

    ///
    /// \brief AssembleConnectivityMatrix
    ///
    /// Assembles StructureFETI::mConnectivityMatrix
    ///
    /// First it loops over StructureFETI::mInterfaces and inserts
    ///
    /// +1 if StructureFETI::mRank is lower than its neighbour
    /// -1 if StructureFETI::mRank is greater than its neighbour
    ///
    /// this guarantees the continuity of the field variable across interfaces.
    ///
    /// Next it loops over StructureFETI::mBoundaryDofIds and inserts +1 to constraint the corresponding DOF
    ///
    void AssembleConnectivityMatrix();

    ///
    /// \brief AssembleBoundaryDofIds
    ///
    void AssembleBoundaryDofIds();


    ///
    /// \brief ApplyVirtualConstraints
    /// \param nodeIdsBoundaries
    /// \param nodeIdsLoads
    void ApplyVirtualConstraints(const std::vector<int>& nodeIdsBoundaries, const std::vector<int>& nodeIdsLoads);


    ///
    /// \brief ImportMeshJson
    ///
    /// Json mesh file generated with <a href="https://github.com/dbeurle/GmshReader">GmshReader</a>
    ///
    /// \param rFileName
    /// \param interpolationTypeId
    ///
    void ImportMeshJson(std::string rFileName, const int interpolationTypeId);


    ///
    /// \brief IsFloating
    /// \return
    ///
    bool IsFloating()
    {
        return mIsFloating;
    }

    ///
    /// \brief SetIsFloating
    /// \param isFloating
    ///
    void SetIsFloating(const bool isFloating)
    {
        mIsFloating = isFloating;
    }

    /**
    *   \brief CalculateRigidBodyModes
    *
    *   \f[
    *   R =
    *   \begin{bmatrix}
    *   -\boldsymbol{K}^{-1}_{11} \boldsymbol{K}_{12}
    *   \\
    *   \boldsymbol{I}
    *   \end{bmatrix}
    *   \f]
    *
    */
    void CalculateRigidBodyModes();

    /**
    *   \brief CalculateRigidBodyModesTotalFETI
    *
    *   In 2D:
    *
    *   \f[
    *   R =
    *   \begin{bmatrix}
    *   1 & 0 & -y_0
    *   \\
    *   0 & 1 &  x_0
    *   \\
    *   \vdots & \vdots &  \vdots
    *   \\
    *   1 & 0 & -y_N
    *   \\
    *   0 & 1 &  x_N
    *   \end{bmatrix}
    *   \f]
    *
    */
    void CalculateRigidBodyModesTotalFETI();

    ///
    /// \brief CheckRigidBodyModes
    /// \param tolerance
    ///
    /// Checks if the calculated rigid body modes are in the null space of \f$ \boldsymbol{K}_s \f$
    ///
    void CheckRigidBodyModes(const StructureOutputBlockMatrix hessian0, const double tolerance = 1.e-9) const;

    ///
    /// \brief CheckProjectionMatrix
    /// \param tolerance
    ///
    /// Checks if the calculated projection matrix satisfies \f$ \boldsymbol{P} -\boldsymbol{P} \boldsymbol{P} =
    /// \boldsymbol{0} \f$
    ///
    void CheckProjectionMatrix(const double tolerance = 1.e-9) const;

    ///
    /// \brief CheckProjectionMatrix
    /// \param tolerance
    ///
    /// Checks if the calculated projection matrix satisfies \f$ \boldsymbol{P}^T -\boldsymbol{G} = \boldsymbol{0} \f$
    ///
    void CheckProjectionOfCoarseGrid(const double tolerance = 1.e-9) const;

    ///
    /// \brief CheckStiffnessPartitioning
    /// \param hessian0
    /// \param tolerance
    ///
    void CheckStiffnessPartitioning(const StructureOutputBlockMatrix hessian0, const double tolerance = 1.e-9) const;
    void CreateDummy1D();
    void CalculateInterfaceRigidBodyModes();


    ///
    /// \brief Constraints all degrees of freedom of nodes in node group
    /// \param nodeGroupId
    ///
    void ApplyConstraintsTotalFeti(const int nodeGroupId);

    ///
    /// \brief Constraints all degrees of freedom
    /// \param dofIds to be constrained
    ///
    void ApplyConstraintsTotalFeti(const std::vector<int>& dofIds);

    ///
    /// \param dofIdAndPrescribedDisplacementMap
    void ApplyPrescribedDisplacements(const std::map<int, double> dofIdAndPrescribedDisplacementMap);

    ///
    void CalculateProjectionMatrix();

    ///
    void CalculateG();

    ///
    void CalculateAndAppyVirtualConstraints();

    ///
    /// \param interfaceCoordinates
    /// \param globalNodeId
    /// \param interfaceId
    /// \param connectivityValue
    void AddNodeIdsToInterface(const std::vector<double>& interfaceCoordinates, int& globalNodeId, int& interfaceId,
                               const double connectivityValue)
    {
        const int groupNodesInterface = GroupCreate(eGroupId::Nodes);
        GroupAddNodeCoordinateRange(groupNodesInterface, eDirection::X, interfaceCoordinates[0] - 1.e-6, interfaceCoordinates[0] + 1.e-6);


        std::vector<int> nodeIdsInterface = GroupGetMemberIds(groupNodesInterface);

//        for (const auto nodeId : nodeIdsInterface)
//            GetLogger() << nodeId << "\n";

        mInterfaces[interfaceId].mValue = connectivityValue;


        for (size_t i = 0; i < nodeIdsInterface.size(); ++i)
        {
            mInterfaces[interfaceId].mNodeIdsMap.emplace(globalNodeId, nodeIdsInterface[i]);
            ++globalNodeId;
        }
    }

    ///
    /// \param meshDimensions
    /// \param numElements
    /// \return
    std::pair<int,int> CreateRectangularMesh2D(const std::vector<double>& meshDimensions, const std::vector<int>& numElements)
    {

        assert(GetDimension() == static_cast<int>(meshDimensions.size()) and "Dimensions mismatch");
        assert(GetDimension() == static_cast<int>(numElements.size()) and "Dimensions mismatch");

        mNumInterfaceNodesTotal = (mNumProcesses - 1) * (numElements[1] + 1);

        double deltaX = meshDimensions[0] / mNumProcesses;

        std::vector<double> endPoints;
        endPoints.push_back((mRank + 1) * deltaX);
        endPoints.push_back(meshDimensions[1]);

        std::vector<double> startPoints;
        startPoints.push_back(mRank * deltaX);
        startPoints.push_back(0);

        // create local meshes
        auto importContainer = NuTo::MeshGenerator::Grid(*this, startPoints, endPoints, numElements,
                                                         NuTo::Interpolation::eShapeType::QUAD2D);




        int numInterfaces;
        int globalNodeId = mRank == 0 ? 0 : (mRank - 1) * (numElements[1] + 1);

        if (mRank == 0 or mRank == mNumProcesses - 1)
        {
            numInterfaces = 1;
        }
        else
        {
            numInterfaces = 2;
        }

        mInterfaces.resize(numInterfaces);
        int interfaceId = 0;
        if (mRank > 0)
        {
            AddNodeIdsToInterface(startPoints, globalNodeId, interfaceId, -1);

            ++interfaceId;
        }

        if (mRank < mNumProcesses - 1)
        {
            AddNodeIdsToInterface(endPoints, globalNodeId, interfaceId, 1);
        }

        return importContainer;
    }

    ///
    /// \return
    const Eigen::VectorXd& GetPrescribedDofVector()
    {
        return mPrescribedDofVector;
    }



    ///
    /// \return
    std::vector<int> CalculateLagrangeMultiplierIds()
    {

        std::vector<int> lagrangeMultiplierDofIds;

        // add lagrange multipliers from interfaces
        for (const auto& interface : mInterfaces)
            for(const auto& nodeIdPair : interface.mNodeIdsMap)
            {
                /// \todo Think about other DOFs. They must be added too!
                auto dofIds = NodeGetDofIds(nodeIdPair.second, NuTo::Node::eDof::DISPLACEMENTS);

                for (const auto& id : dofIds)
                    lagrangeMultiplierDofIds.push_back(id);
            }

        lagrangeMultiplierDofIds.insert(lagrangeMultiplierDofIds.end(), mBoundaryDofIds.begin(), mBoundaryDofIds.end());

        lagrangeMultiplierDofIds.insert(lagrangeMultiplierDofIds.end(), mPrescribedDisplacementDofIds.begin(),
                                        mPrescribedDisplacementDofIds.end());





        return lagrangeMultiplierDofIds;


    }


    /// \brief Assembles vector for multiplicity scaling
    void MultiplicityScaling()
    {

        std::vector<int> multiplicity(mNumInterfaceNodesTotal,0);

        for (const auto& interface : mInterfaces)
            for (const auto& nodePair : interface.mNodeIdsMap)
            {
                const int globalIndex = nodePair.first;

                multiplicity[globalIndex]++;

            }

        MPI_Allreduce(MPI_IN_PLACE, multiplicity.data(), multiplicity.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        Eigen::DiagonalMatrix<double, Eigen::Dynamic> mFetiScalingMatrix(mNumInterfaceNodesTotal);

        for (int i = 0; i < mNumInterfaceNodesTotal; ++i)
            mFetiScalingMatrix.diagonal()[i] = multiplicity[i];

        GetLogger() << "mFetiScalingMatrix (multiplicity scaling) \n" << mFetiScalingMatrix << "\n\n";
    }

protected:
    bool mIsFloating = true;

    NodeList mNodes;
    ElementList mElements;
    InterfaceList mInterfaces;
    std::set<int> mSubdomainBoundaryNodeIds;


    /**
     *  \brief mG Natural coarse space
     *
     *  \f[
     *  G =
     *  \begin{bmatrix}
     *  B^T_1 R_1
     *  &
     *  B^T_2 R_2
     *  &
     *  \dots
     *  &
     *  B^T_{N_s} R_{N_s}
     *  \end{bmatrix}
     *  \f]
     *
     *  \f[
     *  G \in \mathbb{R}^{\text{number of Lagrange multipliers} \times \text{total number of rigid body modes}  }
     *  \f]
     *
     */
    Matrix mG;


    /**
     *   \brief mProjectionMatrix
     *
     *   \f[
     *   P = I - G (G^T G)^{-1} G^T
     *   \f]
     *
     *   \f[
     *   P \in \mathbb{R}^{\text{number of Lagrange multipliers} \times \text{number of Lagrange multipliers}  }
     *   \f]
     *
     */
    Matrix mProjectionMatrix;


    int mNumRigidBodyModes = -1337;
    int mNumRigidBodyModesTotal = -1337;
public:
    int mNumInterfaceNodesTotal = -1337;
protected:
    ///
    /// \brief mBoundaryDofIds
    ///
    /// An array of ids of all the degrees of freedom that are constraint in the subdomain
    ///
    std::vector<int> mBoundaryDofIds;


    Eigen::VectorXd mPrescribedDofVector;

    ///
    /// \brief mPrescribedDisplacementDofIds
    ///
    std::vector<int> mPrescribedDisplacementDofIds;

    ///
    /// \brief mGlobalStartIndexBoundaryDofIds
    ///
    /// This variable is important for the assembly of the connectivity matrix. It is equal to the sum of the boundary
    /// DOFs of all subdomains with lower rank. E.g.
    ///
    /// Subdomain 0:  mBoundaryDofIds.size() = 3: =>  mGlobalStartIndexBoundaryDofIds = 0:
    /// Subdomain 1:  mBoundaryDofIds.size() = 6: =>  mGlobalStartIndexBoundaryDofIds = 3:
    /// Subdomain 2:  mBoundaryDofIds.size() = 2: =>  mGlobalStartIndexBoundaryDofIds = 9:
    ///
    int mGlobalStartIndexBoundaryDofIds = -1337;


    int mGlobalStartIndexPrescribedDisplacementDofIds = -1337;

    ///
    /// \brief mNumTotalBoundaryDofIds
    ///
    /// The total number of degrees of freedom that constraint in all subdomains
    ///
    /// MPI_Allreduce(  &numLocalBoundaryDofIds,
    ///                 &structure.mNumTotalBoundaryDofIds,
    ///                 1,
    ///                 MPI_INT,
    ///                 MPI_SUM,
    ///                 MPI_COMM_WORLD);
    ///
    int mNumTotalBoundaryDofIds = -1337;

    ///
    /// \brief mNumTotalPrescribedDisplacementDofIds
    ///
    /// The total number of degrees of freedom that have a prescribed displacement
    ///
    ///
public:
    int mNumTotalPrescribedDisplacementDofIds = 0;
protected:
    ///
    /// \brief mRigidBodyModes
    ///
    /// A matrix that contains the rigid body modes of the subdomain
    ///
    ///
    /// \f[
    /// R_s \in \mathbb{R}^{\text{number of DOFs} \times \text{number of rigid body modes}}
    /// \f]
    ///
    Matrix mRigidBodyModes;

    StructureOutputBlockMatrix mRigidBodyModesBlockMatrix;

    ///
    /// \brief mInterfaceRigidBodyModes
    ///
    /// \f[
    /// B_s R_s \in \mathbb{R}^{\text{number of Lagrange multipliers} \times \text{number of rigid body modes}}
    /// \f]
    ///
    Matrix mInterfaceRigidBodyModes;

    ///
    /// \brief mConnectivityMatrix
    ///
    /// Sparse matrix that restricts the DOFs of one subdomain to its interface boundary. Contains only \f$ 0, 1, -1 \f$
    ///
    /// \f[
    /// B_s \in \mathbb{R}^{\text{number of Lagrange multipliers} \times \text{number of degrees of freedom}  }
    /// \f]
    ///
    SparseMatrix mConnectivityMatrix;
};
} // namespace NuTo
