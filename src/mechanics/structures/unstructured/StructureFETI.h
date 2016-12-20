
#pragma once

#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/nodes/NodeEnum.h"


#include <eigen3/Eigen/Sparse>
#include <cstring>
#include <fstream>
#include <set>


namespace NuTo
{


//! @author Philip Huschke
//! @date September 2016
//! @brief Class for FETI methods
class StructureFETI: public Structure
{
private:
    struct Node
    {
        Eigen::Vector3d mCoordinates;
        int mId;
    };

    struct Element
    {
        int mId;
        std::vector<int> mNodeIds;
    };

    struct Boundary
    {
        std::map<int, int> mNodeIdsMap;
        std::vector<int> mValues;
    };

    struct Interface
    {
        std::map<int, int> mNodeIdsMap;
        std::vector<int> mNodeIds;
        int mValue;
    };

    using NodeList          = std::vector<Node>;
    using ElementList       = std::vector<Element>;
    using BoundaryList      = std::vector<Boundary>;
    using InterfaceList     = std::vector<Interface>;
    using Matrix            = Eigen::MatrixXd;
    using SparseMatrix      = Eigen::SparseMatrix<double>;
    using SparseMatrixMap   = std::map<NuTo::Node::eDof, SparseMatrix>;


public:
    //! @brief Constructor
    //! @param rDimension   Structural dimension (1,2 or 3)
    StructureFETI(int rDimension);

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
    void FindKeywordInFile(std::ifstream &file, std::string keyword);

    ///
    /// \brief GetRigidBodyModes
    /// \return mRigidBodyModes
    ///
    Matrix&                      GetRigidBodyModes()           {return mRigidBodyModes;}

    ///
    /// \brief GetInterfaceRigidBodyModes
    /// \return mInterfaceRigidBodyModes
    ///
    Matrix&                      GetInterfaceRigidBodyModes()  {return mInterfaceRigidBodyModes;}

    ///
    /// \brief GetConnectivityMatrix
    /// \return mConnectivityMatrix
    ///
    SparseMatrix&                GetConnectivityMatrix()       {return mConnectivityMatrix;}

    ///
    /// \brief AssembleConnectivityMatrix
    ///
    ///
    void AssembleConnectivityMatrix();

    ///
    /// \brief AssembleBoundaryDofIds
    ///
    void AssembleBoundaryDofIds();

    ///
    /// \brief mRigidBodyModes
    ///
    /// \f[
    /// R_s \in \mathbb{R}^{\text{number of degrees of freedom} \times \text{number of rigid body modes}  }
    /// \f]
    ///
    Matrix                      mRigidBodyModes;


    ///
    /// \brief mInterfaceRigidBodyModes
    ///
    /// \f[
    /// B_s R_s \in \mathbb{R}^{\text{number of Lagrange multipliers} \times \text{number of rigid body modes}  }
    /// \f]
    ///
    Matrix                      mInterfaceRigidBodyModes;

    ///
    /// \brief mConnectivityMatrix
    ///
    /// \f[
    /// B_s \in \mathbb{R}^{\text{number of Lagrange multipliers} \times \text{number of degrees of freedom}  }
    /// \f]
    ///
    SparseMatrix                mConnectivityMatrix;



    ///
    /// \brief mG Natural coarse space
    ///
    /// \f[
    /// G =
    /// \begin{bmatrix}
    /// B^T_1 R_1
    /// &
    /// B^T_2 R_2
    /// &
    /// \dots
    /// &
    /// B^T_{N_s} R_{N_s}
    /// \end{bmatrix}
    /// \f]
    ///
    /// \f[
    /// G \in \mathbb{R}^{\text{number of Lagrange multipliers} \times \text{total number of rigid body modes}  }
    /// \f]
    Matrix                      mG;


    ///
    /// \brief mProjectionMatrix
    ///
    /// \f[
    /// P = I - G (G^T G)^{-1} G^T
    /// \f]
    ///
    /// \f[
    /// P \in \mathbb{R}^{\text{number of Lagrange multipliers} \times \text{number of Lagrange multipliers}  }
    /// \f]
    ///
    Matrix                      mProjectionMatrix;


    NodeList                    mNodes;
    ElementList                 mElements;
    BoundaryList                mBoundaries;
    InterfaceList               mInterfaces;
    int                         mNumRigidBodyModes = 0;
    int                         mNumInterfaceNodesTotal;
    std::set<int>               mSubdomainBoundaryNodeIds;

    std::vector<int>            mBoundaryDofIds;
    std::vector<int>            mGlobalBoundaryDofIds;
    int                         mNumTotalBoundaryDofIds;
    int                         mGlobalStartIndexBoundaryDofIds;

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
    bool IsFloating() { return mIsFloating;}

    ///
    /// \brief SetIsFloating
    /// \param isFloating
    ///
    void SetIsFloating(const bool isFloating) {mIsFloating = isFloating;}

    ///
    /// \brief CalculateRigidBodyModes
    ///
    void CalculateRigidBodyModes();

protected:


    bool mIsFloating = true;




};
} //namespace NuTo

