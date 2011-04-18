// $Id$

#ifndef STRUCTUREMULTISCALE_H
#define STRUCTUREMULTISCALE_H

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/export.hpp>
#include <boost/serialization/array.hpp>
#else
#include <boost/array.hpp>
#endif //Serialize

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/constitutive/mechanics/ConstitutiveEngineeringStressStrain.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain2D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress2D.h"
#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/elements/ElementDataEnum.h"
#include "nuto/mechanics/elements/IpDataEnum.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"

namespace NuTo
{
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... standard class for irregular (unstructured) structures, which are used as fine scale models at an integration point
class StructureMultiscale : public Structure
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    //! @param mDimension  Structural dimension (1,2 or 3)
    StructureMultiscale(int mDimension);

#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif //SWIG

    //! @brief ... save the object to a file
    //! @param filename ... filename
    //! @param aType ... type of file, either BINARY, XML or TEXT
    void Save (const std::string &filename, std::string rType )const;

    //! @brief ... restore the object from a file
    //! @param filename ... filename
    //! @param aType ... type of file, either BINARY, XML or TEXT
    void Restore (const std::string &filename, std::string rType );
#endif // ENABLE_SERIALIZATION

    //! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
    //!            in case of restoring from a file with the wrong object type, the file id is printed
    //! @return    class name
    std::string GetTypeId()const
    {
        return std::string("StructureMultiscale");
    }

    //************ Info routine         ***************
    //**  defined in structures/StructureMultiscale.cpp *********
    //*************************************************
    //! @brief ... Info routine that prints general information about the object (detail according to verbose level)
    void Info()const;

#ifndef SWIG
    //! @brief ... calculate the displacement based on the homogeneous strain
    //! @param rCoordinates ... coordinates of the point
    //! @param rCoarseDisplacements ... return array of displacements
    void GetDisplacementsEpsilonHom2D(double rCoordinates[2], double rDisplacements[2], const boost::array<double,2>& rCenter)const;

    //! @brief ... calculate the displacement based on the crack opening
    //! @param rCoordinates ... coordinates of the point
    //! @param rCoarseDisplacements ... return array of displacements
    void GetDisplacementsCrack2D(double rCoordinates[2], double rDisplacements[2])const;

    //! @brief derivative of displacement with respect to homogeneous strain
    //! @param rdX_dEpsilonHom[3] return value, derivative of x-displacement with respect to homogeneous strain (exx, eyy, gxy)
    //! @param rdY_dEpsilonHom[3] return value, derivative of x-displacement with respect to homogeneous strain (exx, eyy, gxy)
    void GetdDisplacementdEpsilonHom(double rCoordinates[2], double rdX_dEpsilonHom[3], double rdY_dEpsilonHom[3], const boost::array<double,2>& rCenter)const;

    //! @brief derivative of displacement with respect to discontinuity (crack opening)
    //! @param rdX_dCrackOpening[2] return value, derivative of x-displacement with respect to crack opening (ux, uy)
    //! @param rdY_dCrackOpening[2] return value, derivative of x-displacement with respect to crack opening (ux, uy)
    void GetdDisplacementdCrackOpening(double rCoordinates[2], double rdX_dCrackOpening[2], double rdY_dCrackOpening[2])const;

    //! @brief second derivative of displacement with respect to alpha and discontinuity (crack opening)
    //! @param rdX_dCrackOpening[2] return value, derivative of x-displacement with respect to alpha and crack opening (ux, uy)
    //! @param rdY_dCrackOpening[2] return value, derivative of y-displacement with respect to alpha and crack opening (ux, uy)
    void Getd2Displacementd2CrackOpening(double rCoordinates[2], double rdX_dAlphaCrackOpening[2], double rdY_dAlphaCrackOpening[2])const;

    //! @brief derivative of displacement with respect to discontinuity (crack opening)
    //! @param rdX_dAlpha[2] return value, derivative of x-displacement with respect to crack orientation (alpha)
    //! @param rdy_dAlpha[2] return value, derivative of x-displacement with respect to crack orientation (alpha)
    void GetdDisplacementdCrackOrientation(double rCoordinates[2], double rdX_dAlpha[1], double rdy_dAlpha[1])const;

    //! @brief second derivative of displacement with respect to discontinuity (crack opening)
    //! @param rdX_dAlpha[2] return value, second derivative of x-displacement with respect to crack orientation (alpha)
    //! @param rdy_dAlpha[2] return value, second derivative of y-displacement with respect to crack orientation (alpha)
    void Getd2Displacementd2CrackOrientation(double rCoordinates[2], double rd2X_d2Alpha[1], double rd2Y_d2Alpha[1])const;
#endif //SWIG

    //! @brief ... the boundary nodes were transformed from pure displacement type nodes to multiscale nodes
    //! the displacements are decomposed into a local displacement field and a global homogeneous/crack displacement
    void TransformMultiscaleNodes(int rGroupBoundaryNodes, int rGroupNodes, int rGroupElements, boost::array<double,2>& rCenter, double& rLength, double& rArea, bool rCrackedDomain);

    //! @brief ... the boundary nodes were transformed from pure displacement type nodes to multiscale nodes
    //! the displacements are decomposed into a local displacement field and a global homogeneous/crack displacement
    void TransformMultiscaleNodes();

#ifndef SWIG
    //! @brief numbers non standard DOFs' e.g. in StructureMultiscale, for standard structures this routine is empty
    void NumberAdditionalGlobalDofs();

    //! @brief ...merge additional dof values
    void NodeMergeAdditionalGlobalDofValues(const NuTo::FullMatrix<double>& rActiveDofValues, const NuTo::FullMatrix<double>& rDependentDofValues);

    //! @brief extract dof values additional dof values
    //! @param rActiveDofValues ... vector of global active dof values (ordering according to global dofs, size is number of active dofs)
    //! @param rDependentDofValues ... vector of global dependent dof values (ordering according to (global dofs) - (number of active dofs), size is (total number of dofs) - (number of active dofs))
    void NodeExtractAdditionalGlobalDofValues(NuTo::FullMatrix<double>& rActiveDofValues, NuTo::FullMatrix<double>& rDependentDofValues) const;

    //! @brief ... based on the global dofs build submatrices of the global coefficent matrix0
    //! @param rMatrixJJ ... submatrix jj (number of active dof x number of active dof)
    //! @param rMatrixJK ... submatrix jk (number of active dof x number of dependent dof)
    void BuildGlobalCoefficientSubMatrices0General(NuTo::SparseMatrix<double>& rMatrixJJ, NuTo::SparseMatrix<double>& rMatrixJK) const;

    //! @brief ... based on the global dofs build submatrices of the global coefficent matrix0
    //! @param rMatrixJJ ... submatrix jj (number of active dof x number of active dof)
    //! @param rMatrixJK ... submatrix jk (number of active dof x number of dependent dof)
    //! @param rMatrixKJ ... submatrix kj (number of dependent dof x number of active dof)
    //! @param rMatrixKK ... submatrix kk (number of dependent dof x number of dependent dof)
    void BuildGlobalCoefficientSubMatrices0General(NuTo::SparseMatrix<double>& rMatrixJJ, NuTo::SparseMatrix<double>& rMatrixJK, NuTo::SparseMatrix<double>& rMatrixKJ, NuTo::SparseMatrix<double>& rMatrixKK) const;

    //! @brief ... based on the global dofs build submatrices of the global coefficent matrix0
    //! @param rMatrixJJ ... submatrix jj (number of active dof x number of active dof)
    //! @param rMatrixJK ... submatrix jk (number of active dof x number of dependent dof)
    void BuildGlobalCoefficientSubMatrices0Symmetric(NuTo::SparseMatrix<double>& rMatrixJJ, NuTo::SparseMatrix<double>& rMatrixJK) const;

    //! @brief ... based on the global dofs build submatrices of the global coefficent matrix0
    //! @param rMatrixJJ ... submatrix jj (number of active dof x number of active dof)
    //! @param rMatrixJK ... submatrix jk (number of active dof x number of dependent dof)
    //! @param rMatrixKK ... submatrix kk (number of dependent dof x number of dependent dof)
    void BuildGlobalCoefficientSubMatrices0Symmetric(NuTo::SparseMatrix<double>& rMatrixJJ, NuTo::SparseMatrix<double>& rMatrixJK, NuTo::SparseMatrix<double>& rMatrixKK) const;

    //! @brief ... for a single element matrix, assemble the global matrix, consider here additionally the global degrees of freedom from StructureMultiscale (alpha, ux,uy) and eventually the macroscopic strain
    //! @param elementMatrix ...stiffness matrix of the element
    //! @param elementMatrixGlobalDofsRow ...corresponding row DOFs
    //! @param elementMatrixGlobalDofsColumn ...corresponding column DOFs
    //! @param elementVector ...gradient matrix of the element (internal force)
    //! @param mappingDofMultiscaleNode ...mapping of all Dofs to Multiscale dofs (either -1 for standard dofs or the entry in dDOF and dDOF2, where additional infos are stored)
    //! @param dDOF ... derivatives of multiscale dofs with respect to alpha, ux, uy, ehomxx, ehomyy, gammahomxy [0..6] scaled derivatives [7..12]
    //! @param dDOF2 ...for each dof, the corresponding second order derivative (alpha^2, alpha ux, alpha uy)
    //! @param rMatrixJJ ... submatrix jj (number of active dof x number of active dof)
    //! @param rMatrixJK ... submatrix jk (number of active dof x number of dependent dof)
    //! @param rMatrixKJ ... submatrix kj (number of dependent dof x number of active dof)
    //! @param rMatrixKK ... submatrix kk (number of dependent dof x number of dependent dof)
    //! @param rCalcMatrixKJ_KK ...true, if rMatrixKJ and rMatrixKK should be calculated (rMatrixJJ and rMatrixJK) are always calculated
    //! @param rCalcMatrixM ...true, if all matrices corresponding to M should be calculated
    void AddElementMatrixToGlobalSubMatricesGeneral(
            NuTo::FullMatrix<double>& elementMatrix,
            std::vector<int>& elementMatrixGlobalDofsRow,
            std::vector<int>& elementMatrixGlobalDofsColumn,
            NuTo::FullMatrix<double>& elementVector,
            std::vector<int>& mappingDofMultiscaleNode,
            std::vector<boost::array<double,6> >& rDOF,
            std::vector<boost::array<double,3> >& rDOF2,
            NuTo::SparseMatrix<double>* rMatrixJJ,
            NuTo::SparseMatrix<double>* rMatrixJK,
            NuTo::SparseMatrix<double>* rMatrixKJ,
            NuTo::SparseMatrix<double>* rMatrixKK,
            bool rCalcMatrixKJ_KK)const;

    //! @brief ... based on the global dofs build sub-vectors of the global internal potential gradient
    //! @param rActiveDofGradientVector ... global internal potential gradient which corresponds to the active dofs
    //! @param rDependentDofGradientVector ... global internal potential gradient which corresponds to the dependent dofs
    void BuildGlobalGradientInternalPotentialSubVectors(NuTo::FullMatrix<double>& rActiveDofGradientVector, NuTo::FullMatrix<double>& rDependentDofGradientVector) const;

    //! @brief ... based on the global dofs build sub-vectors of the global internal potential gradient
    //! @param rActiveDofGradientVector ... global internal potential gradient which corresponds to the active dofs
    //! @param rDependentDofGradientVector ... global internal potential gradient which corresponds to the dependent dofs
    //! @param rEpsilonMGradientVector ... global internal potential gradient which corresponds to the global macroscopic strain (derivative is equal to homogeneous strain)
    void BuildGlobalGradientInternalPotentialSubVectors(NuTo::FullMatrix<double>& rActiveDofGradientVector, NuTo::FullMatrix<double>& rDependentDofGradientVector, NuTo::FullMatrix<double>* rEpsilonMGradientVector) const;

    void CalculateCohesiveForce(NuTo::FullMatrix<double>& rCohesiveForce) const;

    //! @brief calculate the distance of a point to the crack
    //! @param rCoordinates[2] coordinates of the point
    //! @return distance to crack
    double CalculateDistanceToCrack2D(double rCoordinates[2])const;

    //! @brief calculate the derivative of the distance of a point to the crack with respect to the crack orientation alpha
    //! @param rCoordinates[2] coordinates of the point
    //! @return distance to crack
    double CalculatedDistanceToCrack2DdAlpha(double rCoordinates[2])const;

    //! @brief calculate the second derivative of the distance of a point to the crack
    //! @param rCoordinates[2] coordinates of the point
    //! @return second derivative of distance to crack
    double Calculated2DistanceToCrack2Dd2Alpha(double rCoordinates[2])const;

    //! @brief calculate the derivative of the displacements at the nodes with respect to homogeneous strain, crack opening and crack orientation
    //! @param rMappingDofMultiscaleNode return value, for each dof, the corresponding entry in the rDOF vector, for nonmultiscale dofs, there is a -1
    //! @param rDOF return value, for each dof, the corresponding derivatives ehomxx, ehomyy, gammahomxy, ux, uy, alpha  [0..5]
    //! @param rDOF2 return value, for each dof, the corresponding second derivatives (alpha^2, alpha ux, alpha uy)
    void CalculatedDispdGlobalDofs(std::vector<int>& rMappingDofMultiscaleNode, std::vector<boost::array<double,6> >& rDOF, std::vector<boost::array<double,3> >& rDOF2)const;

    //! @brief calculate the derivative of the displacements at the nodes with respect to crack opening without considering this to be balanced by epsilon_hom
    //! @param rMappingDofMultiscaleNode return value, for each dof, the corresponding entry in the rDOF vector, for nonmultiscale dofs, there is a -1
    //! @param rDOF return value, for each dof, the corresponding derivatives (ux, uy)
    void CalculatedDispdCrackOpening(std::vector<int>& rMappingDofMultiscaleNode,std::vector<boost::array<double,2> >& rDOF)const;

    //! @brief set the total current strain
    void SetTotalEngineeringStrain(EngineeringStrain2D& rTotalEngineeringStrain);

    //! @brief set the maximum total strain (used in the automatic increment of the Newton Raphson iteration multiplied by the load factor to give the totalEngineeringStrain)
    void SetDeltaTotalEngineeringStrain(EngineeringStrain2D& rDeltaTotalEngineeringStrain);

    //! @brief set the maximum total strain (used in the automatic increment of the Newton Raphson iteration multiplied by the load factor to give the totalEngineeringStrain)
    void SetPrevTotalEngineeringStrain(EngineeringStrain2D& rPrevTotalEngineeringStrain);

    //! @brief returns the total strain
    const NuTo::EngineeringStrain2D& GetTotalEngineeringStrain()const;

    //! @brief just for testing
    void SetGlobalCrackOpening(boost::array<double,2> ptr)const
    {
        const_cast<StructureMultiscale*>(this)->mCrackOpening[0] = ptr[0];
        const_cast<StructureMultiscale*>(this)->mCrackOpening[1] = ptr[1];
    }

    //! @brief just for testing
    void SetCrackAngle(double alpha)const
    {
        const_cast<StructureMultiscale*>(this)->mCrackAngle = alpha;
    }

    //! @brief return the previous crack angle
    double GetPrevCrackAngle()const;

    //! @brief sets the previous crack angle
    void SetPrevCrackAngle(double rPrevCrackAngle);

    //! @brief return the total strain
    const NuTo::EngineeringStrain2D& GetHomogeneousEngineeringStrain()const;

    //! @brief set the homgeneous strain (just for test purpose)
    void SetHomogeneousEngineeringStrain(NuTo::EngineeringStrain2D rStrain);

    //! @brief renumbers the global dofs in the structure after
    void ReNumberAdditionalGlobalDofs(std::vector<int>& rMappingInitialToNewOrdering);


#endif
    //! @brief returns the constraint equation for the crack angle
    int GetConstraintCrackAngle()const
    {
        return mConstraintCrackAngle;
    }

    //! @brief returns the constraint equation for the total strain
    int GetConstraintTotalStrain()const
    {
        return mConstraintTotalStrain;
    }

    void SetlCoarseScale(double rlCoarseScale)
    {
        mlCoarseScale = rlCoarseScale;
    }
    void SetSquareCoarseScaleModel(bool rSquareCoarseScaleModel)
    {
        mSquareCoarseScaleModel = rSquareCoarseScaleModel;
    }

    double GetlCoarseScale()const
    {
        return mlCoarseScale;
    }

    double GetAreaDamage()const
    {
        return mFineScaleAreaDamage;
    }

    double GetAreaHomogeneous()const
    {
        return mFineScaleAreaHomogeneous;
    }

    int GetDofCrackAngle()const
    {
        return mDOFCrackAngle;
    }

    double GetCrackAngle()const
    {
        return mCrackAngle;
    }

    boost::array<int,2> GetDofGlobalCrackOpening2D()const
    {
        return mDOFCrackOpening;
    }

    boost::array<double,2> GetGlobalCrackOpening2D()const
    {
        return mCrackOpening;
    }
    boost::array<int,3> GetDofGlobalTotalStrain2D()const
    {
        return mDOFGlobalTotalStrain;
    }

    void SetGroupBoundaryNodesElements(int rGroupIdBoundaryNodesDamage, int rGroupIdBoundaryNodesHomogeneous, int rGroupIdNodesDamage, int rGroupIdNodesHomogeneous, int rGroupIdElementsDamage, int rGroupIdElementsHomogeneous)
    {
        mGroupBoundaryNodesDamage = rGroupIdBoundaryNodesDamage;
        mGroupBoundaryNodesHomogeneous = rGroupIdBoundaryNodesHomogeneous;
        mGroupNodesDamage = rGroupIdNodesDamage;
        mGroupNodesHomogeneous = rGroupIdNodesHomogeneous;
        mGroupElementsDamage = rGroupIdElementsDamage;
        mGroupElementsHomogeneous = rGroupIdElementsHomogeneous;
        mBoundaryNodesElementsAssigned = true;
    }

    int GetGroupElementsDamage()const
    {
    	return mGroupElementsDamage;
    }

    int GetGroupElementsHomogeneous()const
    {
    	return mGroupElementsHomogeneous;
    }

    int GetGroupBoundaryNodesDamage()const
    {
    	return mGroupBoundaryNodesDamage;
    }

    int GetGroupBoundaryNodesHomogeneous()const
    {
    	return mGroupBoundaryNodesHomogeneous;
    }

    double GetScalingFactorDamage()const
    {
        return mlCoarseScale/mlFineScaleDamage;
    }

    double GetScalingFactorHomogeneous()const
    {
    	return mlCoarseScale*(mlCoarseScale-mFineScaleAreaDamage/mlFineScaleDamage)/(mFineScaleAreaHomogeneous);
    }

    const std::string& GetIPName()const
    {
        return mIPName;
    }
    //calculate from the existing crack opening and orientation the cracking strain and the homogeneous strain
    //for test purpose in public section
    void CalculateHomogeneousEngineeringStrain();

    void SetCrackTransitionRadius(double rCrackTransitionRadius)
    {
        mCrackTransitionRadius = rCrackTransitionRadius;
    }

    //! @brief calculates the crack angle for elastic solutions
    double CalculateCrackAngleElastic()const;

    //! @brief calculates the difference between the crack angle of the elastic solution and the current angle
    //! attention, the periodicity of the crack angle has to be taken into account
    double CalculateDeltaCrackAngleElastic()const;

    //! @brief calculates the crack angle for elastic solutions (initial value, no scaling with previous crack angle)
    //! @return alpha crack angle in the range 0..Pi (return value)
    double CalculateInitialCrackAngleElastic()const;

    //! @brief add a constraint equation for alpha, which corresponds to an artificial spring
    //! @parameter rPenaltyStiffness penalty stiffness
    //! @parameter rScalingFactor scaling factor
    //! @return id of the constraint
    int ConstraintNonlinearCrackAngle(double rPenaltyStiffness, double rScalingFactor);

    //! @brief set the penalty stiffness for the nonlinear crack angle constraint
    void SetPenaltyStiffnessCrackAngle(double rParameter);

    //! @brief set the scaling factor for the nonlinear crack angle constraint
    void SetPenaltyStiffnessScalingFactorCrackAngle(double rParameter);

    //! @brief add a constraint equation for the tangential crack opening, which corresponds to an artificial spring
    //! @parameter rScalingFactor scaling factor
    //! @parameter rPenaltyStiffness penalty stiffness
    //! @return id of the constraint
    int ConstraintNonlinearTangentialCrackOpening(double rScalingFactor, double rPenaltyStiffness);

    //! @brief set the penalty stiffness for the nonlinear TangentialCrackOpening constraint
    void SetPenaltyStiffnessTangentialCrackOpening(double rParameter);

    //! @brief set the scaling factor for the nonlinear TangentialCrackOpening constraint
    void SetPenaltyStiffnessScalingFactorTangentialCrackOpening(double rParameter);

    //! @brief set the tolerance for the transition between crack angle from principal strain and previous strain
    void SetToleranceElasticCrackAngleHigh(double rParameter);

    //! @brief set the tolerance for the transition between crack angle from principal strain and previous strain
    void SetToleranceElasticCrackAngleLow(double rParameter);

    //! @brief add a constraint equation for the crack opening (normal crack opening non negativ)
    //! @parameter rPenaltyStiffness penalty stiffness for augmented Lagrangian
    //! @return id of the constraint
    int ConstraintLagrangeCrackOpening(double rPenaltyStiffness);

    //! @brief add a constraint equation for the total strain
    //! @parameter rStrain applied strain (rhs)
    //! @return id of the constraint
    int ConstraintLinearGlobalTotalStrain(const EngineeringStrain2D& rStrain);

    //! @brief this routine is only relevant for the multiscale model, since an update on the fine scale should only be performed
    //for an update on the coarse scale
    //as a consequence, in an iterative solution with updates in between the initial state has to be restored after leaving the routine
    //this routine saves the current state before an update in the Newton Raphson iteration is performed
    //this only happens for more than one load step (either prescibed or with automatic load control)
    //if not alreaday saved, the structure is saved to an internal stream (be careful to not include this stream in the serialization routine)
    void SaveStructure()const;

    void RestoreStructure();

    //! @brief set the load factor (load or displacement control) overload this function to use Newton Raphson
    //! @param load factor
    void SetLoadFactor(double rLoadFactor);

    //! @brief do a postprocessing step after each converged load step (for Newton Raphson iteration) overload this function to use Newton Raphson
    void PostProcessDataAfterConvergence(int rLoadStep, int rNumNewtonIterations, double rLoadFactor, double rDeltaLoadFactor)const;

    //! @brief do a postprocessing step after each line search within the load step(for Newton Raphson iteration) overload this function to use Newton Raphson
    void PostProcessDataAfterLineSearch(int rLoadStep, int rNewtonIteration, double rLineSearchFactor, double rLoadFactor)const;

    //! @brief only for debugging, info at some stages of the Newton Raphson iteration
    void NewtonRaphsonInfo(int rVerboseLevel)const;

    //! @brief initializes some variables etc. before the Newton-Raphson routine is executed
    void InitBeforeNewtonRaphson();

    bool GetSavedToStringStream()
    {
        return mSavedToStringStream;
    }

    double GetlFineScaleDamage()const
    {
        return mlFineScaleDamage;
    }

    double GetlFineScaleHomogeneous()const
    {
        return mlFineScaleHomogeneous;
    }

    const boost::array<double,2>& GetCenterDamage()const
    {
        return mCenterDamage;
    }

    const boost::array<double,2>& GetCenterHomogeneous()const
    {
        return mCenterHomogeneous;
    }

    //! @brief calculates the total energy of the system
    //! @return total energy
    double ElementTotalGetTotalEnergy()const;

    //! @brief calculates the total energy of the system
    //! @return total energy
    double ElementGroupGetTotalEnergy(int rGroupId)const;

    //! @brief sets the result directory where the results are written to
    void SetResultDirectoryAfterLineSearch(std::string rResultDirectoryAfterLineSearch)
    {
    	mResultDirectoryAfterLineSearch = rResultDirectoryAfterLineSearch;
    }

    //! @brief sets the result file for the converged solution where the results are written to
    void SetResultFileAfterConvergence(std::string rResultFileAfterConvergence)
    {
    	mResultFileAfterConvergence = rResultFileAfterConvergence;
    }

protected:
    //! @brief ... standard constructor just for the serialization routine
    StructureMultiscale()
    {}

    //! @brief Calculate the derivate of the homogeneous strain with respect to changes of the crack orientation and crack opening
    //! this is due to the constraint equation relating total strain, homogeneous strain and cracking strain
    //! @parameter rbHomAlpha dHom wrt alpha
    //! @paramter rbHomU[0-2] for wrt ux [3-5] for wrt uy
    //! @parameter bHessian depsilondalpha2[0-2], depsilondalphadux[3-5], depsilondalphadux[6-8]
    void GetdEpsilonHomdCrack(double rbHomAlpha[3], double rbHomU[6], double rbHessian[9])const;


    double mCrackAngle;
    double mPrevCrackAngle;
    int mDOFCrackAngle;
    boost::array<double,2> mCrackOpening; //UT UN
    boost::array<int,2> mDOFCrackOpening;
    //! @brief this is the current total strain
    EngineeringStrain2D mEpsilonTot;
    boost::array<int,3> mDOFGlobalTotalStrain;
    //! @brief this is the current homogeneous part of the strain
    EngineeringStrain2D mEpsilonHom;
    //! @brief this is the previous homogeneous part of the strain (last update)
    EngineeringStrain2D mPrevEpsilonHom;
    //! @brief this is the previous average stress (last update)
    EngineeringStress2D mPrevAverageStress;
    //! @brief true for an approximately square coarse scale model, false for a either triangular meshes or meshes not aligned with the axes
    //! in case the coarse scale element is well aligned, the ratio CrackLength/Area can be calculated exactly
    //! whereas for other shapes, an approximation has to be used (additional errors)
    bool mSquareCoarseScaleModel;
    //! @brief for mSquareCoarseScaleModel == edge length of the macroscale element (sqrt of the area of the macroscale element)
    //! for mSquareCoarseScaleModel==false ratio of the crack length divided by the area approximated by (sqrt of the area of the macroscale element)
    double mlCoarseScale;
    //! @brief parameter length of the fine scale model related to the damage model
    double mlFineScaleDamage;
    //! @brief parameter length of the fine scale model related to the homogeneous model
    double mlFineScaleHomogeneous;
    //! @brief length to regularize the Heaviside function
    double mCrackTransitionRadius;
    //! @brief center of the model
    boost::array<double,2> mCenterDamage;
    //! @brief center of the model
    boost::array<double,2> mCenterHomogeneous;
    //! @brief area of the fine scale model
    // this has to be stored, since for a round fine scale model (especially coarse scale) the real area is substantially
    // lower than the area of the circle
    double mFineScaleAreaDamage;
    double mFineScaleAreaHomogeneous;
    //! @brief thickness of the structure
    double mThickness;
    int mConstraintFineScaleX,
        mConstraintFineScaleY,
        mConstraintNormalCrackOpening,
        mConstraintTangentialCrackOpening,
        mConstraintCrackAngle,
        mConstraintTotalStrain;
    bool mBoundaryNodesElementsAssigned;
    int mGroupBoundaryNodesDamage;
    int mGroupBoundaryNodesHomogeneous;
    int mGroupNodesDamage;
    int mGroupNodesHomogeneous;
    int mGroupElementsDamage;
    int mGroupElementsHomogeneous;
    std::string mIPName;
    std::string mResultDirectoryAfterLineSearch;
    std::string mResultFileAfterConvergence;

    //! @brief tolerance for difference between the principal strains, where
    //if difference is bigger, the angle is calculated from the largest principal strain
    double mToleranceElasticCrackAngleHigh;
    //! @brief tolerance for difference between the principal strains, where
    //the previous angle is used, no update, (in between, there is a linear interpolation)
    double mToleranceElasticCrackAngleLow;

    // the following parameters are used for the Newton-Raphson iteration
    //! @brief this is the previous strain (used in the load application of the Newton procedure)
    EngineeringStrain2D mPrevEpsilonTot;
    //! @brief this is the delta strain (used in the load application of the Newton procedure)
    EngineeringStrain2D mDeltaEpsilonTot;
    //! @brief set to true, if during an update in NR the structure had to be saved, otherwise unchanged
    mutable bool mSavedToStringStream;
    //! @brief structure is saved to that string stream
    mutable std::stringstream mSaveStringStream;
};
}
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
BOOST_CLASS_EXPORT_KEY(NuTo::StructureMultiscale)
#endif // SWIG
#endif // ENABLE_SERIALIZATION
#endif // StructureMultiscale_H
