#pragma once

#include <iostream>


#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/QR>


namespace NuTo
{

/// \brief Generalized minimal residual method
template <class T, class Preconditioner = Eigen::DiagonalPreconditioner<double>>
int Gmres(const T& A, const Eigen::VectorXd& rhs, Eigen::VectorXd& x, const int maxNumRestarts, const double tolerance,
          const int krylovDimension)
{

    using MatrixType = Eigen::MatrixXd;
    using VectorType = Eigen::VectorXd;

    const int n = A.rows();

    // Initialize upper Hessenberg matrix
    MatrixType H = MatrixType::Zero(krylovDimension + 1, krylovDimension);

    // Initialize Krylov subspace
    MatrixType V = MatrixType::Zero(n, krylovDimension + 1);

    // container for Givens rotation matrices
    std::vector<Eigen::JacobiRotation<double>> Givens(krylovDimension + 1);

    VectorType e1 = VectorType::Zero(n);
    e1[0] = 1.0;

    // Initialize preconditioner
    Preconditioner precond(A);

    // preconditioned residual
    VectorType r = precond.solve(rhs - A * x);
    double rNorm = r.norm();
    double rhsNorm = rhs.norm();

    if (rhsNorm < 1.e-5)
        rhsNorm = 1.0;

    int numRestarts = 0;


    while (rNorm > tolerance * rhsNorm and numRestarts < maxNumRestarts)
    {

        // the first vector in the Krylov subspace is the normalized residual
        V.col(0) = r / rNorm;

        // initialize the s vector used to estimate the residual
        VectorType s = rNorm * e1;

        // construct orthonormal basis using Gram-Schmidt
        for (int i = 0; i < krylovDimension; ++i)
        {
            VectorType w = precond.solve(A * V.col(i));

            for (int iRow = 0; iRow < i + 1; ++iRow)
            {
                H(iRow, i) = w.transpose() * V.col(iRow);
                w = w - H(iRow, i) * V.col(iRow);
            }

            H(i + 1, i) = w.norm();
            V.col(i + 1) = w / H(i + 1, i);

            // Apply the Givens Rotations to ensure that H is an upper triangular matrix.
            // First apply previous rotations to the current matrix
            for (int iRow = 0; iRow < i; ++iRow)
            {
                H.col(i).applyOnTheLeft(iRow, iRow + 1, Givens[iRow].adjoint());
            }

            // form the i-th rotation matrix
            Givens[i].makeGivens(H(i, i), H(i + 1, i));

            // Apply the new Givens rotation on the
            // new entry in the uppper Hessenberg matrix
            H.col(i).applyOnTheLeft(i, i + 1, Givens[i].adjoint());

            // Finally apply the new Givens rotation on the s vector
            s.applyOnTheLeft(i, i + 1, Givens[i].adjoint());

            rNorm = std::abs(s[i + 1]);

            if (rNorm < tolerance * rhsNorm)
            {
                VectorType y = s.head(i + 1);
                H.topLeftCorner(i + 1, i + 1).triangularView<Eigen::Upper>().solveInPlace(y);
                x = x + V.block(0, 0, n, i + 1) * y;
                return numRestarts;
            }
        }


        // we have exceeded the number of iterations. Update the approximation and start over
        VectorType y = s.head(krylovDimension);
        H.topLeftCorner(krylovDimension, krylovDimension).triangularView<Eigen::Upper>().solveInPlace(y);
        x = x + V.block(0, 0, n, krylovDimension) * y;

        // compute preconditioned residual
        r = precond.solve(rhs - A * x);
        rNorm = r.norm();

        ++numRestarts;
    }

    return numRestarts;
}

} // namespace NuTo
