//
// Created by Kevin Gori on 08/04/2021.
//

#include <algorithm>
#include "GTR.h"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

namespace Arbor {
Matrix4 createQ(double AtoC, double AtoG, double AtoT, double CtoG, double CtoT, double GtoT, double piA, double piC, double piG, double piT) {
    Matrix4 q{
            -(AtoC * piC + AtoG * piG + AtoT * piT), AtoC * piA, AtoG * piA, AtoT * piA,
            AtoC * piC, -(AtoC * piA + CtoG * piG + CtoT * piT), CtoG * piC, CtoT * piC,
            AtoG * piG, CtoG * piG, -(AtoG * piA + CtoG * piC + GtoT * piT), GtoT * piG,
            AtoT * piT, CtoT * piT, GtoT * piT, -(AtoT * piA + CtoT * piC + GtoT * piG)
    };

    double scale = -1.0 / (q[0] * piA + q[5] * piC + q[10] * piG + q[15] * piT);
    std::transform(q.begin(), q.end(), q.begin(), [scale](double val) -> double { return val * scale; });
    return q;
}

void GTR::initialise() {
    m_Q = createQ(m_AtoC, m_AtoG, m_AtoT, m_CtoG, m_CtoT, 1.0, m_A, m_C, m_G, m_T);
    Eigen::Map<const Eigen::Matrix4d> Q(m_Q.data());
    Eigen::Vector4d pi(m_A, m_C, m_G, m_T);
    Eigen::Matrix4d B = pi.array().sqrt().matrix().asDiagonal() * Q * (1.0 / pi.array().sqrt()).matrix().asDiagonal();
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix4d> eigensolver(B);

    auto R = eigensolver.eigenvectors();
    auto evals = eigensolver.eigenvalues();
    auto evecs = (1.0 / pi.array().sqrt()).matrix().asDiagonal() * R;
    auto ivecs = R.transpose() * pi.array().sqrt().matrix().asDiagonal();

    Eigen::Matrix4d::Map(m_Eigenvectors.data()) = evecs;
    Eigen::Matrix4d::Map(m_InverseEigenvectors.data()) = ivecs;
    Eigen::Vector4d::Map(m_Eigenvalues.data()) = evals;
}

const Matrix4 GTR::P(double time) const {
    Vector4 tmp{std::exp(m_Eigenvalues[0] * time),
                std::exp(m_Eigenvalues[1] * time),
                std::exp(m_Eigenvalues[2] * time),
                std::exp(m_Eigenvalues[3] * time)};
    return reconstitute_eigen(m_Eigenvectors, tmp, m_InverseEigenvectors);
}
}
