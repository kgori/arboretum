// Implements the Tamura-Nei substitution model from 1993.
// Several other models are special cases of the TN93 model - HKY85, F84, F81 and K80.
// The TN93 model has an analytical solution for the eigenvalues and vectors of the
// Q Matrix4, so is cheaper to compute than the GTR model.
// Created by Kevin Gori on 06/04/2021.
//

#pragma once

#include "SubstitutionModels.h"
#include "doctest.h"

namespace Arbor {

// Scaling factor for TN93 Q Matrix
double tn93_scale(double pi_a, double pi_c, double pi_g, double pi_t, double alpha_y, double alpha_r, double beta) {
    return 2 * (alpha_y * pi_c * pi_t +
                beta * pi_a * pi_t +
                beta * pi_a * pi_c +
                alpha_r * pi_a * pi_g +
                beta * pi_g * pi_t +
                beta * pi_c * pi_g);
}

// Constructs a TN93 Q matrix from rate and frequency parameters
Matrix4 tn93_q(double pi_a, double pi_c, double pi_g, double pi_t, double alpha_y, double alpha_r, double beta, double scale) {
    double pi_y = pi_t + pi_c;
    double pi_r = pi_a + pi_g;
    Matrix4 q{
            -(alpha_r * pi_g + beta * pi_y), beta * pi_a, alpha_r * pi_a, beta * pi_a,
            beta * pi_c, -(alpha_y * pi_t + beta * pi_r), beta * pi_c, alpha_y * pi_c,
            alpha_r * pi_g, beta * pi_g, -(alpha_r * pi_a + beta * pi_y), beta * pi_g,
            beta * pi_t, alpha_y * pi_t, beta * pi_t, -(alpha_y * pi_c + beta * pi_r)
    };
    for (auto &elem : q) {
        elem /= scale;
    }
    return q;
}

// Computes to eigenvectors of Q
Matrix4 tn93_evecs(double pi_a, double pi_c, double pi_g, double pi_t) {
    double pi_y = pi_t + pi_c;
    double pi_r = pi_a + pi_g;
    Matrix4 eigenvectors{
            1.0, 1.0, 1.0, 1.0,
            -1.0 / pi_r, 1.0 / pi_y, -1.0 / pi_r, 1.0 / pi_y,
            pi_g / pi_r, 0.0, -pi_a / pi_r, 0.0,
            0.0, -pi_t / pi_y, 0.0, pi_c / pi_y
    };
    return eigenvectors;
}

// Computes the inverse eigenvectors of Q
Matrix4 tn93_ivecs(double pi_a, double pi_c, double pi_g, double pi_t) {
    double pi_y = pi_t + pi_c;
    double pi_r = pi_a + pi_g;
    Matrix4 inverse_eigenvectors{
            pi_a, -pi_a * pi_y, 1.0, 0.0,
            pi_c, pi_c * pi_r, 0.0, -1.0,
            pi_g, -pi_g * pi_y, -1.0, 0.0,
            pi_t, pi_t * pi_r, 0.0, 1.0
    };
    return inverse_eigenvectors;
}

// Computes the eigenvalues of Q
Vector4 tn93_evals(double pi_a, double pi_c, double pi_g, double pi_t, double alpha_y, double alpha_r, double beta,
                   double scale) {
    double pi_y = pi_t + pi_c;
    double pi_r = pi_a + pi_g;
    Vector4 eigenvalues{
            0.0, -beta, -(pi_r * alpha_r + pi_y * beta), -(pi_y * alpha_y + pi_r * beta)
    };
    for (auto &val : eigenvalues) val /= scale;
    return eigenvalues;
};


class TN93 : ReversibleModel {
public:
    TN93()
            : m_AlphaY(1.0), m_AlphaR(1.0), m_Beta(1.0), m_A(0.25), m_C(0.25), m_G(0.25), m_T(0.25) {
        double scale = tn93_scale(m_A, m_C, m_G, m_T, m_AlphaY, m_AlphaR, m_Beta);
        initialise();
    }

    TN93(double alpha_y, double alpha_r, double beta, double a, double c, double g, double t)
            : m_AlphaY(alpha_y), m_AlphaR(alpha_r), m_Beta(beta), m_A(a), m_C(c), m_G(g), m_T(t) {
        initialise();
    }

    inline const Matrix4 &Q() const override { return m_Q; };

    inline const Matrix4 P(double time) const override;

protected:
    void initialise() {
        double scale = tn93_scale(m_A, m_C, m_G, m_T, m_AlphaY, m_AlphaR, m_Beta);
        m_Q = tn93_q(m_A, m_C, m_G, m_T, m_AlphaY, m_AlphaR, m_Beta, scale);
        m_Eigenvectors = tn93_evecs(m_A, m_C, m_G, m_T);
        m_InverseEigenvectors = tn93_ivecs(m_A, m_C, m_G, m_T);
        m_Eigenvalues = tn93_evals(m_A, m_C, m_G, m_T, m_AlphaY, m_AlphaR, m_Beta, scale);
    }
    double m_AlphaY, m_AlphaR, m_Beta, m_A, m_C, m_G, m_T;
    Matrix4 m_Q, m_Eigenvectors, m_InverseEigenvectors;
    Vector4 m_Eigenvalues;
};

inline const Matrix4 TN93::P(double time) const {
    Vector4 tmp{std::exp(m_Eigenvalues[0] * time),
                std::exp(m_Eigenvalues[1] * time),
                std::exp(m_Eigenvalues[2] * time),
                std::exp(m_Eigenvalues[3] * time)};
    return reconstitute_eigen(m_Eigenvectors, tmp, m_InverseEigenvectors);
}


TEST_CASE("Testing TN93") {
    TN93 tn93(10.0, 10.0, 1.0, 0.25, 0.25, 0.25, 0.25);

    SUBCASE("QMatrix") {
        Matrix4 q = tn93.Q();
        REQUIRE(q[0] == doctest::Approx(-1.0));
        REQUIRE(q[1] == doctest::Approx(0.08333333));
        REQUIRE(((q[0] + q[4] + q[8] + q[12]) == doctest::Approx(0.0)));
    }

    SUBCASE("PMatrix") {
        Matrix4 p = tn93.P(2.5);
        REQUIRE(p[0] == doctest::Approx(0.36375994));
        REQUIRE((p[0] + p[4] + p[8] + p[12]) == doctest::Approx(1.0));
    }
}
}