//
// Created by Kevin Gori on 07/04/2021.
//

#pragma once

#include "SubstitutionModels.h"

namespace Arbor {

class GTR : public ReversibleModel {
public:
    GTR()
    : m_AtoC(1.0), m_AtoG(1.0), m_AtoT(1.0), m_CtoG(1.0), m_CtoT(1.0),
    m_A(0.25), m_C(0.25), m_G(0.25), m_T(0.25) {
        initialise();
    }

    GTR(double AtoC, double AtoG, double AtoT, double CtoG, double CtoT, double a, double c, double g, double t)
    : m_AtoC(AtoC), m_AtoG(AtoG), m_AtoT(AtoT), m_CtoG(CtoG), m_CtoT(CtoT),
      m_A(a), m_C(c), m_G(g), m_T(t) {
        initialise();
    }

    inline const Matrix4 &Q() const override { return m_Q; };

    const Matrix4 P(double time) const override;

private:
    double m_AtoC, m_AtoG, m_AtoT, m_CtoG, m_CtoT, m_A, m_C, m_G, m_T;
    Matrix4 m_Q, m_Eigenvectors, m_InverseEigenvectors;
    Vector4 m_Eigenvalues;
    void initialise();
};

TEST_CASE("Testing GTR") {
    GTR gtr(2.0, 3.0, 4.0, 5.0, 6.0, 0.1, 0.2, 0.3, 0.4);

    SUBCASE("QMatrix") {
        Matrix4 q = gtr.Q();
        REQUIRE(q[0] == doctest::Approx(-1.21848739));
        REQUIRE(q[1] == doctest::Approx(0.08403361));
        REQUIRE(((q[0] + q[4] + q[8] + q[12]) == doctest::Approx(0.0)));
    }

    SUBCASE("PMatrix") {
        Matrix4 p = gtr.P(1.08);
        REQUIRE(p[0] == doctest::Approx(0.30969459));
        REQUIRE((p[0] + p[4] + p[8] + p[12]) == doctest::Approx(1.0));
    }
}
}