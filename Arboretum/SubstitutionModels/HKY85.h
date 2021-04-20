//
// Created by Kevin Gori on 21/10/2020.
//

#pragma once

#include "TN93.h"
#include "doctest.h"

namespace Arbor {
class HKY85 : public TN93 {
public:
    HKY85()
            : m_Kappa(1.0), TN93() {
        double scale = tn93_scale(m_A, m_C, m_G, m_T, m_AlphaY, m_AlphaR, m_Beta);
        initialise();
    }

    HKY85(double kappa, double a, double c, double g, double t)
            : m_Kappa(kappa), TN93(kappa, kappa, 1.0, a, c, g, t) {
        double scale = tn93_scale(m_A, m_C, m_G, m_T, m_AlphaY, m_AlphaR, m_Beta);
        initialise();
    }

private:
    double m_Kappa;

};

TEST_CASE("Testing HKY85") {
    HKY85 hky85(2.0, 0.1, 0.2, 0.3, 0.4);

    SUBCASE("QMatrix") {
        Matrix4 q = hky85.Q();
        REQUIRE(q[0] == doctest::Approx(-1.30434783));
        REQUIRE(q[1] == doctest::Approx(0.10869565));
        REQUIRE(q[2] == doctest::Approx(0.2173913));
        REQUIRE(q[3] == doctest::Approx(0.10869565));
        REQUIRE(((q[0] + q[4] + q[8] + q[12]) == doctest::Approx(0.0)));
    }

    SUBCASE("PMatrix") {
        Matrix4 p = hky85.P(2.5);
        REQUIRE(p[0] == doctest::Approx(0.12661231));
        REQUIRE((p[0] + p[4] + p[8] + p[12]) == doctest::Approx(1.0));
    }
}
}