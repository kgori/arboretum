//
// Created by Kevin Gori on 07/04/2021.
//

#pragma once

#include "TN93.h"
#include "doctest.h"

namespace Arbor {
class F84 : public TN93 {
public:
    F84()
        : TN93() {}

    F84(double kappa, double a, double c, double g, double t)
        : TN93(1.0 + kappa / (c + t), 1.0 + kappa / (a + g), 1.0, a, c, g, t) {}

private:
    double m_Kappa;
};

TEST_CASE("Testing F84") {
    F84 f84(2.0, 0.1, 0.2, 0.3, 0.4);

    SUBCASE("QMatrix") {
        Matrix4 q = f84.Q();
        REQUIRE(q[0] == doctest::Approx(-1.56521739));
        REQUIRE(q[5] == doctest::Approx(-1.39130435));
        REQUIRE(q[10] == doctest::Approx(-0.7826087));
        REQUIRE(q[15] == doctest::Approx(-0.82608696));
        REQUIRE(((q[0] + q[4] + q[8] + q[12]) == doctest::Approx(0.0)));
    }

    SUBCASE("PMatrix") {
        Matrix4 p = f84.P(0.85);
        REQUIRE(p[0] == doctest::Approx(0.3283379));
        REQUIRE((p[0] + p[4] + p[8] + p[12]) == doctest::Approx(1.0));
    }
}

}