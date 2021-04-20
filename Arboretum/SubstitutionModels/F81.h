//
// Created by Kevin Gori on 07/04/2021.
//

#pragma once

#include "TN93.h"
#include "doctest.h"

namespace Arbor {
class F81 : public TN93 {
public:
    F81()
        : TN93() {}

    F81(double a, double c, double g, double t)
        : TN93(1.0, 1.0, 1.0, a, c, g, t) {}
};

TEST_CASE("Testing F81") {
    F81 f81(0.1, 0.2, 0.3, 0.4);

    SUBCASE("QMatrix") {
        Matrix4 q = f81.Q();
        REQUIRE(q[0] == doctest::Approx(-1.28571429));
        REQUIRE(q[5] == doctest::Approx(-1.14285714));
        REQUIRE(q[10] == doctest::Approx(-1.0));
        REQUIRE(q[15] == doctest::Approx(-0.85714286));
        REQUIRE(((q[0] + q[4] + q[8] + q[12]) == doctest::Approx(0.0)));
    }

    SUBCASE("PMatrix") {
        Matrix4 p = f81.P(0.15);
        REQUIRE(p[0] == doctest::Approx(0.82640597));
        REQUIRE((p[0] + p[4] + p[8] + p[12]) == doctest::Approx(1.0));
    }
}

}