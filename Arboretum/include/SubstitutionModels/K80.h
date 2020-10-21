//
// Created by Kevin Gori on 21/10/2020.
//

#pragma once

#include "SubstitutionModels.h"
#include "../../vendor/doctest.h"

namespace Arbor {
class K80 : ReversibleModel {
public:
    K80() : m_Kappa(1.0) {
        m_Q = create_Q(m_Kappa);
    };

    K80(double kappa) : ReversibleModel(), m_Kappa(kappa) {
        m_Q = create_Q(m_Kappa);
    };

    inline const Arbor::Matrix &Q() const override { return m_Q; };

    inline Arbor::Matrix P(double time) const override;

    inline void update(double kappa) {
        m_Kappa = kappa;
        m_Q = create_Q(kappa);
    }

private:
    inline Arbor::Matrix create_Q(double kappa);

    double m_Kappa = 1.0;
    Arbor::Matrix m_Q;
};

inline Arbor::Matrix K80::create_Q(double kappa) {
    double diag = -1.0;
    double ts = 1.0 / (kappa + 2.0);
    double tv = kappa / (kappa + 2.0);
    Arbor::Matrix q{
            diag, ts, tv, ts,
            ts, diag, ts, tv,
            tv, ts, diag, ts,
            ts, tv, ts, diag
    };
    return q;
}

inline Arbor::Matrix K80::P(double time) const {
    // eigenvalues
    double e1 = -4.0 * time / (m_Kappa + 2.0);
    double e2 = -2.0 * time * (m_Kappa + 1.0) / (m_Kappa + 2.0);

    // exponentiate the eigenvalues
    e1 = exp(e1);
    e2 = exp(e2);

    // compute repeated elements of P
    double p1 = 0.25 * (1.0 + e1 + 2 * e2);
    double p2 = 0.25 * (1.0 - e1);
    double p3 = 0.25 * (1.0 + e1 - 2 * e2);

    // Construct P
    return Arbor::Matrix{
            p1, p2, p3, p2,
            p2, p1, p2, p3,
            p3, p2, p1, p2,
            p2, p3, p2, p1
    };
}

TEST_CASE("Testing K80") {
    K80 k80(10.0);

    SUBCASE("QMatrix") {
        Arbor::Matrix q = k80.Q();
                REQUIRE(q[0] == doctest::Approx(-1.0));
                REQUIRE(q[1] == doctest::Approx(0.08333333));
                REQUIRE(((q[0] + q[1] + q[2] + q[3]) == doctest::Approx(0.0)));
    }

    SUBCASE("PMatrix") {
        Arbor::Matrix p = k80.P(2.5);
                REQUIRE(p[0] == doctest::Approx(0.36375994));
                REQUIRE((p[0] + p[1] + p[2] + p[3]) == doctest::Approx(1.0));
    }
}
}