//
// Created by Kevin Gori on 07/04/2021.
//

#include <array>
#include "doctest.h"
#pragma once

namespace Arbor {
// Aliases for 4x4 matrix and length-4 vector
using Matrix4 = std::array<double, 16>;
using Vector4 = std::array<double, 4>;

// Matrix multiply two 4x4 Matrices
Matrix4 matmul(const Matrix4 &a, const Matrix4 &b);

/* Perform the multiplication LVL-1, where L is a 4x4 matrix of eigenvectors,
 * V is a diagonal matrix of eigenvalues (represented by a 4-vector), and L-1
 * is the inverse of L.
 */
Matrix4 reconstitute_eigen(const Matrix4 &evecs, const Vector4 &diag, const Matrix4 &ivecs);

// Utility to print a matrix to std::cout
void print_matrix(const Matrix4 &m);
}

TEST_CASE("Testing reconstitute eigen") {
    Arbor::Matrix4 m1{1.0, 2.0, 3.0, 4.0,
                      5.0, 6.0, 7.0, 8.0,
                      9.0, 10.0, 11.0, 12.0,
                      13.0, 14.0, 15.0, 16.0};
    Arbor::Vector4 v1{8.0, 8.5, 9.0, 9.5};
    Arbor::Matrix4 m2{0.1, 0.2, 0.3, 0.4,
                      0.5, 0.6, 0.7, 0.8,
                      0.9, 0.10, 0.11, 0.12,
                      0.13, 0.14, 0.15, 0.16};

    Arbor::Matrix4 r = Arbor::reconstitute_eigen(m1, v1, m2);

    REQUIRE(r[0] == doctest::Approx(83));
    REQUIRE(r[1] == doctest::Approx(92));
    REQUIRE(r[2] == doctest::Approx(101));
    REQUIRE(r[3] == doctest::Approx(110));
    REQUIRE(r[4] == doctest::Approx(185));
    REQUIRE(r[5] == doctest::Approx(208));
    REQUIRE(r[6] == doctest::Approx(231));
    REQUIRE(r[7] == doctest::Approx(254));
    REQUIRE(r[8] == doctest::Approx(35.18));
    REQUIRE(r[9] == doctest::Approx(45.36));
}
