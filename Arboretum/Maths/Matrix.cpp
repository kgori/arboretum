//
// Created by Kevin Gori on 07/04/2021.
//

#include <iostream>
#include "Matrix.h"

namespace Arbor {
Matrix4 matmul(const Matrix4 &a, const Matrix4 &b) {
    Matrix4 result{a[0] * b[0] + a[4] * b[1] + a[8] * b[2] + a[12] * b[3],
                   a[1] * b[0] + a[5] * b[1] + a[9] * b[2] + a[13] * b[3],
                   a[2] * b[0] + a[6] * b[1] + a[10] * b[2] + a[14] * b[3],
                   a[3] * b[0] + a[7] * b[1] + a[11] * b[2] + a[15] * b[3],
                   a[0] * b[4] + a[4] * b[5] + a[8] * b[6] + a[12] * b[7],
                   a[1] * b[4] + a[5] * b[5] + a[9] * b[6] + a[13] * b[7],
                   a[2] * b[4] + a[6] * b[5] + a[10] * b[6] + a[14] * b[7],
                   a[3] * b[4] + a[7] * b[5] + a[11] * b[6] + a[15] * b[7],
                   a[0] * b[8] + a[4] * b[9] + a[8] * b[10] + a[12] * b[11],
                   a[1] * b[8] + a[5] * b[9] + a[9] * b[10] + a[13] * b[11],
                   a[2] * b[8] + a[6] * b[9] + a[10] * b[10] + a[14] * b[11],
                   a[3] * b[8] + a[7] * b[9] + a[11] * b[10] + a[15] * b[11],
                   a[0] * b[12] + a[4] * b[13] + a[8] * b[14] + a[12] * b[15],
                   a[1] * b[12] + a[5] * b[13] + a[9] * b[14] + a[13] * b[15],
                   a[2] * b[12] + a[6] * b[13] + a[10] * b[14] + a[14] * b[15],
                   a[3] * b[12] + a[7] * b[13] + a[11] * b[14] + a[15] * b[15]};
    return result;
}

Matrix4 reconstitute_eigen(const Matrix4 &evecs, const Vector4 &v, const Matrix4 &ivecs) {
    Matrix4 result{evecs[0] * v[0] * ivecs[0] + evecs[4] * v[1] * ivecs[1] + evecs[8] * v[2] * ivecs[2] +
                   evecs[12] * v[3] * ivecs[3],
                   evecs[1] * v[0] * ivecs[0] + evecs[5] * v[1] * ivecs[1] + evecs[9] * v[2] * ivecs[2] +
                   evecs[13] * v[3] * ivecs[3],
                   evecs[2] * v[0] * ivecs[0] + evecs[6] * v[1] * ivecs[1] + evecs[10] * v[2] * ivecs[2] +
                   evecs[14] * v[3] * ivecs[3],
                   evecs[3] * v[0] * ivecs[0] + evecs[7] * v[1] * ivecs[1] + evecs[11] * v[2] * ivecs[2] +
                   evecs[15] * v[3] * ivecs[3],
                   evecs[0] * v[0] * ivecs[4] + evecs[4] * v[1] * ivecs[5] + evecs[8] * v[2] * ivecs[6] +
                   evecs[12] * v[3] * ivecs[7],
                   evecs[1] * v[0] * ivecs[4] + evecs[5] * v[1] * ivecs[5] + evecs[9] * v[2] * ivecs[6] +
                   evecs[13] * v[3] * ivecs[7],
                   evecs[2] * v[0] * ivecs[4] + evecs[6] * v[1] * ivecs[5] + evecs[10] * v[2] * ivecs[6] +
                   evecs[14] * v[3] * ivecs[7],
                   evecs[3] * v[0] * ivecs[4] + evecs[7] * v[1] * ivecs[5] + evecs[11] * v[2] * ivecs[6] +
                   evecs[15] * v[3] * ivecs[7],
                   evecs[0] * v[0] * ivecs[8] + evecs[4] * v[1] * ivecs[9] + evecs[8] * v[2] * ivecs[10] +
                   evecs[12] * v[3] * ivecs[11],
                   evecs[1] * v[0] * ivecs[8] + evecs[5] * v[1] * ivecs[9] + evecs[9] * v[2] * ivecs[10] +
                   evecs[13] * v[3] * ivecs[11],
                   evecs[2] * v[0] * ivecs[8] + evecs[6] * v[1] * ivecs[9] + evecs[10] * v[2] * ivecs[10] +
                   evecs[14] * v[3] * ivecs[11],
                   evecs[3] * v[0] * ivecs[8] + evecs[7] * v[1] * ivecs[9] + evecs[11] * v[2] * ivecs[10] +
                   evecs[15] * v[3] * ivecs[11],
                   evecs[0] * v[0] * ivecs[12] + evecs[4] * v[1] * ivecs[13] + evecs[8] * v[2] * ivecs[14] +
                   evecs[12] * v[3] * ivecs[15],
                   evecs[1] * v[0] * ivecs[12] + evecs[5] * v[1] * ivecs[13] + evecs[9] * v[2] * ivecs[14] +
                   evecs[13] * v[3] * ivecs[15],
                   evecs[2] * v[0] * ivecs[12] + evecs[6] * v[1] * ivecs[13] + evecs[10] * v[2] * ivecs[14] +
                   evecs[14] * v[3] * ivecs[15],
                   evecs[3] * v[0] * ivecs[12] + evecs[7] * v[1] * ivecs[13] + evecs[11] * v[2] * ivecs[14] +
                   evecs[15] * v[3] * ivecs[15]};
    return result;
}

void print_matrix(const Matrix4 &m) {
    std::cout << "[\n" << " " << m[0] << ", " << m[4] << ", " << m[8] << ", " << m[12] << '\n'
              << " " << m[1] << ", " << m[5] << ", " << m[9] << ", " << m[13] << '\n'
              << " " << m[2] << ", " << m[6] << ", " << m[10] << ", " << m[14] << '\n'
              << " " << m[3] << ", " << m[7] << ", " << m[11] << ", " << m[15] << '\n'
              << "]" << std::endl;
}
}

