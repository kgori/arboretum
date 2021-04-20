//
// Created by Kevin Gori on 07/04/2021.
//

#include <iostream>
#include "Matrix.h"

namespace Arbor {


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

#ifdef __AVX__
    #include <immintrin.h>
    #ifdef __FMA__
        /* FMA implementation of Matrix4 * Vector4 */
        void matmulvec_impl(const Matrix4& m, const Vector4& vin, Vector4& vout) {
            __m256d sum = _mm256_set1_pd(0.0);
            __m256d column;
            __m256d scalar;

            for (int i = 0; i < 4; ++i) {
                column = _mm256_load_pd(&m[i*4]);
                scalar = _mm256_set1_pd(vin[i]);
                sum = _mm256_fmadd_pd(column, scalar, sum);
            }

            _mm256_store_pd(&vout[0], sum);
        }

        void matmulvec_impl(const Matrix4& m, const double* vin, double* vout) {
            __m256d sum = _mm256_set1_pd(0.0);
            __m256d column;
            __m256d scalar;

            for (int i = 0; i < 4; ++i) {
                column = _mm256_load_pd(&m[i*4]);
                scalar = _mm256_set1_pd(vin[i]);
                sum = _mm256_fmadd_pd(column, scalar, sum);
            }

            _mm256_store_pd(&vout[0], sum);
        }

        Vector4 matmulvec(const Matrix4& m, const Vector4& v) {
            Vector4 result{0, 0, 0, 0};
            matmulvec_impl(m, v, result);
            return result;
        }
    #else
        /* AVX (no FMA) implementation of Matrix4 * Vector4 */
        void matmulvec_impl(const Matrix4& m, const Vector4& vin, Vector4& vout) {
            __m256d sum = _mm256_set1_pd(0.0);
            __m256d column;
            __m256d scalar;

            for (int i = 0; i < 4; ++i) {
                column = _mm256_load_pd(&m[i*4]);
                scalar = _mm256_set1_pd(vin[i]);
                sum = _mm256_add_pd(_mm256_mul_pd(column, scalar), sum);
            }

            _mm256_store_pd(&vout[0], sum);
        }

        void matmulvec_impl(const Matrix4& m, const double* vin, double* vout) {
            __m256d sum = _mm256_set1_pd(0.0);
            __m256d column;
            __m256d scalar;

            for (int i = 0; i < 4; ++i) {
                column = _mm256_load_pd(&m[i*4]);
                scalar = _mm256_set1_pd(vin[i]);
                sum = _mm256_add_pd(_mm256_mul_pd(column, scalar), sum);
            }

            _mm256_store_pd(&vout[0], sum);
        }

        Vector4 matmulvec(const Matrix4& m, const Vector4& v) {
            Vector4 result{0, 0, 0, 0};
            matmulvec_impl(m, v, result);
            return result;
        }


#endif
    /* AVX/FMA common implemetation of elementwise Vector4 * Vector4 multiplication */
    Vector4 operator*(const Vector4& v1, const Vector4& v2) {
        Vector4 result{0, 0, 0, 0};
        __m256d vv1 = _mm256_load_pd(&v1[0]);
        __m256d vv2 = _mm256_load_pd(&v2[0]);
        __m256d vresult = _mm256_mul_pd(vv1, vv2);
        _mm256_store_pd(&result[0], vresult);
        return result;
    }

    Matrix4 matmul(const Matrix4 &a, const Matrix4 &b) {
        Matrix4 result{0};
        for (int i = 0; i < 4; ++i) {
            matmulvec_impl(a, &b[i*4], &result[i*4]);
        }
        return result;
    }
#else
    Vector4 matmulvec(const Matrix4& m, const Vector4& v) {
        Vector4 result{
                m[0] * v[0] + m[4] * v[1] + m[8]  * v[2] + m[12] * v[3];
                m[1] * v[0] + m[5] * v[1] + m[9]  * v[2] + m[13] * v[3];
                m[2] * v[0] + m[6] * v[1] + m[10] * v[2] + m[14] * v[3];
                m[3] * v[0] + m[7] * v[1] + m[11] * v[2] + m[15] * v[3];
        };
        return result;
    }

    Vector4 operator*(const Vector4& v1, const Vector4& v2) {
        Vector4 result{
                v1[0] * v2[0],
                v1[1] * v2[1],
                v1[2] * v2[2],
                v1[3] * v2[3]
        };
        return result;
    }

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
#endif

void print_matrix(const Matrix4 &m) {
    std::cout << "[\n" << " " << m[0] << ", " << m[4] << ", " << m[8] << ", " << m[12] << '\n'
              << " " << m[1] << ", " << m[5] << ", " << m[9] << ", " << m[13] << '\n'
              << " " << m[2] << ", " << m[6] << ", " << m[10] << ", " << m[14] << '\n'
              << " " << m[3] << ", " << m[7] << ", " << m[11] << ", " << m[15] << '\n'
              << "]" << std::endl;
}

void print_vector(const Vector4 &v) {
    std::cout << "[" << v[0] << ", " << v[1] << ", " << v[2] << ", " << v[3] << "]" << std::endl;
}
}
