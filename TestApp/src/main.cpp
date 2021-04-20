#include <iostream>
#include <numeric>
#include "Arboretum.h"

int main() {
    Arbor::Hello hello;
    std::string greeting = Arbor::Hello::greet("TestApp");
    std::cout << greeting << std::endl;

    Arbor::K80 k80(10);
    auto q_matrix = k80.Q();
    std::cout << "K80(10.0) Q...\n";
    Arbor::print_matrix(q_matrix);
    k80.update(4.5);
    q_matrix = k80.Q();
    std::cout << "K80(4.5) Q...\n";
    Arbor::print_matrix(q_matrix);

    auto p_matrix = k80.P(2.0);
    std::cout << "K80(4.5)::P(2.0)...\n";
    Arbor::print_matrix(p_matrix);

    Arbor::TN93 tn93(4.5, 4.5, 1.0, 0.25, 0.25, 0.25, 0.25);
    q_matrix = tn93.Q();
    std::cout << "TN93(4.5) Q...\n";
    Arbor::print_matrix(q_matrix);

    p_matrix = tn93.P(2.0);
    std::cout << "TN93(4.5)::P(2.0)...\n";
    Arbor::print_matrix(p_matrix);

    Arbor::GTR gtr(2.0, 3.0, 4.0, 5.0, 6.0, 0.1, 0.2, 0.3, 0.4);
    q_matrix = gtr.Q();
    std::cout << "GTR(2,3,4,5,6,0.1,0.2,0.3,0.4) Q...\n";
    Arbor::print_matrix(q_matrix);

#ifdef __SSE__
    std::cout << "SSE available\n";
#endif

#ifdef __SSE2__
    std::cout << "SSE2 available\n";
#endif

#ifdef __SSE3__
    std::cout << "SSE3 available\n";
#endif

#ifdef __AVX__
    std::cout << "AVX available\n";
#endif

#ifdef __AVX2__
    std::cout << "AVX2 available\n";
#endif

#ifdef __FMA__
    std::cout << "FMA available\n";
#endif


    using Arbor::Vector4;
    using Arbor::Matrix4;

    Arbor::K80 subsmod{2.0};
    auto p_matrix_inner = subsmod.P(0.1);
    auto p_matrix_outer = subsmod.P(0.2);

    std::cout << "Probability matrix for inner branches" << std::endl;
    Arbor::print_matrix(p_matrix_inner);
    std::cout << "Probability matrix for outer branches" << std::endl;
    Arbor::print_matrix(p_matrix_outer);

    Vector4 const /* T */ l1{0, 0, 0, 1};
    Vector4 const /* C */ l2{0, 1, 0, 0};
    Vector4 const /* A */ l3{1, 0, 0, 0};
    Vector4 const /* C */ l4{0, 1, 0, 0};
    Vector4 const /* C */ l5{0, 1, 0, 0};

    Vector4 n7 = (p_matrix_outer * l1) * (p_matrix_outer * l2);
    Vector4 n6 = (p_matrix_inner * n7) * (p_matrix_outer * l3);
    Vector4 n8 = (p_matrix_outer * l4) * (p_matrix_outer * l5);
    Vector4 n0 = (p_matrix_inner * n6) * (p_matrix_inner * n8);

    Arbor::print_vector(n0);
    Arbor::print_vector(n6);
    Arbor::print_vector(n8);
    Arbor::print_vector(n7);

    return 0;
}
