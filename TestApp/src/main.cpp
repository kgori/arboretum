#include <iostream>
#include <numeric>
#include "Arboretum.h"
#include "SubstitutionModels/K80.h"

int main() {
    Arbor::Hello hello;
    std::string greeting = Arbor::Hello::greet("TestApp");
    std::cout << greeting << std::endl;

    Arbor::K80 k80(10);
    auto q_matrix = k80.Q();
    std::cout << "Q... " << q_matrix[2] << "\n";
    k80.update(4.5);
    q_matrix = k80.Q();
    std::cout << "Q... " << q_matrix[2] << "\n";

    auto p_matrix = k80.P(2.0);
    std::cout << "P... ["
        << p_matrix[0] << ", "
        << p_matrix[1] << ", "
        << p_matrix[2] << ", "
        << p_matrix[3] << "]"
        << "\n";
    return 0;
}
