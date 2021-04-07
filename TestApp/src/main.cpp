#include <iostream>
#include <numeric>
#include "Arboretum.h"
#include "SubstitutionModels/SubstitutionModels.h"
#include "SubstitutionModels/K80.h"
#include "SubstitutionModels/TN93.h"

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


    return 0;
}
