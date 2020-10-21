//
// Created by Kevin Gori on 19/10/2020.
//

#pragma once
#include <vector>
#include <cmath>

#include "../../vendor/doctest.h"

template<typename T>
std::initializer_list<T> make_init_list (std::initializer_list<T> && l) {
    return l;
}

namespace Arbor {

using Matrix = std::vector<double>;

class ReversibleModel {
public:

protected:
    ReversibleModel() = default;

    virtual inline Matrix P(double time) const = 0;
    virtual inline const Matrix& Q() const = 0;
};

}
