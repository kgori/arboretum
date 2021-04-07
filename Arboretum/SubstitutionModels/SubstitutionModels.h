//
// Created by Kevin Gori on 19/10/2020.
//

#pragma once
#include <array>
#include <cmath>
#include "Maths/Matrix.h"

namespace Arbor {

class ReversibleModel {
public:

protected:
    ReversibleModel() = default;

    virtual inline const Matrix4 P(double time) const = 0;
    virtual inline const Matrix4& Q() const = 0;
};

}
