//
// Created by Kevin Gori on 19/10/2020.
//

#pragma once
#include <string>
#include "../vendor/doctest.h"

namespace Arbor {
class Hello {
public:
    Hello();

    ~Hello();

    static std::string greet(const std::string& message);
    static std::string greet(const char* message);
};
}


TEST_CASE("Hello") {
    CHECK(Arbor::Hello::greet("Test Case") == "Hello, Test Case, from Arboretum!");
}
