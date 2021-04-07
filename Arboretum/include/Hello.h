//
// Created by Kevin Gori on 19/10/2020.
//

#pragma once
#include <string>
#include "doctest.h"

namespace Arbor {

class Hello {
public:
    static std::string greet(const std::string &message);
    static std::string greet(const char *message);
};


TEST_CASE ("Hello") {
    CHECK(Hello::greet("Test Case") == "Hello, Test Case, from Arboretum!");
}

}