//
// Created by Kevin Gori on 19/10/2020.
//

#include "Hello.h"

namespace Arbor {

std::string Hello::greet(const std::string& message) {
    return std::string("Hello, " + message + ", from Arboretum!");
}

std::string Hello::greet(const char* message) {
    return std::string("Hello, ") + std::string(message) + std::string(", from Arboretum!");
}
}
