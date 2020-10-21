add_executable(Arboretum_test
        test/main.cpp)

target_include_directories(Arboretum_test PRIVATE
        ${PROJECT_SOURCE_DIR}/Arboretum/include
        ${PROJECT_SOURCE_DIR}/Arboretum/vendor)

target_link_libraries(Arboretum_test PUBLIC Arboretum)

add_library(Arboretum SHARED
        src/Hello.cpp
        include/SubstitutionModels/SubstitutionModels.h
        include/SubstitutionModels/K80.h)

target_compile_definitions(Arboretum PRIVATE DOCTEST_CONFIG_DISABLE)

target_include_directories(Arboretum PRIVATE
        ${PROJECT_SOURCE_DIR}/Arboretum/include
        ${PROJECT_SOURCE_DIR}/Arboretum/vendor)

