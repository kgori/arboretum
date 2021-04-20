find_package(Eigen3 CONFIG REQUIRED)

add_library(Arboretum SHARED
        src/Hello.cpp
        SubstitutionModels/SubstitutionModels.h
        SubstitutionModels/K80.h
        SubstitutionModels/TN93.h
        SubstitutionModels/HKY85.h
        SubstitutionModels/DNAModels.h
        Maths/Matrix.h
        Maths/Matrix.cpp
        SubstitutionModels/F81.h
        SubstitutionModels/F84.h
        SubstitutionModels/GTR.h
        SubstitutionModels/GTR.cpp)

target_compile_definitions(Arboretum PRIVATE DOCTEST_CONFIG_DISABLE)
target_compile_options(Arboretum PUBLIC "-mavx2" "-mfma")

target_include_directories(Arboretum PRIVATE
        ${EIGEN3_INCLUDE_DIR}
        ${PROJECT_SOURCE_DIR}/Arboretum
        ${PROJECT_SOURCE_DIR}/Arboretum/include
        ${PROJECT_SOURCE_DIR}/Arboretum/vendor)

add_executable(Arboretum_test
        test/main.cpp)

target_include_directories(Arboretum_test PRIVATE
        ${EIGEN3_INCLUDE_DIR}
        ${PROJECT_SOURCE_DIR}/Arboretum
        ${PROJECT_SOURCE_DIR}/Arboretum/include
        ${PROJECT_SOURCE_DIR}/Arboretum/vendor)

target_link_libraries(Arboretum_test PUBLIC Arboretum)
target_link_libraries(Arboretum_test PRIVATE Eigen3::Eigen)
