add_executable(TestApp
        src/main.cpp)
target_compile_definitions(TestApp PRIVATE DOCTEST_CONFIG_DISABLE)
target_link_libraries(TestApp PUBLIC Arboretum)
target_include_directories(TestApp PUBLIC
        "${PROJECT_SOURCE_DIR}/Arboretum/include")