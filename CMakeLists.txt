#CMake version and project name
cmake_minimum_required(VERSION 3.15)
project(SimpleSolverProject)

# Set the C++ standard (C++17 is a good modern choice)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Create the Solver Library ---
file(GLOB_RECURSE SOLVER_SOURCES "src/*.cpp" "src/*.h")
list(FILTER SOLVER_SOURCES EXCLUDE REGEX ".*main\\.cpp$")

add_library(simple_solver ${SOLVER_SOURCES})

# public headers.
target_include_directories(simple_solver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

# This is the main application that runs the simulation.
add_executable(app src/main.cpp)


# Link the application against the solver library.
target_link_libraries(app PRIVATE simple_solver)

# Install Rules
install(TARGETS app DESTINATION bin)
install(TARGETS simple_solver DESTINATION lib)
install(DIRECTORY include/ DESTINATION include)