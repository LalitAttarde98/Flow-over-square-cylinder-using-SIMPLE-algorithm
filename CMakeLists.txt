# Set the minimum CMake version and project name
cmake_minimum_required(VERSION 3.15)
project(SimpleSolverProject)

# Set the C++ standard (C++17 is a good modern choice)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# --- 1. Create the Solver Library ---
# This library will contain all the simulation logic.
# We add all .cpp files from the src/ directory, *except* main.cpp
file(GLOB_RECURSE SOLVER_SOURCES "src/*.cpp" "src/*.h")
list(FILTER SOLVER_SOURCES EXCLUDE REGEX ".*main\\.cpp$")

add_library(simple_solver ${SOLVER_SOURCES})

# Tell the library where to find its own public headers.
# This also makes it so any app that links this library
# automatically knows about the 'include' directory.
target_include_directories(simple_solver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

# --- 2. Create the Executable ---
# This is the main application that runs the simulation.
add_executable(app src/main.cpp)


# Link the application against the solver library.
# The 'app' needs the logic from 'simple_solver'.
target_link_libraries(app PRIVATE simple_solver)

# --- 3. (Optional) Install Rules ---
# This tells CMake how to install your app and library
install(TARGETS app DESTINATION bin)
install(TARGETS simple_solver DESTINATION lib)
install(DIRECTORY include/ DESTINATION include)