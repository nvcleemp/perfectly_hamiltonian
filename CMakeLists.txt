cmake_minimum_required(VERSION 3.15)
project(perfectly_hamiltonian C)

set(CMAKE_C_STANDARD 99)

include_directories ("${PROJECT_SOURCE_DIR}/planegraph")
add_subdirectory (planegraph)

include_directories ("${PROJECT_SOURCE_DIR}/cprogutil")
add_subdirectory (cprogutil)

add_executable(reg_planar_perfectly_hamiltonian reg_planar_perfectly_hamiltonian.c)
target_link_libraries (reg_planar_perfectly_hamiltonian planegraph cprogutil)

add_executable(5reg_planar_no_cyclic_6_edge_cut 5reg_planar_no_cyclic_6_edge_cut.c)
target_link_libraries (5reg_planar_no_cyclic_6_edge_cut planegraph)

add_executable(5reg_planar_has_cyclic_7_edge_cut 5reg_planar_has_7_edge_cut.c)
target_link_libraries (5reg_planar_has_cyclic_7_edge_cut planegraph)
