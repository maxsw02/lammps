#!/bin/bash

# Add libtorch
cat >> cmake/CMakeLists.txt << "EOF2"

find_package(Torch REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")
target_link_libraries(lammps PUBLIC "${TORCH_LIBRARIES}")
EOF2

echo "Done!"
