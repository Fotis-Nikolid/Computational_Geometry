#!/bin/bash
cgal_create_CMakeLists -s to_polygon
cmake -DCGAL_DIR=/usr/lib/CGAL
make
./to_polygon -i $1 -o $2 -algorithm $3 -edge_selection $4 