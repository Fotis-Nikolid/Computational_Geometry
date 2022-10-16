#!/bin/bash
cgal_create_CMakeLists -s project
cmake -DCGAL_DIR=/usr/lib/CGAL
make