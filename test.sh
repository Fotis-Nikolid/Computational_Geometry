#!/bin/bash
cgal_create_CMakeLists -s test
cmake -DCGAL_DIR=/usr/lib/CGAL
make
./test $1 $2 $3