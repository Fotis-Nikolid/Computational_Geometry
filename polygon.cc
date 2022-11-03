#ifndef POLYGON_CC
#define POLYGON_CC
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <iostream>
#include <vector>
#include "polygon.h"

//receives a vector of points, and an algorithm with associated criteria, and then creates the polygon based on the algorithm chosen
template<class Kernel> Polygon<Kernel>::Polygon(std::vector<Point_2> Points, std::string alg, std::string Edge_Selection, std::string Sorting)
{
    if(Edge_Selection != "1" && Edge_Selection != "2" && Edge_Selection != "3")
    {
        std::cerr << "An -edge_selection argument must be given ( 1:random area, 2:minimum area, 3:max area )" << std::endl;
        exit(1);
    }
    if(alg == "convex_hull")
    {
        Hull<Kernel> alg;
        std::list<Point_2> list(Points.begin(),Points.end());
        dt_Area = alg.solve(pol, list, *(Edge_Selection.c_str()));
    }
    else if(alg == "incremental")
    {
        Incremental<Kernel> inc(Points, Sorting, *(Edge_Selection.c_str()));
        dt_Area = inc.getPolygonArea();
        pol = inc.getPolygon();
    }
    else
    {
        std::cerr << "Algorithm " + alg + " not implemented" << std::endl;
        exit(1);
    }
}

template<class Kernel>
int Polygon<Kernel>::Size() {
    return pol.vertices().size();
}

template<class Kernel>
CGAL::Polygon_2<Kernel> Polygon<Kernel>::get_Polygon() {
    return this->pol;
}

template<class Kernel>
float Polygon<Kernel>::Area() {
    return dt_Area;
}
template<class Kernel>
bool Polygon<Kernel>::Simple() 
{
    return pol.is_simple();    
}

#endif