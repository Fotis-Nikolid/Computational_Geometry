#if 0
#ifndef POLYGON_CC
#define POLYGON_CC
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <iostream>
#include <vector>
#include "polygon.h"

//receives a vector of points, and an algorithm with associated criteria, and then creates the polygon based on the algorithm chosen
template<class Kernel> Polygon<Kernel>::Polygon(std::vector<Point_2> Points, std::string Algorithm,std::string Criteria,std::string Step_Choice,std::string Initialization,int L,double threshold,int K,int Attempts)
{
    if(Criteria!="max" && Criteria!="min")
    {
        std::cerr << "An criteria argument must be given (max or min)" << std::endl;
        exit(1);
    }
    if (Initialization!="incremental" && Initialization!="convex_hull") {
        std::cerr << "An Initialization algorithm must be ginen (incremental or convex_hull)" << std::endl;
        exit(1);
    }
    if(Algorithm == "local_search")
    {
        LocalSearch<Kernel> loc(Points, Criteria,Initialization);

        initial_area = loc.getPolygonArea();

        failed = !loc.MinimizePolygon(L, threshold, K);

        pol = loc.getPolygon();

        dt_Area = loc.getPolygonArea();
    }
    else if(Algorithm == "simulated_annealing")
    {
        Simulated_Annealing<Kernel> algorithm;
        if (Step_Choice=="local" || Step_Choice=="global") {
            if (algorithm.solve(pol,Points,Criteria,Step_Choice,Initialization,L,Attempts,initial_area)) {
                dt_Area=abs(pol.area());
            }
            else {
                failed=true;
            }

        }
        else if (Step_Choice=="subdivision") {
            if (algorithm.sub_division(pol,Points,Criteria,Initialization,L,initial_area)) {
                dt_Area=abs(pol.area());
            }
            else {
                failed=true;
            }
        }
    }
    else
    {
        std::cerr << "Algorithm " + Algorithm + " not implemented" << std::endl;
        exit(1);
    }
}

template<class Kernel>
int Polygon<Kernel>::Size() {
    return pol.vertices().size();
}

template<class Kernel>
bool Polygon<Kernel>::Success() {
    return !failed;
}


template<class Kernel>
CGAL::Polygon_2<Kernel> Polygon<Kernel>::get_Polygon() {
    return this->pol;
}

template<class Kernel>
double Polygon<Kernel>::Area() {
    return dt_Area;
}
template<class Kernel>
double Polygon<Kernel>::Init_Area() {
    return initial_area;
}
template<class Kernel>
bool Polygon<Kernel>::Simple() 
{
    return pol.is_simple();    
}

#endif
#endif