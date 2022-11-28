#ifndef POLYGON_CC
#define POLYGON_CC
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <iostream>
#include <vector>
#include "polygon.h"

//receives a vector of points, and an algorithm with associated criteria, and then creates the polygon based on the algorithm chosen
template<class Kernel> Polygon<Kernel>::Polygon(std::vector<Point_2> Points, std::string Algorithm, std::string Criteria,std::string Step_Choice="local",Attempts=1,Iterations=100)
{
    if(Criteria!="max" && Criteria!="min")
    {
        std::cerr << "An criteria argument must be given (max or min)" << std::endl;
        exit(1);
    }
    if(Algorithm == "local_step")
    {

    }
    else if(Algorithm == "simulated_annealing")
    {
        Simulated_Annealing algorithm;
        if (Step_Choice=="local" || Step_Choice=="global") {
            dt_Area=algorithm.solve(pol,Points,Criteria,Step_Choice,Iterations,Attempts);
        }
        else if (Step_Choice=="subdivision") {
            dt_Area=algorithm.expanded_solve(pol,Points,Criteria,Iterations,Attempts);
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