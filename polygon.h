#ifndef POLYGON_H
#define POLYGON_H
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <iostream>
#include <vector>
#include "simulated_annealing.h"
#include "local_search.h"

template <class Kernel> class Polygon {
    typedef CGAL::Polygon_2<Kernel> Polygon_2;
    typedef CGAL::Point_2<Kernel>  Point_2;
    typedef CGAL::Segment_2<Kernel> Segment_2;

    private:
        Polygon_2 pol;
        float dt_Area;
        float initial_area;
        
    public:
        //Given a set of Points and a choice of (Agorithm,Criteria) creates a simple polygon
        Polygon(std::vector<Point_2>, std::string algorithm, std::string step_choice,std::string criteria,int L,int threashold,int Attempts=1);        
        int Size();
        Polygon_2 get_Polygon();
        float Area();
        float Init_Area();
        bool Simple();
};

#include "polygon.cc"

#endif
