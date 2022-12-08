#ifndef POLYGON_H
#define POLYGON_H
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <iostream>
#include <vector>
#include "simulated_annealing.h"
#include "localsearch.h"

template <class Kernel> class Polygon {
    typedef CGAL::Polygon_2<Kernel> Polygon_2;
    typedef CGAL::Point_2<Kernel>  Point_2;
    typedef CGAL::Segment_2<Kernel> Segment_2;

    private:
        bool failed=false;
        Polygon_2 pol;
        double dt_Area;
        double initial_area;
        
    public:
        //Given a set of Points and a choice of (Agorithm,Criteria) creates a simple polygon
        Polygon(std::vector<Point_2>, std::string Algorithm,std::string Criteria,std::string Step_Choice,int L = 0,double threshold = 0.0, int K = 0,int Attempts=1);        
        int Size();
        Polygon_2 get_Polygon();
        double Area();
        double Init_Area();
        bool Simple();
        bool Success();
};

#include "polygon.cc"

#endif
