#ifndef POLYGON_H
#define POLYGON_H
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <iostream>
#include <vector>
#include "hull.h"
#include "incremental.h"

template <class Kernel>
class Polygon {
    typedef CGAL::Polygon_2<Kernel> Polygon_2;
    typedef CGAL::Point_2<Kernel>  Point_2;
    typedef CGAL::Segment_2<Kernel> Segment_2;

    private:
        Polygon_2 pol;
        float dt_Area;
        
    public:
        //Given a set of Points and a choice of (Agorithm,Criteria) creates a simple polygon
        void Algorithm(std::vector<Point_2>, std::string);
        void Algorithm(std::vector<Point_2>, std::string, std::string);
        int Size();
        Polygon_2 get_Polygon();
        float Area();
        float Ratio();
        bool Simple();
};


template<class Kernel>
void Polygon<Kernel>::Algorithm(std::vector<Point_2> Points, std::string Edge_Selection)
{
    Hull<Kernel> alg;
    std::list<Point_2> list(Points.begin(),Points.end());
    dt_Area = alg.solve(pol, list, *(Edge_Selection.c_str()));
}

template<class Kernel>
void Polygon<Kernel>::Algorithm(std::vector<Point_2> Points, std::string Edge_Selection, std::string Sorting)
{
    Incremental<Kernel> inc(Points, Sorting, *(Edge_Selection.c_str()));
    dt_Area = inc.getPolygonArea();
    pol = inc.getPolygon();
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
float Polygon<Kernel>::Ratio() {
    
}

template<class Kernel>
bool Polygon<Kernel>::Simple() 
{
    return pol.is_simple();    
}

#endif
