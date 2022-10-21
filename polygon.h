#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <iostream>
#include <vector>


template <class Kernel>
class Polygon {
    typedef CGAL::Polygon_2<Kernel> Polygon_2;
    typedef CGAL::Point_2<Kernel>  Point_2;
    typedef CGAL::Segment_2<Kernel> Segment_2;

    private:
        Polygon_2 Polygon;
        float dt_Area;
        
    public:
        //Given a set of Points and a choice of (Agorithm,Criteria) creates a simple polygon 
        void Hull_Based(std::vector<Point_2> Points,std::string Edge_Selection);
        void Increment_Based(std::vector<Point_2> Points,std::string Edge_Selection,std::string Sorting);
        int Size();
        std::vector<Point_2> Points();
        std::vector<Segment_2> Edges();
        float Area();
        float Ratio();
};

