#ifndef LOCALSEARCH_H
#define LOCALSEARCH_H
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <vector>

template<class Kernel> class LocalSearch
{
    typedef CGAL::Polygon_2<Kernel> Polygon_2;
    typedef CGAL::Point_2<Kernel>  Point_2;
    typedef CGAL::Segment_2<Kernel> Segment_2;
    typedef CGAL::Triangle_2<Kernel> Triangle_2;

    private:
        Polygon_2 Polygon;
        bool (*compare)(Polygon_2, Polygon_2);

        bool visible_points(const Point_2, const Point_2);
        void relocate_edges(Polygon_2&, int, int, int, int);
        double swap_L_with_edge(const int, const int);
        bool solve_specific_K(const int, const int, const double);
    
    public:
        LocalSearch(const std::vector<Point_2>, const std::string, const int, const int, const double);
};

#include "localsearch.cpp"

#endif