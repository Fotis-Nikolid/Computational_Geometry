#ifndef INCREMENTAL_H
#define INCREMENTAL_H
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <vector>

template<class Kernel> class Incremental 
{
    
    typedef CGAL::Polygon_2<Kernel> Polygon_2;
    typedef CGAL::Point_2<Kernel>  Point_2;
    typedef CGAL::Segment_2<Kernel> Segment_2;
    typedef CGAL::Triangle_2<Kernel> Triangle_2;

    private: 
        Polygon_2 Convex_Hull_Polygon;
        Polygon_2 Real_Polygon;

        void Initialize(std::vector<Point_2>&, const char);
        void Sort(std::vector<Point_2>&, const std::string);
        bool red_visible(const Segment_2, const Point_2);
        bool visible(const Segment_2, const Point_2);

        class RedEdgesBoundaries
        {
            public:
                Point_2 first_vertex;
                Point_2 second_vertex;
        };

        RedEdgesBoundaries find_red_edges_boundaries_and_recreate_convex_hull(const Point_2);

        void construct_new_polygon(const RedEdgesBoundaries, const Point_2, const char);

    public:
        Incremental(const std::vector<Point_2>, const std::string, const char);
        float getPolygonArea();
        Polygon_2 getPolygon();
};

#include "incremental.cpp"

#endif