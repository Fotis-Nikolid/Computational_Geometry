#ifdef INCREMENTAL_H
#define INCREMENTAL_H
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <limits>
#include <cstdlib>
#include <ctime>

template<class Kernel>
class Incremental 
{
    
    typedef CGAL::Polygon_2<Kernel> Polygon_2;
    typedef CGAL::Point_2<Kernel>  Point_2;
    typedef CGAL::Segment_2<Kernel> Segment_2;
    typedef CGAL::Triangle_2<Kernel> Triangle_2;

    private: 
        Polygon_2 Convex_Hull_Polygon;
        Polygon_2 Real_Polygon;

        void Initialize(std::vector<Point_2>&, char);
        static bool comp_x_less(Point_2 p1, Point_2 p2);
        static bool comp_x_more(Point_2 p1, Point_2 p2);
        static bool comp_y_less(Point_2 p1, Point_2 p2);
        static bool comp_y_more(Point_2 p1, Point_2 p2);
        void Sort(std::vector<Point_2>&, std::string);
        bool red_visible(Segment_2, Point_2);
        bool visible(Segment_2, Point_2);

        class RedEdgesBoundaries
        {
            public:
                Point_2 first_vertex;
                Point_2 second_vertex;
        };

        RedEdgesBoundaries find_red_edges_boundaries_and_recreate_convex_hull(Point_2);

        void construct_new_polygon(RedEdgesBoundaries, Point_2, char);

    public:
        Incremental(const std::vector<Point_2>, std::string, char);
        float getPolygonArea();
        Polygon_2 getPolygon();
};

#endif