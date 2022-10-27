#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <iostream>
#include <vector>

template<class Kernel>
class Incremental 
{
    
    typedef CGAL::Polygon_2<Kernel> Polygon_2;
    typedef Kernel::Point_2  Point_2;
    typedef Kernel::Segment_2 Segment;

    private: 
        Polygon_2 Convex_Hull_Polygon;
        Polygon_2 Real_Polygon;

        void Initialize(std::vector<Point_2>, char);
        /*bool comp_x_less(Point_2 p1, Point_2 p2);
        bool comp_x_more(Point_2 p1, Point_2 p2);
        bool comp_y_less(Point_2 p1, Point_2 p2);
        bool comp_y_more(Point_2 p1, Point_2 p2);*/
        void Sort(std::vector<Point_2>, std::string);
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
        Incremental(const vector<Polygon_2>, std::string, std::string);
        float getPolygonArea();
        Polygon_2 getPolygon();
};