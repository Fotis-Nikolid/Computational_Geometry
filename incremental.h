#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <iostream>
#include <vector>

//flags
#define F1A 0
#define F1B 1
#define F2A 2
#define F2B 3

template<class Kernel>
class Incremental {
    
    typedef CGAL::Polygon_2<Kernel> Polygon_2;
    typedef Kernel::Point_2  Point_2;
    typedef Kernel::Segment_2 Segment;

    private: 
        Polygon_2 Convex_Hull_Polygon;
        Polygon_2 Real_Polygon;
        float Area;

        void Initialize(std::vector<Point_2>);
        bool comp_x_less(Point_2 p1, Point_2 p2);
        bool comp_x_more(Point_2 p1, Point_2 p2);
        bool comp_y_less(Point_2 p1, Point_2 p2);
        bool comp_y_more(Point_2 p1, Point_2 p2);
        void Sort(std::vector<Point_2>,int);

    public:
        Incremental(vector<Polygon_2>, int);
        float getPolygonArea();
        Polygon_2 getPolygon();
};