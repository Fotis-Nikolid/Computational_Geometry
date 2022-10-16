#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/convex_hull_2.h>
#include <iostream>
#include <vector>
#include <list>

    

template<class Kernel>
class Hull {
    
    typedef CGAL::Polygon_2<Kernel> Polygon_2;
    typedef CGAL::Point_2<Kernel>  Point_2;
    typedef CGAL::Segment_2<Kernel> Segment_2;
    typedef CGAL::Triangle_2<Kernel> Triangle_2;

    private: 
        float Area;
        
        //internal helper functions
        bool overlaps_point(Triangle_2 triangle,std::vector<Point_2> Points);
        bool is_visible(Segment_2 Edge,Point_2 n_point,Polygon_2 Polygon);
        double Edge_Selection(Polygon_2 Polygon,Point_2 n_point,std::list<Point_2> remaining_points,std::string criteria);
    public:
        double solve(Polygon_2& Polygon,std::list<Point_2>& Points,std::string Criteria);
        

};