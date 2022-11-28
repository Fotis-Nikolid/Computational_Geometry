#ifndef HULL_H
#define HULL_H
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/convex_hull_2.h>
#include <iostream>
#include <vector>
#include <list>
#include <time.h>
#include <cstdlib>


    

template<class Kernel>
class Hull {
    
    typedef CGAL::Polygon_2<Kernel> Polygon_2;
    typedef CGAL::Point_2<Kernel>  Point_2;
    typedef CGAL::Segment_2<Kernel> Segment_2;
    typedef CGAL::Triangle_2<Kernel> Triangle_2;

    private:  
        Edge_2 first_edge;
        Edge_2 last_edge; 
        //data structs       
        std::unordered_map<Segment_2,Point_2> pairings;
        Point_2 last_inserted;
        //internal helper functions
        bool is_visible(Segment_2 Edge,Point_2 n_point,Polygon_2 Polygon);
        double Edge_Selection(Polygon_2& Polygon,std::list<Point_2>& remaining_points,char criteria);
    public:
        double solve(Polygon_2& Polygon,std::list<Point_2> Points,char Criteria);
        

};

#include "hull.cc"

#endif