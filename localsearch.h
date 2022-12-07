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

    private:
        Polygon_2 Polygon;
        bool (*compare)(const Polygon_2&, const Polygon_2&);

        bool visible_points(const Point_2&, const Point_2&, const std::unordered_map<Segment_2, bool>&);
        
        //(Polygon to set changes, index of the first vertex of L, index of the last, index of the edge[0], L)
        //remove L from the current potion and put it after edge[0]
        void RelocateEdges(Polygon_2&, int, int, int, int);

        //(Polygon to save the minimized polygon, index of edge[0] to replace, L)
        //returns the difference between the old and new area
        //finds the best L to put in the edge potion
        double ReplaceEdgeWithBest_L(Polygon_2*, const int, const int);

        bool solve(const int, const double, const int);
    
    public:
        //(points, min or max)
        LocalSearch(const std::vector<Point_2>&, const std::string&);

        //(L, threshold, K)
        bool MinimizePolygon(const int, const double, const int K = 0);

        float getPolygonArea();
        Polygon_2 getPolygon();
};

#include "localsearch.cpp"
#endif