#ifndef KDTREE_H
#define KDTREE_H
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <vector>

template<class Kernel> class kdTree
{
    typedef CGAL::Point_2<Kernel>  Point_2;
    typedef std::vector<Point_2> PointsVector;
    typedef typename PointsVector::iterator PointsIterator;

    private:
        class Node;

        Node* root;

        PointsIterator find_median(PointsIterator, PointsIterator, const int);
        Node* insert(PointsIterator, PointsIterator, int depth = 0);
        void points_inside_bounds(PointsVector&, Node* ,const int, const int, const int, const int, int depth = 0);
    
    public:
        kdTree();
        kdTree(PointsVector);
        ~kdTree();
        void find_points_inside_bounds(PointsVector&, const int, const int, const int, const int);//upper x , y lower x , y (in this order)
};

#include "kd_tree.cpp"

#endif